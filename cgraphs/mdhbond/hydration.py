#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
#    Author: Malte Siemers, Freie Universität Berlin
#
#    If you use this software or anything it produces for work to be published,
#    please cite:
#
#    Malte Siemers, Michalis Lazaratos, Konstantina Karathanou,
#    Federico Guerra, Leonid Brown, and Ana-Nicoleta Bondar.
#    HWire: A graph-based algorithm to analyze dynamic H-bond networks
#    in membrane proteins, Journal of Chemical Theory and Computation, 2019.


from . import helpfunctions as _hf
from .basic import BasicFunctionality
import numpy as _np
from scipy import spatial as _sp
from scipy.optimize import curve_fit as _curve_fit
from scipy.special import gamma as _gamma
from matplotlib.ticker import MaxNLocator
#import matplotlib
#matplotlib.use('TKAgg', warn=False)
import matplotlib.pyplot as _plt

_np.int = int

class HydrationAnalysis(BasicFunctionality):

    def __init__(self, selection=None, structure=None, trajectories=None,
                 start=None, stop=None, step=1, restore_filename=None):

        super(HydrationAnalysis, self).__init__(selection=selection, structure=structure, trajectories=trajectories,
             start=start, stop=stop, step=step, restore_filename=restore_filename)


    def set_presence_in_hull(self):
        result = {}
        frame_count = 0
        frames = self.nb_frames

        for ts in self._universe.trajectory[self._trajectory_slice]:
            water_coordinates = self._water.positions
            select_coordinates = self._mda_selection.positions

            hull = _sp.Delaunay(select_coordinates)
            local_water_index = (hull.find_simplex(water_coordinates) != -1).nonzero()[0]

            frame_res = [self._water_ids[i] for i in local_water_index]

            for water in frame_res:
                try:
                    result[water][frame_count] = True
                except:
                    result[water] = _np.zeros(frames, dtype=bool)
                    result[water][frame_count] = True
            frame_count+=1
        self.initial_results = result
        self.filtered_results = result

    def set_presence_around(self, around_distance, per_residue=False):
        result = {}
        frame_count = 0
        frames = self.nb_frames
        _all_ids =  _hf.MDA_info_list(self._mda_selection) + self._water_ids
        _first_water_id = len(self._mda_selection)

        for ts in self._universe.trajectory[self._trajectory_slice]:

            water_coordinates = self._water.positions
            selection_coordinates = self._mda_selection.positions

            selection_tree = _sp.cKDTree(selection_coordinates)
            water_tree = _sp.cKDTree(water_coordinates, leafsize=32)
            local_water_index = []
            [local_water_index.extend(l) for l in selection_tree.query_ball_tree(water_tree, float(around_distance))]
            local_water_index = _np.unique(local_water_index)

            if per_residue:
                local_water_coordinates = water_coordinates[local_water_index]
                local_water_tree = _sp.cKDTree(local_water_coordinates)

                local_pairs = [(i, local_water_index[j]+_first_water_id) for i, bla in enumerate(selection_tree.query_ball_tree(local_water_tree, around_distance)) for j in bla]
                check = [_all_ids[j] + ':' + _all_ids[i] in result for (i,j) in local_pairs]
                frame_res = [_all_ids[j] + ':' + _all_ids[i] if check[ii] else _all_ids[i] + ':' + _all_ids[j] for ii, (i, j) in enumerate(local_pairs)]
            else:
                frame_res = [self._water_ids[i] for i in local_water_index]

            for water in frame_res:
                try:
                    result[water][frame_count] = True
                except:
                    result[water] = _np.zeros(frames, dtype=bool)
                    result[water][frame_count] = True
            frame_count+=1
        self.initial_results = result
        self.filtered_results = result

    def filter_occupancy(self, min_occupancy, use_filtered=True):
        if use_filtered: results = self.filtered_results
        else: results = self.initial_results
        if len(results) == 0: raise AssertionError('nothing to filter!')
        filtered_result = {key:results[key] for key in results if _np.mean(results[key])>min_occupancy}
        if self.applied_filters['occupancy'] != None: self.applied_filters['occupancy']=max(self.applied_filters['occupancy'], min_occupancy)
        else: self.applied_filters['occupancy'] = min_occupancy
        self.filtered_results = filtered_result

    def compute_mean_square_displacement(self, nb_blocks=1, use_filtered=True):
        if use_filtered: results = self.filtered_results
        else: results = self.initial_results
        water_index = _np.in1d(self._water_ids, [key for key in results])
        index_trans = _np.nonzero(water_index)[0]
        water_involved = self._water[water_index]
        water_coords = _np.empty((self.nb_frames, water_index.sum(), 3))
        for i,t in enumerate(self._universe.trajectory[self._trajectory_slice]):
            water_coords[i] = water_involved.positions
        all_msd = []
        block_length = int(self.nb_frames/nb_blocks)
        for run in range(nb_blocks):
            s = slice(run*block_length,(run+1)*block_length,1)
            msds=_np.zeros((len(water_involved), block_length))
            for water in range(len(water_involved)):
                ts = results[self._water_ids[index_trans[water]]][s]
                c = water_coords[s][:,water,:]
                msds[water][ts] = _hf.msd_fft(c)[ts]
            msds = msds.mean(0)
            all_msd.append(msds)
        all_msd = _np.array(all_msd)
        self.msd = all_msd.mean(0), all_msd.std(0)/_np.sqrt(nb_blocks)
        return self.msd

    def compute_diffusion_coefficient(self, fit_begin=0.1, fit_end=0.9, nb_blocks=1, use_filtered=True):
        msd, err = self.msd
        block_length = len(msd)
        cut_begin = int(block_length*fit_begin)
        cut_end = int(block_length*fit_end)
        def func(x, m, c):
            return m*x + c
        popt, pcov = _curve_fit(func, _np.arange(block_length)[cut_begin:cut_end], msd[cut_begin:cut_end])
        return popt[0]/6.

    def draw_msd(self, scatter_size=0.5, filename=None, return_figure=False):
        msd, err = self.msd
        fig, ax = _plt.subplots()
        block_length = len(msd)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        _plt.xlabel('Frames' , fontsize = 16)
        _plt.ylabel('Mean Square Displacement' , fontsize = 16)
        _plt.scatter(_np.arange(block_length), msd, s=scatter_size)
        if return_figure:
            _plt.close()
            return fig
        self._save_or_draw(filename, return_figure=return_figure)

    def _compute_mean_square_displacement_slow(self, nb_blocks=1, use_filtered=True):
        if use_filtered: results = self.filtered_results
        else: results = self.initial_results

        water_index = _np.in1d(self._water_ids, [key for key in results])
        water_involved = self._water[water_index]
        res_vals = _np.array([results[key] for key in results]).T
        water_coords = _np.empty((self.nb_frames, water_index.sum(), 3), dtype=_np.float16)
        for i,t in enumerate(self._universe.trajectory[self._trajectory_slice]):
            water_coords[i] = water_involved.positions
        all_msd = []
        block_length = int(self.nb_frames/nb_blocks)
        for run in range(nb_blocks):
            s = slice(run*block_length,(run+1)*block_length,1)
            msd = _np.zeros(int(block_length / 2) -1)
            for t0 in range(int(block_length / 2)):
                self._universe.trajectory[s][t0]
                c0 = water_involved.positions
                for dt in range(1, int(block_length / 2)):
                    local_index = res_vals[t0] & res_vals[t0+dt]
                    if local_index.sum() == 0:
                        msd[dt-1] = 0.
                        continue
                    self._universe.trajectory[s][t0+dt]
                    d = c0[local_index] - water_involved.positions[local_index]
                    msd[dt-1] += (d**2).sum(1).mean()
            msd /= t0+1
            all_msd.append(msd)
        all_msd = _np.array(all_msd)
        return all_msd.mean(0), all_msd.std(0)/_np.sqrt(nb_blocks)


    def compute_mean_residence_time(self, nb_blocks=1, filter_artifacts=False, per_residue=False, return_average_time=False, use_filtered=True):
        if use_filtered: results = self.filtered_results
        else: results = self.initial_results

        def func(x, tau, lamda):
            return _np.exp(-(x/tau) ** lamda)
        all_residence_time = []
        all_average_time = []
        block_length = int(self.nb_frames/nb_blocks)
        xdata = _np.arange(1,block_length+1)
        for run in range(nb_blocks):
            s = slice(run*block_length,(run+1)*block_length,1)
            intervals = {key:_hf.intervals_binary(results[key][s]) for key in results}
            residence_time = _np.zeros(block_length)
            average_time = _np.array([0.,0.])

            if per_residue:
                residence_time = {}
                del_res = []
                average_time = {}
                for pair in results:
                    segn, resn, resi, _, _, _ = _hf.deconst_key(pair, True)
                    residue_key = segn+'-'+resn+'-'+str(resi)
                    work_intervals = intervals[pair]
                    interval_lengths = _np.diff(work_intervals, 1).astype(_np.int).flatten()
                    if filter_artifacts: interval_lengths = interval_lengths[(interval_lengths<block_length)]
                    try: average_time[residue_key] += _np.array([interval_lengths.sum(),interval_lengths.size])
                    except KeyError:
                        average_time[residue_key] = _np.array([interval_lengths.sum(),interval_lengths.size])
                        residence_time[residue_key] = _np.zeros(block_length)
                    for l in interval_lengths:
                        residence_time[residue_key][:l] += _np.arange(1, l+1)[::-1]

                for residue_key in residence_time:
                    residence_time[residue_key] /= _np.arange(1,block_length+1)[::-1]
                    residence_time[residue_key] /= residence_time[residue_key][0]
                    if _np.isnan(residence_time[residue_key]).any():
                        del_res.append(residue_key)
                    average_time[residue_key] = average_time[residue_key][0]/average_time[residue_key][1]
                    try:
                        (tau, lamda), pcov = _curve_fit(func, xdata, residence_time[residue_key], p0=(10.0, 0.5))
                        residence_time[residue_key] = tau/lamda * _gamma(1/lamda)
                    except:
                        if average_time[residue_key] <= 2.0:
                            residence_time[residue_key] = 1.0
                        else:
                            try:
                                (tau, lamda), pcov = _curve_fit(func, xdata, residence_time[residue_key], p0=(block_length, 0.5))
                                residence_time[residue_key] = tau/lamda * _gamma(1/lamda)
                            except:
                                residence_time[residue_key] = block_length
                for key in del_res:
                    del residence_time[key]

            else:
                for water in results:
                    work_intervals = intervals[water]
                    interval_lengths = _np.diff(work_intervals, 1).astype(_np.int).flatten()
                    if filter_artifacts: interval_lengths = interval_lengths[(interval_lengths<block_length)]
                    average_time += _np.array([interval_lengths.sum(),interval_lengths.size])
                    for l in interval_lengths:
                        residence_time[:l] += _np.arange(1, l+1)[::-1]
                residence_time /= _np.arange(1,block_length+1)[::-1]
                residence_time /= residence_time[0]
                average_time = average_time[0]/average_time[1]
                try:
                    (tau, lamda), pcov = _curve_fit(func, xdata, residence_time, p0=(10.0, 0.5))
                    residence_time = tau/lamda * _gamma(1/lamda)
                except:
                    if average_time <= 2.0:
                        residence_time = 1.0
                    else:
                        residence_time = block_length
            all_residence_time.append(residence_time)
            all_average_time.append(average_time)
        if per_residue:
            for key in residence_time:
                temp_residence_time = []
                temp_average_time = []
                for residence_dict in all_residence_time:
                    try:temp_residence_time.append(residence_dict[key])
                    except:pass
                for average_dict in all_average_time:
                    try:temp_average_time.append(average_time[key])
                    except:pass
                residence_time[key] = (_np.mean(temp_residence_time), _np.std(temp_residence_time)/_np.sqrt(nb_blocks))
                average_time[key] = (_np.mean(temp_average_time), _np.std(temp_average_time)/_np.sqrt(nb_blocks))
        else:
            residence_time = (_np.mean(all_residence_time), _np.std(all_residence_time)/_np.sqrt(nb_blocks))
            average_time = (_np.mean(all_average_time), _np.std(all_average_time)/_np.sqrt(nb_blocks))

        if return_average_time: return residence_time, average_time
        else: return residence_time

    def draw_mean_residence_time_per_residue(self, resids_to_plot, nb_blocks=1, filter_artifacts_above=None, filename=None, return_figure=False):
        residence_time =  self.compute_mean_residence_time(nb_blocks=nb_blocks, per_residue=True)
        segnames, resnames, resids, residence, errors = [], [], [], [], []
        for key in residence_time:
            segname, resname, resid = key.split('-')
            segnames.append(segname)
            resnames.append(resname)
            resids.append(int(resid))
            residence.append(residence_time[key][0])
            errors.append(residence_time[key][1])
        index = _np.argsort(resids)
        segnames, resnames, resids, residence_time, errors = _np.array(segnames)[index], _np.array(resnames)[index], _np.array(resids)[index], _np.array(residence)[index], _np.array(errors)[index]
        if filter_artifacts_above != None:
            f = residence_time>filter_artifacts_above
            residence_time[f]=0
            errors[f]=0
        index = _np.in1d(resids, resids_to_plot)
        segnames, resnames, resids, residence_time, errors = _np.array(segnames)[index], _np.array(resnames)[index], _np.array(resids)[index], _np.array(residence_time)[index], _np.array(errors)[index]
        u_segnames = _np.unique(segnames)
        per_segment = {seg:[] for seg in u_segnames}
        per_segment_error = {seg:[] for seg in u_segnames}
        for i in range(len(residence_time)): per_segment[segnames[i]].append(residence_time[i])
        for i in range(len(residence_time)): per_segment_error[segnames[i]].append(errors[i])
        labels = _np.array([resnames[i]+str(resids[i]) for i in range(len(resnames))])
        _, idx = _np.unique(labels, return_index=True)
        labels = labels[_np.sort(idx)]
        fig, ax = _plt.subplots()
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        _plt.xlabel('Residue' , fontsize = 16)
        _plt.ylabel('Mean Residence Time [Frames]' , fontsize = 16)
        if u_segnames.size == 1: _plt.bar(_np.arange(len(residence_time)), residence_time, yerr=errors)
        else:
            width = 0.8/u_segnames.size
            centers = _np.linspace(-1,1,u_segnames.size)*(width/2)
            for i in range(len(u_segnames)):
                per_segment[u_segnames[i]] = _np.array(per_segment[u_segnames[i]])
                per_segment[u_segnames[i]][per_segment[u_segnames[i]]>self.nb_frames]=0
                _plt.bar(_np.arange(len(per_segment[u_segnames[i]])) + centers[i], per_segment[u_segnames[i]], yerr=per_segment_error[u_segnames[i]], width=width, label=u_segnames[i])
            _plt.legend()
        _plt.xticks(_np.arange(-1, len(labels), 1.0), ['']+list(labels), rotation=45)
        if return_figure:
            _plt.close()
            return fig
        self._save_or_draw(filename, return_figure=return_figure)


    def draw_mean_residence_time_histogram(self, nb_bins=50, filter_lower_than=None, filename=None, return_figure=False):
        residence_time =  self.compute_mean_residence_time(per_residue=True)
        residence_time = _np.array([residence_time[a] for a in residence_time])
        if filter_lower_than != None: residence_time = residence_time[residence_time<filter_lower_than]
        fig, ax = _plt.subplots()
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        _plt.xlabel('Mean Residence Time [Frames]' , fontsize = 16)
        _plt.ylabel('Count' , fontsize = 16)
        _plt.hist(residence_time, nb_bins)
        if return_figure:
            _plt.close()
            return fig
        self._save_or_draw(filename, return_figure=return_figure)

    def compute_water_count(self, use_filtered=True):
        if use_filtered: water_presence = self.filtered_results
        else: water_presence = self.initial_results
        presence = _np.array([val for val in water_presence.values()]).T
        return presence.sum(1)

    def draw_water_count(self, use_filtered=True, filename=None, return_figure=False):
        water_count = self.compute_water_count(use_filtered)
        fig, ax = _plt.subplots()
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        _plt.xlabel('Frame' , fontsize = 16)
        _plt.ylabel('Count' , fontsize = 16)
        _plt.plot(water_count)
        if return_figure:
            _plt.close()
            return fig
        self._save_or_draw(filename, return_figure=return_figure)

    def compute_influx_count(self, return_mean_and_std = False, use_filtered=True):
        if use_filtered: water_presence = self.filtered_results
        else: water_presence = self.initial_results
        presence = _np.array([val for val in water_presence.values()]).T
        results = _np.zeros((2,self.nb_frames))
        #base = presence[0].sum()
        for frame in range(self.nb_frames-1):
            results[0][frame] = (_np.logical_not(presence[frame]) & presence[frame+1]).sum()
            results[1][frame] = (presence[frame] & _np.logical_not(presence[frame+1])).sum()
        if return_mean_and_std: return results, results.mean(1), results.std(1)
        else: return results

    def draw_influx_count(self, interval=1, use_filtered=True, filename=None, return_figure=False):
        results = self.compute_influx_count(use_filtered=use_filtered)
        rest = (results[0].size%interval) * -1
        if rest != 0: results = results.T[:rest].T
        influx = results[0].reshape((-1,interval)).mean(1)
        outflux = results[1].reshape((-1,interval)).mean(1)*-1
        xdata = _np.arange(len(influx))*interval
        fig, ax = _plt.subplots()
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        _plt.bar(xdata,influx, label='influx', width=interval)
        _plt.bar(xdata,outflux, label='efflux', width=interval)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        _plt.legend()
        if return_figure:
            _plt.close()
            return fig
        self._save_or_draw(filename, return_figure=return_figure)

