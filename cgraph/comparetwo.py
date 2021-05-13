from . import helperfunctions as _hf
import shutil
import numpy as np
from .proteingraphanalyser import ProteinGraphAnalyser
from .conservedgraph import ConservedGraph
import matplotlib.pyplot as plt


class CompareTwo(ProteinGraphAnalyser):
    def __init__(self, type_option, pdb1=None, pdb2=None, psf1=None, psf2=None, dcd1=None, dcd2=None, target_folder=''):
        self.type_option = type_option
        if self.type_option == 'pdb':
            self.pdb1_code = _hf.retrieve_pdb_code(pdb1, '.pdb')
            self.pdb2_code = _hf.retrieve_pdb_code(pdb2, '.pdb')
            self.compare_folder = _hf.create_directory(target_folder+'/workfolder/compare_'+self.pdb1_code+'_'+self.pdb2_code)+'/'
            shutil.copy(pdb1, self.compare_folder+self.pdb1_code+'.pdb')
            shutil.copy(pdb2, self.compare_folder+self.pdb2_code+'.pdb')

            ProteinGraphAnalyser.__init__(self, pdb_root_folder=self.compare_folder, target_folder=target_folder, reference_pdb=pdb1)
            self.logger.info('COMPARE STRUCTURES '+ self.pdb1_code + ' WITH ' + self.pdb2_code)
            ProteinGraphAnalyser.align_structures(self)
        # elif self.type_option == 'dcd' and len(dcd_files) and len(psf_file):
        #     pass

    def plot_graph_comparison(self, color1='blue', color2='green', label_nodes=True, label_edges=True, xlabel='PCA projected xy plane', ylabel='Z coordinates ($\AA$)',occupancy=None):
        if len(self.graph_coord_objects.items()) != 2: self.logger.warning('There are '+str(len(self.graph_coord_objects.items()))+' structures selected. Graph comparison is possible for exactly two structures.')
        else:
            ConservedGraph.get_conserved_graph(self)
            pos1 = self._get_node_positions(self.graph_coord_objects[self.pdb1_code], pca=False)
            pos2 = self._get_node_positions(self.graph_coord_objects[self.pdb2_code], pca=False)

            all_pos = {}
            for i, pos in enumerate([pos1, pos2]):
                for key, value in pos.items():
                    if key.split('-')[1]=='HOH': key=str(i+1)+'-HOH-'+key.split('-')[2]
                    if key not in all_pos.keys(): all_pos.update({ key:value })

            node_pca_pos = _hf.calculate_pca_positions(all_pos)

            plot_name = 'H-bond' if self.graph_type == 'hbond' else 'Water wire'
            fig, ax = _hf.create_plot(title=plot_name+' graph comparison of structure '+self.pdb1_code+' with '+self.pdb2_code, xlabel=xlabel, ylabel=ylabel)
            graph1 = self.graph_coord_objects[self.pdb1_code]['graph']
            graph2 = self.graph_coord_objects[self.pdb2_code]['graph']

            for e in graph1.edges:
                e0 = e[0] if e[0].split('-')[1] != 'HOH' else '1-HOH-'+e[0].split('-')[2]
                e1 = e[1] if e[1].split('-')[1] != 'HOH' else '1-HOH-'+e[1].split('-')[2]
                color = 'gray' if _hf.is_conserved_edge(np.array([[e2[0], e2[1]] for e2 in graph2.edges]), e0, e1) else color1

                edge_line = [node_pca_pos[e0], node_pca_pos[e1]]
                x=[edge_line[0][0], edge_line[1][0]]
                y=[edge_line[0][1], edge_line[1][1]]
                ax.plot(x, y, color=color, marker='o', linewidth=2, markersize=15, markerfacecolor=color, markeredgecolor=color)

            for e in graph2.edges:
                e0 = e[0] if e[0].split('-')[1] != 'HOH' else '2-HOH-'+e[0].split('-')[2]
                e1 = e[1] if e[1].split('-')[1] != 'HOH' else '2-HOH-'+e[1].split('-')[2]

                if not _hf.is_conserved_edge(np.array([[e2[0], e2[1]] for e2 in graph1.edges]), e0, e1):
                    edge_line = [node_pca_pos[e0], node_pca_pos[e1]]
                    x=[edge_line[0][0], edge_line[1][0]]
                    y=[edge_line[0][1], edge_line[1][1]]
                    ax.plot(x, y, color=color2, marker='o', linewidth=2, markersize=15, markerfacecolor=color2, markeredgecolor=color2)

            for n in graph1.nodes:
                n = n if n.split('-')[1] != 'HOH' else '1-HOH-'+n.split('-')[2]
                color = 'gray' if n in graph2.nodes else color1
                if n.split('-')[1] == 'HOH':
                    ax.scatter(node_pca_pos[n][0], node_pca_pos[n][1],color='#db5c5c', s=150, zorder=5, edgecolors=color)
                else: ax.scatter(node_pca_pos[n][0], node_pca_pos[n][1], s=200, color=color, zorder=5)

            for n in graph2.nodes:
                n = n if n.split('-')[1] != 'HOH' else '2-HOH-'+n.split('-')[2]
                if n not in graph1.nodes:
                    if n.split('-')[1] == 'HOH':
                        ax.scatter(node_pca_pos[n][0], node_pca_pos[n][1], color='#db5c5c', s=150, zorder=5, edgecolors=color2)
                    else: ax.scatter(node_pca_pos[n][0], node_pca_pos[n][1], s=200, color=color2, zorder=5)

            if label_nodes:
                for n in all_pos.keys():
                    values = node_pca_pos[n]
                    # if n.split('-')[1] == 'HOH': ax.annotate('W'+str(int(n.split('-')[2])), (values[0]+0.2, values[1]-0.25), fontsize=12)
                    if n.split('-')[1] == 'HOH': ax.annotate('W', (values[0]+0.3, values[1]-0.25), fontsize=12)
                    else: ax.annotate(str(n.split('-')[0])+'-'+str(_hf.amino_d[n.split('-')[1]])+str(int(n.split('-')[2])), (values[0]+0.2, values[1]-0.25), fontsize=12)
            ax.text(0.97, 0.99, self.pdb1_code, color=color1, fontsize=20, transform=ax.transAxes, ha='center', va='center')
            ax.text(0.97, 0.97, self.pdb2_code, color=color2, fontsize=20, transform=ax.transAxes, ha='center', va='center')

            plt.tight_layout()
            is_label = '_labeled' if label_nodes else ''
            if self.graph_type == 'hbond':
                plt.savefig(self.compare_folder+'compare_H-bond_graph'+self.pdb1_code+'_with_'+self.pdb2_code+is_label+'.png')
                plt.savefig(self.compare_folder+'compare_H-bond_graph'+self.pdb1_code+'_with_'+self.pdb2_code+is_label+'.eps', format='eps')
            elif self.graph_type == 'water_wire':
                waters = '_max_'+str(self.max_water)+'_water_bridges' if self.max_water > 0 else ''
                occ = '_min_occupancy_'+str(occupancy) if occupancy  else ''
                plt.savefig(self.compare_folder+'compare'+waters+occ+'_graph'+self.pdb1_code+'_with_'+self.pdb2_code+is_label+'.png')
                plt.savefig(self.compare_folder+'compare'+waters+occ+'_graph'+self.pdb1_code+'_with_'+self.pdb2_code+is_label+'.eps', format='eps')
            plt.close()

