import helperfunctions as _hf
from proteingraphanalyser import ProteinGraphAnalyser

class ConservedGraph(ProteinGraphAnalyser):
    def __init__(self, pdb_root_folder, target_folder=''):
        ProteinGraphAnalyser.__init__(self, pdb_root_folder, target_folder='')