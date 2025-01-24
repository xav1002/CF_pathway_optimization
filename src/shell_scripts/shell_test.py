import sys
sys.path.insert(0, '..\\Classes')
sys.path.insert(0, '..\\Classes\\WholeCellConsortiumModel')
sys.path.insert(0, '..\\Classes\\MetabolicPathwayNetworkGraph')
sys.path.insert(0, '..\\Classes\\MicorbialConsortiumKineticModel')
sys.path.insert(0, '..\\Classes\\Components')
sys.path.insert(0, '..\\..\\Lib\\site-packages')

from WholeCellConsortiumModel import WholeCellConsortiumModel

wccm = WholeCellConsortiumModel()
wccm.generate_whole_network('test')
test_graph_2 = wccm.seek_optimal_network('test',['C00048'],['C00469','C00011'],2,2)
wccm.visualize_graph('test')