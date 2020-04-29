import numpy as np
from sergio import sergio
np.random.seed(1)
sim = sergio(number_genes=100, number_bins = 1, number_sc = 3000, 
             noise_params = 1, decays=20, sampling_state=1,
             noise_type='dpd')
sim.build_graph('targetFile.py', 'regFile.py', shared_coop_state=10)
np.random.seed(1)
sim.simulate()
np.random.seed(1)
expr = sim.getExpressions()
expr = np.concatenate(expr, axis = 1)
np.random.seed(1)
count_matrix = sim.convert_to_UMIcounts(expr)
np.savetxt('simulationOutput.csv',count_matrix, delimiter=',')