import numpy as np
import os
PRD = os.environ['PRD']
SUBJ_ID = os.environ['SUBJ_ID']
curr_K = np.int(os.environ['curr_K'])
number_tracks = np.int(os.environ['number_tracks'])

# save connectivity and tract length matrices
weights = np.loadtxt(os.path.join(PRD, 'connectivity', 'weights_' + str(curr_K) + '.csv'))
tract_lengths = np.loadtxt(os.path.join(PRD, 'connectivity', 'tract_lengths_' + str(curr_K) + '.csv'))
weights = weights + weights.transpose() - np.diag(np.diag(weights)) # to avoid twice diag
tract_lengths = tract_lengths + tract_lengths.transpose() # because diagonal nul 
np.savetxt(os.path.join(PRD, SUBJ_ID, 'connectivity_'+str(curr_K), 'weights.txt'), weights, fmt='%d')
np.savetxt(os.path.join(PRD, SUBJ_ID, 'connectivity_'+str(curr_K), 'tract_lengths.txt'), tract_lengths, fmt='%.3f')


