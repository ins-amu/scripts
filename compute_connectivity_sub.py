import numpy as np
import os
PRD = os.environ['PRD']
SUBJ_ID = os.environ['SUBJ_ID']
act = os.environ['act']
curr_K = np.int(os.environ['curr_K'])
number_tracks = np.int(os.environ['number_tracks'])

# weights, tract_lengths = np.zeros((70*curr_K, 70*curr_K)), np.zeros((70*curr_K, 70*curr_K))
# for i in range(number_tracks):
#     fw = 'weights_' + str(curr_K) + '_' +  str(i+1) + '.csv'
#     ft = 'tract_lengths_' + str(curr_K) + '_' + str(i+1) + '.csv'
#     weights[:, :] += np.loadtxt(os.path.join(PRD, 'connectivity', fw))
#     tract_lengths[:, :] += np.loadtxt(os.path.join(PRD, 'connectivity', ft))
# 
# save connectivity and tract length matrices
weights = np.loadtxt(os.path.join(PRD, 'connectivity', 'weights_' + str(curr_K) + '.csv'))
tract_lengths = np.loadtxt(os.path.join(PRD, 'connectivity', 'tract_lengths_' + str(curr_K) + '.csv'))
weights = weights + weights.transpose() - np.diag(np.diag(weights)) # to avoid twice diag
tract_lengths = tract_lengths + tract_lengths.transpose() # because diagonal nul 
if act=='yes':
    weights = np.loadtxt(os.path.join(PRD, 'connectivity', 'weights_act_' + str(curr_K) + '.csv'))
    tract_lengths = np.loadtxt(os.path.join(PRD, 'connectivity', 'tract_lengths_act_' + str(curr_K) + '.csv'))
    weights = weights + weights.transpose() - np.diag(np.diag(weights)) # to avoid twice diag
    tract_lengths = tract_lengths + tract_lengths.transpose() # because diagonal nul 
    np.savetxt(os.path.join(PRD, SUBJ_ID, 'connectivity_'+str(curr_K), 'weights_act.txt'), weights, fmt='%d')
    np.savetxt(os.path.join(PRD, SUBJ_ID, 'connectivity_'+str(curr_K), 'tract_lengths_act.txt'), tract_lengths, fmt='%.3f')
else:
    weights = np.loadtxt(os.path.join(PRD, 'connectivity', 'weights_' + str(curr_K) + '.csv'))
    tract_lengths = np.loadtxt(os.path.join(PRD, 'connectivity', 'tract_lengths_' + str(curr_K) + '.csv'))
    weights = weights + weights.transpose() - np.diag(np.diag(weights)) # to avoid twice diag
    tract_lengths = tract_lengths + tract_lengths.transpose() # because diagonal nul 
    np.savetxt(os.path.join(PRD, SUBJ_ID, 'connectivity_'+str(curr_K), 'weights.txt'), weights, fmt='%d')
    np.savetxt(os.path.join(PRD, SUBJ_ID, 'connectivity_'+str(curr_K), 'tract_lengths.txt'), tract_lengths, fmt='%.3f')


