import numpy as np
import sys

# save connectivity and tract length matrices
weights = np.loadtxt(sys.argv[1])
tract_lengths = np.loadtxt(sys.argv[2])
weights = weights + weights.transpose() - np.diag(np.diag(weights)) # to avoid twice diag
tract_lengths = tract_lengths + tract_lengths.transpose() # because diagonal nul 
# tck2connectome produces the number of steps in the tracking file
# so we have to multiply by the step size in mm
# which is 0.5
tract_lengths *= 0.5 
np.savetxt(sys.argv[3], weights, fmt='%d')
np.savetxt(sys.argv[4], tract_lengths, fmt='%.3f')


