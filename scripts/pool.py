import numpy as np
from sys import argv

raw = np.loadtxt(argv[1])
pool = 100
low = int(np.floor(np.min(raw[:,0])/pool + 0.5))
high = int(np.ceil(np.max(raw[:,0])/pool + 0.5))

result = []
for i in range(low,high):
    result.append( np.average(raw[(raw[:,0] >= i*pool-pool/2) &
                                  (raw[:,0] <= i*pool+pool/2), :], axis=0))

np.savetxt(argv[2], np.column_stack((result, result[::-1])))
