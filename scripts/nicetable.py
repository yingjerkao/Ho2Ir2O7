import numpy as np
from sys import argv

raw = np.loadtxt(argv[1])

dircode = argv[2]
if dircode.startswith('100'):
    DIR = np.array([1,0,0])
elif dircode.startswith('111'):
    DIR = np.array([1,1,1])/3**0.5
else:
    DIR = None

SMP = np.array([0,1,0,1,0]*2)
DMP = np.array([1,0,0,0,1]*2)
MPC = np.array([-2,-1,0,1,2]+[0]*5)

Mfactor = 1/3**0.5/1728
Dfactor = 1.0/1728

if DIR is not None:
    mag = Mfactor * raw[:,1:4]@DIR
    cols = 8
else:
    mag = Mfactor * raw[:,1:4]
    cols = 12
single = Dfactor * raw[:,4:14]@SMP
double = Dfactor * raw[:,4:14]@DMP
charge = Dfactor * raw[:,4:14]@MPC

todos = np.column_stack((raw[:,0]/1e4,
                         mag, -np.flipud(mag),
                         single, np.flip(single),
                         double, np.flip(double),
                         -charge, np.flip(charge)))
np.savetxt(argv[3], todos, fmt = "%+5.2f" + " %10.7f"*cols)
