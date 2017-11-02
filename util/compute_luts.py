import os
import numpy as np

"""compute lut_in and lut_out for labelconvert from reference_table """

PRD = os.environ['PRD']
PARCEL = os.environ['PARCEL']

lut_in_names = np.loadtxt(open(os.path.join('share', 'reference_table_' + PARCEL + ".csv"), "r"), delimiter=",", skiprows=1, usecols=(1,), dtype='str')
lut_in_vals = np.loadtxt(open(os.path.join('share', 'reference_table_' + PARCEL + ".csv"), "r"), delimiter=",", skiprows=1, usecols=(0, 2, 3, 4, 8), dtype='int')

f = open(os.path.join(PRD, 'connectivity/lut_in.txt'), 'w')
for i, row in enumerate(lut_in_vals):
    f.write(str(row[0]) + ' ')
    f.write(lut_in_names[i]+' ')
    for j in range(1,5):
        f.write(str(row[j]) + ' ')
    f.write('\n')
f.close()



lut_out_names = np.loadtxt(open(os.path.join('share', 'reference_table_' + PARCEL + ".csv"), "r"), delimiter=",", skiprows=1, usecols=(1,), dtype='str')
lut_out_vals = np.loadtxt(open(os.path.join('share', 'reference_table_' + PARCEL + ".csv"), "r"), delimiter=",", skiprows=1, usecols=(5, 2, 3, 4, 8), dtype='int')

f = open(os.path.join(PRD, 'connectivity/lut_out.txt'), 'w')
for i, row in enumerate(lut_out_vals):
    f.write(str(row[0]) + ' ')
    f.write(lut_out_names[i]+' ')
    for j in range(1,5):
        f.write(str(row[j]) + ' ')
    f.write('\n')
f.close()

#lut_out = np.loadtxt(open(os.path.join('share', 'reference_table_' + PARCEL + ".csv"), "r"), delimiter=",", skiprows=1, usecols=(0, 5), dtype='int')
#np.savetxt(os.path.join(PRD, 'connectivity', 'lut_out.txt'), lut_out, fmt='%d %d')
