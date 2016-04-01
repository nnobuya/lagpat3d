#! /usr/bin/env python

import ConfigParser
import numpy as np

#### read config file

config_fl = ConfigParser.SafeConfigParser()
config_fl.read('./in/config.in')

n_adopt = config_fl.getint('set_calc', 'n_adopt')
select  = config_fl.getint('set_calc', 'select')

print('  - adopts ' + str(n_adopt) + ' particles in each YeS cell')
if select == 1: print('  - re-arrange the order of particles by masse factors')


n_pt = []; f_pt = []; ye = []; et = []; fac = []

for line in open('ej_pt_list.dat'):
    dat = line.split()

    ndat = (len(dat) - 6) /2
    n_tmp0 = []; f_tmp0 =[]
    for i in range(ndat):
        n_tmp0.append(dat[6 + 2*i])
        f_tmp0.append(dat[6 + 2*i + 1])

    ### re-arranges order
    if select == 0:
        n_tmp = n_tmp0
        f_tmp = f_tmp0
    if select == 1:
        n_tmp = []; f_tmp = []
        isort = np.argsort(f_tmp0)[-1::-1]
        for i1 in isort:
            n_tmp.append(n_tmp0[i1])
            f_tmp.append(f_tmp0[i1])

    n_pt.append(n_tmp)
    f_pt.append(n_tmp)
    ye.append(float(dat[2]))
    et.append(float(dat[3]))
    fac.append(float(dat[4]))


total = sum(fac)
fac   = [ r1 /total for r1 in fac ]



out = open('./pt_list.dat', 'w')

out.write('#' + '{0:>13}{1:>10}{2:>10}{3:>20}'.
          format('No', 'Ye', 'S', 'Fac') + '\n')

n = 0; fac_total = 0.0
for i in range(len(n_pt)):

    nout = min(len(n_pt[i]), n_adopt)

    sep = []
    for j in range(nout):
        sep.append(float(f_pt[i][j]))

    total = sum(sep)
    sep   = [ r1 /total for r1 in sep ]

    for j in range(nout):
        n += 1
        out.write('{0:>7d}{1:>7}'.format(n, n_pt[i][j]))
        out.write('{0:10.4f}{1:10.2f}'.format(ye[i],et[i]))
        out.write('{0:20.10e}'.format(fac[i] *sep[j]))
        out.write('\n')

out.close()



exit('terminated')
