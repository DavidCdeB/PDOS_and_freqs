# 
# QHA -son- program: This program is going to be called by the QHA -Master- program
# David Carrasco de Busturia, 13 October 2017 
# Please read the documentation and istructions on: https://github.com/DavidCdeB/QHA
# This program is under the GNU General Public License v3.0. 
# davidcarrascobustur@gmail.com
# d.carrasco-de-busturia@imperial.ac.uk

import re
import os
import glob
from itertools import islice
import numpy as np
import subprocess
import itertools
import sys
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


# Extracting the number of k-points ( name of variable = "N_k" ):
path='./'
template = os.path.join(path, '125.845303_RESTART.out')

ALL_FREQ = []

for fname in glob.glob(template):
  print fname
  f = open(fname, 'r')
# real_part = False

# def read_real_parts(fname):
  for line in f:

        if 'MODES         EIGV          FREQUENCIES     IRREP' in line:

            print line
            print f.next()

            while True:
               target = f.next()
               aux = target.split()
               if not aux:
                  break

               first_No = aux[0]
               second_No = aux[1]
               freq = aux[3]
               print 'freq = ', freq
               print ' first_No original = ', first_No

               first_No = first_No.translate(None, '-')  # remove the '-'

               print ' first_No = ', first_No
               print 'second_No = ', second_No

               factor_freq = abs(int(second_No) - int(first_No)) + 1
               print 'factor_freq = ', factor_freq

               freqs = [freq] * factor_freq
               print 'freqs = ', freqs

               ALL_FREQ.append(freqs)


print 'ALL_FREQ = ', ALL_FREQ
ALL_FREQ = [item for sublist in ALL_FREQ for item in sublist]
print 'ALL_FREQ = ', ALL_FREQ
print 'len(ALL_FREQ) = ', len(ALL_FREQ)

thing = '0.0000'
while thing in ALL_FREQ: ALL_FREQ.remove(thing)
print 'ALL_FREQ = ', ALL_FREQ
print 'len(ALL_FREQ) = ', len(ALL_FREQ)


print ' len(ALL_FREQ) = ', len(ALL_FREQ)
ALL_FREQ = sorted(ALL_FREQ, key=float)  

output_array_2 = np.vstack((ALL_FREQ))#.T
np.savetxt('All_freq.dat', output_array_2, header="FREQS (CM^-1)", fmt="%s")
sys.exit()
os.system('''grep -v "^#"  *PHONDOS | grep -v "^@ " |awk '{ print $1 }'  > freqs_from_PHONDOS.dat''')
os.system('''sort -k1 -n freqs_from_PHONDOS.dat > freqs_from_PHONDOS_sorted.dat''')

freqs_from_PHONDOS_sorted = np.loadtxt('freqs_from_PHONDOS_sorted.dat').T
print 'freqs_from_PHONDOS_sorted  = ', freqs_from_PHONDOS_sorted
output_array_2 = np.vstack((freqs_from_PHONDOS_sorted))
np.savetxt('freqs_from_PHONDOS_sorted_decimals.dat', output_array_2, header="FREQS (CM^-1)", fmt="%.4f")

print '''

All_freq.dat -> All frecuencies from the output file
freqs_from_PHONDOS_sorted_decimals.dat -> Frecuencies (energy) from the PHONDOS file

'''

print 'len(ALL_FREQ) = ', len(ALL_FREQ)
print 'len(freqs_from_PHONDOS_sorted) = ', len(freqs_from_PHONDOS_sorted)

double_check = raw_input("""

Press <ok> to run 
vimdiff  freqs_from_PHONDOS_sorted_decimals.dat  All_freq.dat

""")

os.system('''vimdiff  freqs_from_PHONDOS_sorted_decimals.dat  All_freq.dat''')
