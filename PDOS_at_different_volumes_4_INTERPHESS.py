#
import matplotlib.pyplot as plt
import sys
import numpy as np
import os

folder_222 = '/home/david/Trabajo/structures/PDOS_at_different_volumes/Calcite_I/pob_TZVP/PBE-D3/Only_Veq/3_times_landau_supercell'

folder_222_INTERPHESS_222 = 'INTERPHESS_222'

files_222 = ["125.845303_RESTART_PDOS.PHONDOS"] 
files_222_INTER = ["125.845303_RESTART_PDOS_INTERPHESS_222.PHONDOS"]

vols_all = ["125.845303"]

vols = vols_all#[0:3]
files_222_0_3 = files_222#[0:3]
files_222_0_3_INTER = files_222_INTER#[0:3]

fig = plt.figure()
for indx_vols, indx_files in zip(range(1, len(vols)+1), range(len(vols))):
        print 'indx_vols init = ' , indx_vols
        print 'indx_files init = ' , indx_files
        ax = fig.add_subplot(len(vols), 1, indx_vols)
        ax.text(.85,.80, '\nPBE-D3, pob-TVPZ\nV$_{eq}$ = ' + vols[indx_vols-1] + ' $\AA^{3}$', fontsize=10, #alpha=1.0,
            bbox=dict(alpha=1.0, facecolor='white', pad=1.0), #, facecolor='red', edgecolor='red'),
            horizontalalignment='center',
            transform=ax.transAxes)
        x_222, y_222 = np.loadtxt(os.path.join(folder_222, files_222_0_3[indx_files]), skiprows = 4).T
        x_222_INTER, y_222_INTER = np.loadtxt(os.path.join(folder_222, folder_222_INTERPHESS_222, files_222_0_3_INTER[indx_files]), skiprows = 4).T
        ax.plot(x_222, y_222, color='lightskyblue', label='SCELPHONO 2 times Landau supercell') # black lines OK
        ax.plot(x_222_INTER, y_222_INTER, color='dodgerblue', label='SCELPHONO 2 times Landau supercell') # black lines OK
        ax.grid()

        if indx_vols == 1 and indx_files == 0:
           print 'enter 1st if loop'
           print 'indx_vols == 1 and indx_files == 0: ', indx_files 
           ax.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=4)
           ax.set_ylabel('PDOS (states/cm$^{-1}$/cell)')
           ax.tick_params(labelbottom='on') 


        if indx_vols == len(vols) and indx_files == len(vols)-1:
           print 'enter 2nd if loop'
           print 'indx_files if indx_vols == len(vols)+1 and indx_files == len(vols): = ' , indx_files
           ax.set_xlabel('Energy (cm$^{-1}$)')

        else:
            ax.tick_params(labelbottom='on')

fig.savefig("PDOS_0-3.pdf", bbox_inches='tight', pad_inches=0.3)#, tight_layout() )#, bbox_inches=bbox)

#plt.show()

# 1st Step: Calculating the cumulated entropy:

nu_127_SCEL_222  = np.loadtxt('/home/david/Trabajo/structures/PDOS_at_different_volumes/Calcite_I/pob_TZVP/PBE-D3/Only_Veq/3_times_landau_supercell/All_freq.dat', skiprows = 1).T  
nu_127_SCEL_222_INTER  = np.loadtxt('/home/david/Trabajo/structures/PDOS_at_different_volumes/Calcite_I/pob_TZVP/PBE-D3/Only_Veq/3_times_landau_supercell/INTERPHESS_222/All_freq.dat', skiprows = 1).T  


################### CONSTANTS   ###############
global KB, h, c, T

# KB = boltmann cte, KB = 1.38064852(79)x10-23 J/K
KB = 1.38064852E-23

# h = plank constant, h = 6.626070040(81)x10-34 J s
h = 6.626070040E-34

# T = temperature, T = 298.15 K
T = 298.15

# c = speed of light, c = 2.99792458E8 m/s
c = 2.99792458E+8

def S(nu):
 S = -KB * np.log(1 - np.exp(-h * nu * 1E+2 * c  / (KB*T)))     + ( (h/T) * ( nu * 1E+2 * c * ( (np.exp(h *  nu  * 1E+2 * c  / (KB*T)) - 1)**(-1) )    )  )

 # Conversion: S above this line is in mHartree/K:
 return S  *((1/4.3597482)*1E+18 * 1E+3)

S_127_SCEL_222 = S(nu_127_SCEL_222)
S_127_SCEL_222_INTER = S(nu_127_SCEL_222_INTER)

for i in xrange(1,len(S_127_SCEL_222)):
    S_127_SCEL_222[i] = S_127_SCEL_222[i] + S_127_SCEL_222[i-1]

for i in xrange(1,len(S_127_SCEL_222_INTER)):
    S_127_SCEL_222_INTER[i] = S_127_SCEL_222_INTER[i] + S_127_SCEL_222_INTER[i-1]


#Now, the normalization of the entropy to the number of K points = 64

S_127_SCEL_222 = S_127_SCEL_222 / 16.0
S_127_SCEL_222_INTER = S_127_SCEL_222_INTER / 128.0

output_array_127_SCEL_222 = np.vstack((nu_127_SCEL_222, S_127_SCEL_222)).T
output_array_127_SCEL_222_INTER = np.vstack((nu_127_SCEL_222_INTER, S_127_SCEL_222_INTER)).T

np.savetxt('/home/david/Trabajo/structures/PDOS_at_different_volumes/Calcite_I/pob_TZVP/PBE-D3/Only_Veq/3_times_landau_supercell/nu_S.dat', output_array_127_SCEL_222, header="nu\tS", fmt="%0.12g")
np.savetxt('/home/david/Trabajo/structures/PDOS_at_different_volumes/Calcite_I/pob_TZVP/PBE-D3/Only_Veq/3_times_landau_supercell/INTERPHESS_222/nu_S.dat', output_array_127_SCEL_222_INTER, header="nu\tS", fmt="%0.12g")

folder_222_All_freq = '/home/david/Trabajo/structures/PDOS_at_different_volumes/Calcite_I/pob_TZVP/PBE-D3/Only_Veq/3_times_landau_supercell'
folder_222_All_freq_INTER = '/home/david/Trabajo/structures/PDOS_at_different_volumes/Calcite_I/pob_TZVP/PBE-D3/Only_Veq/3_times_landau_supercell/INTERPHESS_222'

fig = plt.figure()
for indx_vols, indx_files in zip(range(1, len(vols)+1), range(len(vols))):
        print 'indx_vols init = ' , indx_vols
        print 'indx_files init = ' , indx_files
        ax = fig.add_subplot(len(vols), 1, indx_vols)
        ax.text(.85,.65, 'Calcite I,\n2 times Landau supercell\nPBE-D3, pob-TVPZ\nV$_{eq}$ = ' + vols[indx_vols-1] + ' $\AA^{3}$', fontsize=7, #alpha=1.0,
            bbox=dict(alpha=1.0, facecolor='white', pad=3.0), #, facecolor='red', edgecolor='red'),
            horizontalalignment='center',
            transform=ax.transAxes)
        x_222, y_222 = np.loadtxt(os.path.join(folder_222, files_222_0_3[indx_files]), skiprows = 4).T
        x_222_INTER, y_222_INTER = np.loadtxt(os.path.join(folder_222, folder_222_INTERPHESS_222, files_222_0_3_INTER[indx_files]), skiprows = 4).T


        nu_SCEL_222, S_SCEL_222 = np.loadtxt(os.path.join(folder_222_All_freq, 'nu_S.dat'), skiprows = 1).T
        nu_SCEL_222_INTER, S_SCEL_222_INTER = np.loadtxt(os.path.join(folder_222_All_freq_INTER, 'nu_S.dat'), skiprows = 1).T

        lns1 = ax.plot(x_222, y_222, color='lightskyblue', label='No interpolation') 
        lns1INTER = ax.plot(x_222_INTER, y_222_INTER, color='dodgerblue', label='INTERPHESS 2 2 2') 

        ax2 = ax.twinx()
        lns2 = ax2.plot(nu_SCEL_222, S_SCEL_222, 'lightskyblue', linestyle='--') #, label='No interpolation')
        lns2INTER = ax2.plot(nu_SCEL_222_INTER, S_SCEL_222_INTER, 'dodgerblue', linestyle='--') #, label='INTERPHESS 2 2 2')

        ax.grid()


#       lns = lns1+lns2+lns1INTER+lns2INTER
        lns = lns1+lns1INTER
        labs = [l.get_label() for l in lns]
        ax.legend(lns, labs, loc='right', fontsize=9)

        ax.set_ylabel('PDOS (states/cm$^{-1}$/cell)')
        ax2.set_ylabel('Entropy (mHartree/(cell$\cdot$K))\nat T = 298.15 K')
        ax.tick_params(labelbottom='on') 


        if indx_vols == len(vols) and indx_files == len(vols)-1:
           print 'enter 2nd if loop'
           print 'indx_files if indx_vols == len(vols)+1 and indx_files == len(vols): = ' , indx_files
           ax.set_xlabel('Energy (cm$^{-1}$)')

        else:
            ax.tick_params(labelbottom='on')


#lns = lns1+lns2+lns3
#labs = [l.get_label() for l in lns]
#ax.legend(lns, labs, loc=0)


fig.savefig("PDOS_0_and_8_plus_cum_entropy_T_298K.pdf", bbox_inches='tight', pad_inches=0.3)#, tight_layout() )#, bbox_inches=bbox)


plt.show()
sys.exit






ooooooooooooooooooo

vols_0 = vols_all[0]
vols_8 = vols_all[8]
vols = [vols_0, vols_8]

print 'vols = ', vols
files_222_0_8 = [files_222[0], files_222[8]]
print 'files_222_0_8 = ', files_222_0_8
files_444_0_8 = [files_444[0], files_444[8]]
print 'files_444_0_8 = ', files_444_0_8

fig = plt.figure()
print 'len(vols) = ', len(vols)
for indx_vols, indx_files in zip(range(1, len(vols)+1), range(len(vols))):
        print 'indx_vols init = ' , indx_vols
        print 'indx_files init = ' , indx_files
        ax = fig.add_subplot(len(vols), 1, indx_vols)
#       ax.set_title('V = ' + vols[indx_vols-1], fontsize=10) # This title stays in the center, and outside the graph
        ax.text(.85,.80,'V = ' + vols[indx_vols-1] + ' $\AA^{3}$', fontsize=10, #alpha=1.0,
            bbox=dict(alpha=1.0, facecolor='white', pad=1.0), #, facecolor='red', edgecolor='red'),
#           bbox=dict(facecolor='none', edgecolor='blue', pad=10.0, alpha=1.0),
#           alpha=1.0,
            horizontalalignment='center',
            transform=ax.transAxes)
        x_222, y_222 = np.loadtxt(os.path.join(folder_222, files_222_0_8[indx_files]), skiprows = 4).T
        x_444, y_444 = np.loadtxt(os.path.join(folder_444, files_444_0_8[indx_files]), skiprows = 4).T

        nu_SCEL_444, S_SCEL_444 = np.loadtxt(os.path.join(folder_444_All_freq, vols[indx_files], 'nu_S.dat'), skiprows = 1).T
        nu_SCEL_222, S_SCEL_222 = np.loadtxt(os.path.join(folder_222_All_freq, vols[indx_files], 'nu_S.dat'), skiprows = 1).T

        ax.plot(x_222, y_222, color='m', label='PDOS SCELPHONO 2x2x2') # black lines OK
        ax.plot(x_444, y_444, color='k', label='PDOS SCELPHONO 4x4x4') # magenta lines OK

        ax2 = ax.twinx()
        ax2.plot(nu_SCEL_222, S_SCEL_222, 'cornflowerblue', linestyle='--', label='Entropy SCELPHONO 2x2x2')
        ax2.plot(nu_SCEL_444, S_SCEL_444, 'grey', linestyle='--', label='Entropy SCELPHONO 4x4x4')

        ax.grid()
#       ax.set_ylim(0.01, 1.01) 
#       ylabels = [0.0, 0.5, 1.0] 
#       ax.set_yticklabels(ylabels)#, #,rotation=0),
#                 verticalalignment='baseline')#,
#                 horizontalalignment='left')

        if indx_vols == 1 and indx_files == 0:
           print 'enter 1st if loop'
           print 'indx_vols == 1 and indx_files == 0: ', indx_files 
           ax.legend(bbox_to_anchor=(0,1.2,1,0.6), loc="lower left",
                mode="expand",
                borderaxespad=0, 
                ncol=4)
           ax2.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand",
                borderaxespad=0, 
                ncol=4)

           ax.set_ylabel('PDOS (states/cm$^{-1}$/cell)')
           ax2.set_ylabel('Entropy (mHartree/(cell$\cdot$K))\nat T = 298.15 K')
           ax.tick_params(labelbottom='on') 


        if indx_vols == len(vols) and indx_files == len(vols)-1:
           print 'enter 2nd if loop'
           print 'indx_files if indx_vols == len(vols)+1 and indx_files == len(vols): = ' , indx_files
           ax.set_xlabel('Energy (cm$^{-1}$)')

        else:
            ax.tick_params(labelbottom='on')

fig.savefig("PDOS_0_and_8_plus_cum_entropy_T_298K.pdf", bbox_inches='tight', pad_inches=0.3)#, tight_layout() )#, bbox_inches=bbox)


################### CONSTANTS   ###############
global KB, h, c, T

# KB = boltmann cte, KB = 1.38064852(79)x10-23 J/K
KB = 1.38064852E-23

# h = plank constant, h = 6.626070040(81)x10-34 J s
h = 6.626070040E-34

# T = temperature, T = 2000.0 K
T = 2000.0

# c = speed of light, c = 2.99792458E8 m/s
c = 2.99792458E+8

def S(nu):
 S = -KB * np.log(1 - np.exp(-h * nu * 1E+2 * c  / (KB*T)))     + ( (h/T) * ( nu * 1E+2 * c * ( (np.exp(h *  nu  * 1E+2 * c  / (KB*T)) - 1)**(-1) )    )  )

 # Conversion: S above this line is in mHartree/K:
 return S  *((1/4.3597482)*1E+18 * 1E+3)

S_127_SCEL_444 = S(nu_127_SCEL_444)
S_116_SCEL_444 = S(nu_116_SCEL_444)

S_127_SCEL_222 = S(nu_127_SCEL_222)
S_116_SCEL_222 = S(nu_116_SCEL_222)

for i in xrange(1,len(S_127_SCEL_444)):
    S_127_SCEL_444[i] = S_127_SCEL_444[i] + S_127_SCEL_444[i-1]

for i in xrange(1,len(S_116_SCEL_444)):
    S_116_SCEL_444[i] = S_116_SCEL_444[i] + S_116_SCEL_444[i-1]

for i in xrange(1,len(S_127_SCEL_222)):
    S_127_SCEL_222[i] = S_127_SCEL_222[i] + S_127_SCEL_222[i-1]

for i in xrange(1,len(S_116_SCEL_222)):
    S_116_SCEL_222[i] = S_116_SCEL_222[i] + S_116_SCEL_222[i-1]

#Now, the normalization of the entropy to the number of K points = 64
S_127_SCEL_444 = S_127_SCEL_444 / 64.0
S_116_SCEL_444 = S_116_SCEL_444 / 64.0

S_127_SCEL_222 = S_127_SCEL_222 / 8.0
S_116_SCEL_222 = S_116_SCEL_222 / 8.0

output_array_127_SCEL_444 = np.vstack((nu_127_SCEL_444, S_127_SCEL_444)).T
output_array_116_SCEL_444 = np.vstack((nu_116_SCEL_444, S_116_SCEL_444)).T

output_array_127_SCEL_222 = np.vstack((nu_127_SCEL_222, S_127_SCEL_222)).T
output_array_116_SCEL_222 = np.vstack((nu_116_SCEL_222, S_116_SCEL_222)).T

np.savetxt('/home/david/Trabajo/structures/PDOS_at_different_volumes/Calcite_I/CI_SCELPHONO_444_Volumes_outputs/127.054446/nu_S_at_T_2000K.dat', output_array_127_SCEL_444, header="nu\tS at T = 2000.0K", fmt="%0.12g")
np.savetxt('/home/david/Trabajo/structures/PDOS_at_different_volumes/Calcite_I/CI_SCELPHONO_444_Volumes_outputs/116.573346/nu_S_at_T_2000K.dat', output_array_116_SCEL_444, header="nu\tS at T = 2000.0K", fmt="%0.12g")

np.savetxt('/home/david/Trabajo/structures/PDOS_at_different_volumes/Calcite_I/CI_SCELPHONO_222_Volumes_outputs/127.054446/nu_S_at_T_2000K.dat', output_array_127_SCEL_222, header="nu\tS at T = 2000.0K", fmt="%0.12g")
np.savetxt('/home/david/Trabajo/structures/PDOS_at_different_volumes/Calcite_I/CI_SCELPHONO_222_Volumes_outputs/116.573346/nu_S_at_T_2000K.dat', output_array_116_SCEL_222, header="nu\tS at T = 2000.0K", fmt="%0.12g")

folder_444_All_freq = '/home/david/Trabajo/structures/PDOS_at_different_volumes/Calcite_I/CI_SCELPHONO_444_Volumes_outputs'
folder_222_All_freq = '/home/david/Trabajo/structures/PDOS_at_different_volumes/Calcite_I/CI_SCELPHONO_222_Volumes_outputs'

vols_0 = vols_all[0]
vols_8 = vols_all[8]
vols = [vols_0, vols_8]

print 'vols = ', vols
files_222_0_8 = [files_222[0], files_222[8]]
print 'files_222_0_8 = ', files_222_0_8
files_444_0_8 = [files_444[0], files_444[8]]
print 'files_444_0_8 = ', files_444_0_8

fig = plt.figure()
print 'len(vols) = ', len(vols)
for indx_vols, indx_files in zip(range(1, len(vols)+1), range(len(vols))):
        print 'indx_vols init = ' , indx_vols
        print 'indx_files init = ' , indx_files
        ax = fig.add_subplot(len(vols), 1, indx_vols)
#       ax.set_title('V = ' + vols[indx_vols-1], fontsize=10) # This title stays in the center, and outside the graph
        ax.text(.85,.80,'V = ' + vols[indx_vols-1] + ' $\AA^{3}$', fontsize=10, #alpha=1.0,
            bbox=dict(alpha=1.0, facecolor='white', pad=1.0), #, facecolor='red', edgecolor='red'),
#           bbox=dict(facecolor='none', edgecolor='blue', pad=10.0, alpha=1.0),
#           alpha=1.0,
            horizontalalignment='center',
            transform=ax.transAxes)
        x_222, y_222 = np.loadtxt(os.path.join(folder_222, files_222_0_8[indx_files]), skiprows = 4).T
        x_444, y_444 = np.loadtxt(os.path.join(folder_444, files_444_0_8[indx_files]), skiprows = 4).T

        nu_SCEL_444, S_SCEL_444 = np.loadtxt(os.path.join(folder_444_All_freq, vols[indx_files], 'nu_S_at_T_2000K.dat'), skiprows = 1).T
        nu_SCEL_222, S_SCEL_222 = np.loadtxt(os.path.join(folder_222_All_freq, vols[indx_files], 'nu_S_at_T_2000K.dat'), skiprows = 1).T

        ax.plot(x_222, y_222, color='m', label='PDOS SCELPHONO 2x2x2') # black lines OK
        ax.plot(x_444, y_444, color='k', label='PDOS SCELPHONO 4x4x4') # magenta lines OK

        ax2 = ax.twinx()
        ax2.plot(nu_SCEL_222, S_SCEL_222, 'pink', linestyle='--', label='Entropy SCELPHONO 2x2x2')
        ax2.plot(nu_SCEL_444, S_SCEL_444, 'grey', linestyle='--', label='Entropy SCELPHONO 4x4x4')

        ax.grid()
#       ax.set_ylim(0.01, 1.01) 
#       ylabels = [0.0, 0.5, 1.0] 
#       ax.set_yticklabels(ylabels)#, #,rotation=0),
#                 verticalalignment='baseline')#,
#                 horizontalalignment='left')

        if indx_vols == 1 and indx_files == 0:
           print 'enter 1st if loop'
           print 'indx_vols == 1 and indx_files == 0: ', indx_files 
           ax.legend(bbox_to_anchor=(0,1.2,1,0.6), loc="lower left",
                mode="expand",
                borderaxespad=0, 
                ncol=4)
           ax2.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand",
                borderaxespad=0, 
                ncol=4)

           ax.set_ylabel('PDOS (states/cm$^{-1}$/cell)')
           ax2.set_ylabel('Entropy (mHartree/(cell$\cdot$K))\nat T = 2000.0 K')
           ax.tick_params(labelbottom='on') 


        if indx_vols == len(vols) and indx_files == len(vols)-1:
           print 'enter 2nd if loop'
           print 'indx_files if indx_vols == len(vols)+1 and indx_files == len(vols): = ' , indx_files
           ax.set_xlabel('Energy (cm$^{-1}$)')

        else:
            ax.tick_params(labelbottom='on')

fig.savefig("PDOS_0_and_8_plus_cum_entropy_at_T_2000K.pdf", bbox_inches='tight', pad_inches=0.3)#, tight_layout() )#, bbox_inches=bbox)

plt.show()
