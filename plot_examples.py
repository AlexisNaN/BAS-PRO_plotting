import numpy as np
import interpret_tools
import interpret_plots
import interpret_cdf
import os
import sys
from math import log10
from os import listdir
from os.path import isfile, join
from spacepy import pycdf


wd_current = "./"

plot_profile = True
plot_pad = True
plot_anis = True


#scan for solution directories:
dir_solution_dict = {}
dir_counter = 0
dir_content = os.listdir(wd_current)
dir_content_ordered = sorted(dir_content,reverse=True)
for x in dir_content_ordered:
    if os.path.isdir(join(wd_current, x)):
        if 'solution' in x:
            dir_solution_dict[dir_counter] = x
            dir_counter = dir_counter + 1


#get the user to select a solution to plot:
print("Please select a solution to plot (0 most recent):")
idx_sol = interpret_tools.userselectkey(dir_solution_dict,allowmulti=False)[0]
prefix = dir_solution_dict[idx_sol][:15]

#get the user to select whether the dynamic or static CDF should be plotted:
print("Plot dynamic solution? (takes longer)")
dir_yesno = {0:'no, static', 1:'yes, dynamic'}
dyn_yesno = interpret_tools.userselectkey(dir_yesno,allowmulti=False)[0]
if dyn_yesno > 0:
    dyn = True
else:
    dyn = False
print()


#create a CDF of the solution, if it doesn't aleady exist:
fname_cdf = interpret_cdf.convert_to_cdf(wd_current, auto_choose=[idx_sol, dyn_yesno], overwrite = False)
print()
#

#print a summary of the CDF:
print("{} overview:".format(fname_cdf))
cdf = pycdf.CDF(fname_cdf)
keys = list(cdf.keys())
for key in keys:
    attrs = list(cdf[key].attrs)
    print("", key)
    print("", "", "shape:", np.shape(cdf[key]))
    for attr in attrs:
        print("", "", "{}:".format(attr), cdf[key].attrs[attr])
    print()
print()



print("Plotting solution code", prefix)
# create a new subdir for the plots:
wd_local_plots = join(wd_current, "plots", prefix)
plotver = 0
while(os.path.isdir("{}_{}".format(wd_local_plots,plotver))):
    plotver+=1
wd_local_plots = "{}_{}".format(wd_local_plots,plotver)
print("","creating",wd_local_plots)
os.mkdir(wd_local_plots)
# wd_local_plots_pdfs = join(wd_local_plots, "pdfs")
# print("","creating",wd_local_plots_pdfs)
# os.mkdir(wd_local_plots_pdfs)
print()



######################################################################################################################
######################################################################################################################

#
#make plots of psd profiles at fixed mus:
#

if plot_profile:
    print("Plotting profiles...")

    #specify mu and energy in an array with the dimensions of the plot (in terms of nrows, ncols):
    nt = 25
    iKs = [0, 5]
    mus = [[20, 50], [100, 400]] #MeV/G, for plot_f_vs_L_mu_panels
    energies = [[5, 10], [15, 19]] #MeV, for plot_f_vs_L_E_panels
    
    fname_plot = join(wd_local_plots, prefix + '_f_vs_L_mu')#_nK-'+str(len(iKs)).zfill(4))
    interpret_plots.plot_f_vs_L_mu_panels(fname_cdf, mus, iKs, fname_plot, nt = nt)

    fname_plot = join(wd_local_plots, prefix + '_j_vs_L_E')#_nK-'+str(len(iKs)).zfill(4))
    interpret_plots.plot_f_vs_L_E_panels(fname_cdf, energies, iKs, fname_plot, nt = nt)

    print()

######################################################################################################################
######################################################################################################################

#
#make plots of pitch angle distributions at constant energies:
#

if plot_pad:
    print("Plotting PADs...")

    nt = 5
    energies = [2.5, 7.5, 12.5] #MeV
    Lshells = [1.2, 1.3, 1.5, 1.7, 2.0]

    figname = join(wd_local_plots, prefix + '_PAD')
    interpret_plots.plot_padist_panels(fname_cdf, energies, Lshells, figname, nt = nt)

    print()

######################################################################################################################
######################################################################################################################

#
#make plots of pitch angle n at constant energies:
#

if plot_anis:
    print("Plotting anistropy... (may take a long time)")
    
    #specify energies in an array with the dimensions of the plot (in terms of nrows, ncols):
    nt = 2
    energies = [[5, 10], [15, 19]]
    L_axis_plot = [1.2, 1.225, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2.0]
       
    fname_plot = join(wd_local_plots, prefix + "_anisotropy")
    interpret_plots.plot_n_vs_L_panels(fname_cdf, energies, L_axis_plot, fname_plot, nt = nt)

    print()
