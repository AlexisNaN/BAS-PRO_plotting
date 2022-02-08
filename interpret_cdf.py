from spacepy import pycdf
import datetime
import numpy as np
import interpret_tools
import os
from os import listdir
from os.path import isfile, join, exists
import sys
from datetime import datetime
from math import pi

#global variables specifying CDF data variables:
lab_name = 'Name'
lab_dynamic = 'Dynamic'
lab_auth = 'Author' 
lab_date = 'CreateDate'
lab_axmu = 'axis_mu'
lab_axK = 'axis_K'
lab_axL = 'axis_L'
lab_axt = 'axis_t'
lab_axtd = 'axis_t_date'
lab_map = 'map_KL-aeq'
lab_f = 'f'
lab_j = 'j'
lab_axphi = 'axis_phi'
lab_en = 'energy'


def get_phi_2015(L_):
	"""Get the average dipole field strength around Earth's equator and dipole moment. Use like so: B0,m = get_B0_m(2000.0)"""
	#from pyIGRF.loadCoeffs import get_coeffs
	#g, h = get_coeffs(2015)
	g10 = -29441.46 #g[1][0]
	g11 = -1501.77 #g[1][1]
	h11 = 4795.99 #h[1][1]

	B0_2 = g10**2 + g11**2 + h11**2
	B0_ = B0_2**0.5
	B0_ = B0_*1e-9 # = 2.986731323946967e-05 T

	RE = 6.3712e6 #m
	#mu0 = 1.25663706e-6 #

	phi = 2 * pi * B0_ * (RE ** 2) / L_ #T m2
	return phi


def convert_to_cdf(wd_current, auto_choose = [-1, -1], overwrite = False):

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
	if auto_choose[0] >= 0:
		idx_sol = auto_choose[0]
	else:
		print("Please select a solution to convert (0 most recent):")
		idx_sol = interpret_tools.userselectkey(dir_solution_dict,allowmulti=False)[0]

	#move the current working directory to that solution:
	solname = dir_solution_dict[idx_sol]
	prefix = solname[:15]


	#get the user to select a solution to plot:
	if auto_choose[1] >= 0:
		dyn_yesno = auto_choose[1]
	else:
		print("Plot dynamic solution?")
		dir_yesno = {0:'no', 1:'yes'}
		dyn_yesno = interpret_tools.userselectkey(dir_yesno,allowmulti=False)[0]

	dyn = False
	ext_txt = ".txt"
	ext_cdf = ".cdf"
	if dyn_yesno > 0:
		dyn = True
		ext_txt = "_dyn.txt"
		ext_cdf = "_dyn.cdf"


	#check if a previous CDF file exists:
	fname_cdf = join(wd_current, solname, prefix + '_solution' + ext_cdf)
	file_exists = os.path.exists(fname_cdf)
	print("Converting solution code {} to CDF".format(prefix))
	if file_exists:
		if overwrite:
			print("","warning: CDF file already exists, overwriting:")
			print("", fname_cdf)
			os.remove(fname_cdf)
		else:
			print("","warning: CDF file already exists, not updating")
			return fname_cdf			
		#print("Error: CDF file already exists:", fname_cdf)
		#sys.exit(1)
	else:
		print("","file will be output to:")
		print("", fname_cdf)





	fname_axis_K = join(wd_current, solname, prefix + "_axis_K"+ext_txt)
	Kaxisdata = interpret_tools.get_axis_file(fname_axis_K)
	solution_files_eachK = {}
	for idx_K, Kvalue in enumerate(Kaxisdata):
	    fileset = {}
	    fileset['K'] = Kvalue
	    fileset['axis mu'] = join(wd_current, solname, prefix + "_axis_mu"+ext_txt)
	    fileset['axis K'] = join(wd_current, solname, prefix + "_axis_K" + ext_txt)
	    fileset['axis L'] = join(wd_current, solname, prefix + "_axis_L" + ext_txt)
	    fileset['axis t'] = join(wd_current, solname, prefix + "_axis_t"+ext_txt)
	    fileset['config'] = join(wd_current, solname, prefix + "_config"+ext_txt)
	    fileset['map'] = join(wd_current, solname, prefix + "_map_iK-aeq"+ext_txt)
	    kslice_prefix = "_iK-" + str(idx_K+1).zfill(4) + "_2D_"
	    fileset['f'] = join(wd_current, solname, prefix + kslice_prefix + "f"+ext_txt)
	    fileset['en'] = join(wd_current, solname, prefix + kslice_prefix + "en"+ext_txt)
	    fileset['data'] = join(wd_current, solname, prefix + kslice_prefix + "f_data"+ext_txt)
	    fileset['sdev'] = join(wd_current, solname, prefix + kslice_prefix + "f_sdev"+ext_txt)
	    fileset['prog'] = join(wd_current, solname, prefix + "_progress"+ext_txt)
	    solution_files_eachK[idx_K] = fileset


	#read in axis files:
	axis_mu = interpret_tools.get_axis_file(solution_files_eachK[0]['axis mu'])
	axis_K = interpret_tools.get_axis_file(solution_files_eachK[0]['axis K'])
	axis_L = interpret_tools.get_axis_file(solution_files_eachK[0]['axis L'])
	axis_t = interpret_tools.get_axis_file(solution_files_eachK[0]['axis t'])
	if not dyn: axis_t = np.array([axis_t[0]])
	map_aeq = interpret_tools.get_sol_file(solution_files_eachK[0]['map'])
	# convert timestamps to datetime objects (CDF_EPOCH):
	axis_t_date = []
	for ts in axis_t:
		axis_t_date.append(datetime.fromtimestamp(ts))
	axis_t_date = np.array(axis_t_date)
	axis_phi = np.array([get_phi_2015(L_) for L_ in axis_L])


	#open CDF file:
	cdf = pycdf.CDF(fname_cdf, '')
	cdf.attrs[lab_name] = solname
	cdf.attrs[lab_dynamic] = dyn
	cdf.attrs[lab_auth] = 'A. R. Lozinski [alezin33@bas.ac.uk]'
	cdf.attrs[lab_date] = datetime.now()
	# copy axes:
	cdf[lab_axmu] = axis_mu
	cdf[lab_axmu].attrs['units'] = 'log10(mu/ (1 MeV/G) )'
	cdf[lab_axK] = axis_K
	cdf[lab_axK].attrs['units'] = 'G^0.5 RE'
	cdf[lab_axL] = axis_L
	cdf[lab_axt] = axis_t
	cdf[lab_axt].attrs['units'] = 'seconds since Jan 01 1970 (UTC)'
	cdf[lab_axtd] = axis_t_date
	cdf[lab_axtd].attrs['units'] = 'datetime CDF object'
	cdf[lab_map] = map_aeq
	cdf[lab_map].attrs['units'] = 'degrees'
	cdf[lab_axphi] = axis_phi
	cdf[lab_axphi].attrs['units'] = 'T m2'


	#create data grids:
	data_f = np.zeros((np.size(axis_t), np.size(axis_mu), np.size(axis_K), np.size(axis_L)))
	energygrid = np.zeros((1, np.size(axis_mu), np.size(axis_K), np.size(axis_L)))
	print("","data grid dimensions (t, mu, K, L): ", np.shape(data_f))

	for idx_K, Kvalue in enumerate(Kaxisdata):
		soldata_fixedK = interpret_tools.get_sol_file(solution_files_eachK[idx_K]['f'])
		engrid_fixedK = interpret_tools.get_sol_file(solution_files_eachK[idx_K]['en'])
		energygrid[0, :, idx_K, :] = engrid_fixedK[:, :]
		for idx_t in range(len(axis_t)):
			data_f[idx_t, :, idx_K, :] = soldata_fixedK[idx_t * np.size(axis_mu) : (idx_t + 1) * (np.size(axis_mu))]

	cdf[lab_f] = data_f
	cdf[lab_f].attrs['desc'] = 'distribution function: relativistic phase space density multiplied by proton rest mass cubed'
	cdf[lab_f].attrs['units'] = 'km-6 s3'

	cdf[lab_en] = energygrid
	cdf[lab_en].attrs['desc'] = 'energy at each data coordinate, assumed constant in time (ignoring secular variation)'
	cdf[lab_en].attrs['units'] = 'MeV'


	#translate pdf to flux too:
	data_j = np.zeros((np.size(axis_t), np.size(axis_mu), np.size(axis_K), np.size(axis_L)))
	for idx_t in range(len(axis_t)):
		# for idx_mu in range(len(axis_mu)):
		# 	for idx_K in range(len(axis_K)):
		# 		for idx_L in range(len(axis_L)):
		# 			energy = energygrid[0, idx_mu, idx_K, idx_L]
		# 			if energy <= 0:
		# 				continue
		# 			data_j[idx_t, idx_mu, idx_K, idx_L] = interpret_tools.f2j(energy, data_f[idx_t, idx_mu, idx_K, idx_L])
		data_j[idx_t, energygrid[0]>0] = interpret_tools.f2j(energygrid[0, energygrid[0]>0], data_f[idx_t, energygrid[0]>0])
	cdf[lab_j] = data_j
	cdf[lab_j].attrs['desc'] = 'unidirectional differential proton flux'
	cdf[lab_j].attrs['units'] = 'cm-2 s-1 str-1 MeV-1'

	cdf.close()
	print("","done")
	return fname_cdf
	#test points:
	#print(data_f[0, 350, 20, :][data_f[0, 350, 20, :]>0])







#wd_current = "/Users/goldfish/Dropbox/Work/Projects/BAS-PRO/output/"
#convert_to_cdf(wd_current, auto_choose=[0,0])

# interpret_cdf.lab_name
# interpret_cdf.lab_dynamic
# interpret_cdf.lab_auth
# interpret_cdf.lab_date
# interpret_cdf.lab_axmu
# interpret_cdf.lab_axK
# interpret_cdf.lab_axL
# interpret_cdf.lab_phi
# interpret_cdf.lab_axt
# interpret_cdf.lab_axtd
# interpret_cdf.lab_map
# interpret_cdf.lab_f
# interpret_cdf.lab_j
# interpret_cdf.lab_en




