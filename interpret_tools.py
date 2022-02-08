import interpret_cdf
import colorsys
import sys
import numpy as np
from math import sqrt, asin, pi


mass0_proton = 1.6726219e-27
MeV2J = 1.60218e-13
c_ = 299792458

def userselectkey(filedict, allowmulti=False):
    #user selects keys from a dictionary of items supplied as argument
    #input is sanitised and returned as a list of keys
    
    for key in filedict.keys():
        print(key, '...', filedict[key])

    filedict_selectkeys = [] #save the selection as keys
    
    #ask the user to selet results files to plot:
    decided = False
    more = False
    while not decided:
        choice = input("> ")
        #sanitise the input and re-ask if necessary:
        try:
            if choice[-1] == ",":
                more = True * allowmulti
                choice = choice[:-1]
            else:
                more = False
            choice = int(choice)
            if choice in filedict.keys():
                if choice in filedict_selectkeys:
                    print ("  already selected")
                else:
                    filedict_selectkeys.append(choice)
                if not more:
                    decided = True
            else:
                print("  out of range")
        except:
            print("  invalid")
            pass
            
    return filedict_selectkeys

def userinputfloat(allowblank=False, allowmulti=False):
    floatselectkeys = [] #save the selection as keys
    
    decided = False
    more = False
    while not decided:
        choice = input("> ")
                
        if choice == "":
            if allowblank:
                decided = True
                floatselectkeys.append(False)
            else:
                print("  a value must be specified")
        else:
            #sanitise the input and re-ask if necessary:
            try:
                if choice[-1] == ",":
                    more = True * allowmulti
                    choice = choice[:-1]
                else:
                    more = False

                choice = float(choice)
                
                if choice in floatselectkeys:
                    print ("  already selected")
                else:
                    floatselectkeys.append(choice)
                    
                if not more:
                    decided = True
            except:
                print("  invalid")
                pass
    return floatselectkeys


def get_N_HexCol(N):
    HSV_tuples = [(x*1.0/N, 0.5, 0.5) for x in range(N)]
    RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
    return RGB_tuples

def getlastline(fname):
    with open(fname, 'r') as f:
        lines = f.read().splitlines()
    last_line = lines[-1]
    return last_line.strip("\n")

# def bytes2float(bytestring):
#     return float(((bytestring).decode("utf-8")).strip(' ').strip('\n'))
# def bytes2str(bytestring):
#     return str(((bytestring).decode("utf-8")).strip(' ').strip('\n'))

# def eqL2B(Larray):
#     B0 = 3.12e-5
#     T2G = 1.0e4
#     return T2G*B0*(np.power(np.reciprocal(Larray),3))

def getcol(n):
    plotnshuffle = 0 #shift the linestyle/colour scheme/label of the set of lines
    n = plotnshuffle + n
    #ncols = [x for x in get_N_HexCol(len(filedict_plotkey))]
    ncols = ['blue',"#3cb44b",'deepskyblue']
    return ncols[n%len(ncols)]

def get_axis_file(fname):
    with open(fname) as fi:
        ax_ = fi.readlines()
        for idxL in range(0, len(ax_)):
            ax_[idxL] = float(ax_[idxL].strip('/n'))
    ax_ = np.array(ax_)
    return ax_

def get_sol_file(fname):
    lines = []
    with open(fname) as fo:
        lines = fo.readlines()

    sol_f = []
    for idx in range(0, len(lines)):
        sol_f.append([float(x) for x in lines[idx].strip('\n').split()])

    sol_f = np.array(sol_f)
    return sol_f

def get_gamma(ke_MeV):
    #calculate gamma of a proton given its kinetic energy in MeV

    ke_J = ke_MeV * MeV2J

    erest_J = mass0_proton*c_*c_

    gamma = 1 + ke_J/(erest_J)
    return gamma

def get_pr(ke_MeV):
    # !calculate relativistic momentum in SI units (kg m s-1) of a proton given its kinetic energy in MeV


    erest_J = mass0_proton*c_*c_
    
    gamma = get_gamma(ke_MeV)

    mr = mass0_proton * gamma
    etot_J = mr*c_*c_
    
    p_ = np.sqrt((etot_J**2) - (erest_J**2))/c_ #relativistic momentum kg m/s
    return p_

def f2j(energy, psd):
    # !
    # ! input units of energy: MeV
    # ! input units of phase space density: m-6 s3 kg-3
    # !
    # !

    #change units to m-6 s3: 
    temp = psd / 1e18
    # units: m-6 s3
    
    
    temp = temp / (mass0_proton**3.)
    # units are m-6 s3 kg-3

    #get momentum squared:
    p2 = (get_pr(energy) **2.)
    # units: kg2 m2 s-2, or: ** kg J **

    flux = temp * p2
    # units: m-2 s-1 str-1 J-1

    flux = flux / 1e4
    # units: cm-2 s-1 str-1 J-1

    flux = flux * MeV2J
    # units: cm-2 s-1 str-1 MeV-1
    return flux

def j2f(energy, flux):
    # !
    # ! input units of energy: MeV
    # ! input units of flux: cm-2 s-1 str-1 MeV-1
    # !
    # !


    #get momentum squared:
    p2 = (get_pr(energy) **2.)
    # units: kg2 m2 s-2, or: ** kg J **


    flux = flux / MeV2J
    # units: cm-2 s-1 str-1 J-1
    
    flux = flux * 1e4
    # units: m-2 s-1 str-1 J-1
    

    temp = flux / p2
    
    temp = temp * (mass0_proton**3.)

    psd = temp * 1e18
    # units: m-6 s3 kg-3
    
    return psd


def get_colour_from_time(time, ax_t, cmap):
    if (len(ax_t) > 1):
        frac = (time-ax_t[0])/(ax_t[-1]-ax_t[0])
    else:
        frac=0
    return cmap(frac)

def f_sinn(x, A, b, c, n):
    # python does not like a negative number to a decimal power
    # p0 should be something like [10, 0, 25, 4] in practise
    d2r = np.pi / 180.
    sinn = np.abs(np.sin((x+b)*d2r))
    return A * np.power(sinn,n) + c

def f_sinn_simple(x, A, n):
    # python does not like a negative number to a decimal power
    # p0 should be something like [10, 0, 25, 4] in practise
    d2r = np.pi / 180.
    sinn = np.abs(np.sin((x)*d2r))
    return A * np.power(sinn,n)

def get_lc(Lb): #centred dipole loss cone approximation for 2015
    RE = 6.3712e6
    atm_height_std = 100000
    B0 = 2.986731323946967e-05

    ra = (RE + atm_height_std)/RE #~Earth's surface + atm_height_dipolelc m
    
    if ra >= Lb:
        return 90
    else:
        Ba = (B0/(ra**3)) * (4 - 3*ra/Lb)**(0.5)
        dipole_lc = asin(sqrt((B0 / Lb**3)/Ba)) * 180 / pi
        return dipole_lc

def getpad_(cdf, L_extract, en, dynamic, time_plot):
    #
    #
    #   WARNING: when interpolating between L output by the model, we can't determine where the loss cone is!
    #    so the PAD will be returned with a point at (0,0) but then shoot up to the first point outside the lc
    #
    #
    #returns pitch angle distribution alpha and corresponding f, straight from the solution grid

    ax_t = cdf[interpret_cdf.lab_axt]
    ax_mu = cdf[interpret_cdf.lab_axmu]
    ax_K = cdf[interpret_cdf.lab_axK]
    ax_L = cdf[interpret_cdf.lab_axL]
    map_alpha = cdf[interpret_cdf.lab_map]

    sol_f1d_allK = []
    sol_alpha_allK = []

    for idx_K in range(len(ax_K)):
        K_now = ax_K[idx_K]
        sol_en = cdf[interpret_cdf.lab_en][0, :, idx_K, :]
        sol_f = cdf[interpret_cdf.lab_f][:, :, idx_K, :]
        #find the minimum L that is outside the loss cone at the current K:
        idxL_outsidelc = np.argwhere(map_alpha[:,idx_K]>0)[0][0]

        #check if L is out of range for interpolation at this K:
        if L_extract < np.min(ax_L[idxL_outsidelc:]) or L_extract > np.max(ax_L[idxL_outsidelc:]):
            continue

        #get alpha from the map file, but ignore fill values:
        sol_alpha = np.interp(L_extract, ax_L[idxL_outsidelc:], map_alpha[:,idx_K][idxL_outsidelc:])

        #check the energy array too: we may have defined pitch angle at every K for nL (model artifact)

        if np.sum(sol_en[:,idxL_outsidelc:len(ax_L)]<0.) > 0:
        #if there are any elements below 0 in sol_en:
            continue

        #show a warning if energy is out of range for interpolation:
        for idxL in range(idxL_outsidelc, len(ax_L)):
            if (en > sol_en[-1, idxL] or en < sol_en[0, idxL]):
                print("Error: energy {:.2f}MeV is out of bounds at alpha={:.2f} (iK={})".format(en,sol_alpha,idx_K+1))
                sys.exit(1)



        if dynamic: #interpolate to the current t we need:
            #get idx_t_0 and idx_t_1 surrounding the time we want to plot:
            idx_t_0 = -1
            idx_t_1 = idx_t_0
            if time_plot < ax_t[0] or time_plot > ax_t[-1]:
                print("","Error: time_plot is out of range on K idx", idx_K)
                sys.exit(1)
            for idx_t, time_sol in enumerate(ax_t):
                if time_plot >= time_sol:
                    idx_t_0 = idx_t
                    if time_plot == time_sol:
                        idx_t_1 = idx_t_0 #we 'interpolate' across 1 grid
                    else:
                        idx_t_1 = idx_t + 1
                else:
                    break

            #idx_extent_f_required = [idx_t_0*len(ax_mu), (1+idx_t_1)*len(ax_mu)] #from the first row of idx_t_0 to last of idx_t_1
            #idx_ofst = idx_extent_f_required[0]

            #sol_f = unpack_fileset_fonly_part(fileset, idx_extent_f_required)
            

            #get f at every L at the energy under investigation:
            # t0
            sol_f1d_t_0 = []
            for idxL in range(idxL_outsidelc, len(ax_L)):
                #sol_f1d_t_0.append(float(np.interp(en, sol_en[:, idxL], sol_f[idx_t_0*len(ax_mu)-idx_ofst:(1+idx_t_0)*len(ax_mu)-idx_ofst, idxL])))
                sol_f1d_t_0.append(np.interp(en, sol_en[:, idxL], sol_f[idx_t_0, :, idxL]))

            # get f at the L under investigation:
            sol_f1d_t_0 = np.interp(L_extract, ax_L[idxL_outsidelc:], sol_f1d_t_0)

            if not (ax_t[idx_t_0] == ax_t[idx_t_1] ): #skip interpolating from the second surrounding time
                # t1
                sol_f1d_t_1 = []
                for idxL in range(idxL_outsidelc, len(ax_L)):
                    #sol_f1d_t_1.append(float(np.interp(en, sol_en[:, idxL], sol_f[idx_t_1*len(ax_mu)-idx_ofst:(1+idx_t_1)*len(ax_mu)-idx_ofst, idxL])))
                    sol_f1d_t_1.append(np.interp(en, sol_en[:, idxL], sol_f[idx_t_1, :, idxL]))
                # get f at the L under investigation:
                sol_f1d_t_1 = np.interp(L_extract, ax_L[idxL_outsidelc:], sol_f1d_t_1)

                # interpolate to t:
                ax_t_surround = [ax_t[idx_t_0], ax_t[idx_t_1]]
                f_surround = [sol_f1d_t_0, sol_f1d_t_1]
                sol_f1d = np.interp(time_plot, ax_t_surround, f_surround)
            else:
                sol_f1d = sol_f1d_t_0
        else:
            idx_t = 0
            #get f at every L at the energy under investigation:
            sol_f1d = []
            for idxL in range(idxL_outsidelc, len(ax_L)):
                #sol_f1d.append(float(np.interp(en, sol_en[:, idxL], sol_f[idx_t*len(ax_mu):(1+idx_t)*len(ax_mu), idxL])))
                sol_f1d.append(np.interp(en, sol_en[:, idxL], sol_f[idx_t, :, idxL]))
                
            #get f at the L under investigation:
            sol_f1d = np.interp(L_extract, ax_L[idxL_outsidelc:], sol_f1d)

       
        #print (sol_alpha, sol_f1d_t_0, sol_f1d_t_1)
        sol_f1d_allK.append(sol_f1d)
        sol_alpha_allK.append(sol_alpha)


    #add zero to the beginning of the array and reverse it to have K=0 last
    sol_f = np.array([0.]+sol_f1d_allK[::-1])
    sol_alpha = np.array([0.]+sol_alpha_allK[::-1])


    return sol_alpha, sol_f
