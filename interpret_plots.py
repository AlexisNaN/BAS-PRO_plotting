import numpy as np
import interpret_tools
import interpret_cdf
from spacepy import pycdf
from datetime import datetime
import sys
from scipy import optimize as opti
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import matplotlib.dates as mdates
import matplotlib.patheffects as path_effects

# font sizes and styles:--------------------------------+
titleSize = 16
yaxisSize = 16
xaxisSize = 16
ytickSize = 14
xtickSize = 14
legSize = 14
labSize = 15

# set the default font family like so:
matplotlib.rc('font', family='Arial')

# we can also define a new font
import matplotlib.font_manager as font_manager

cmusfont = {'fontname': 'Arial'}
cmusfontFM = font_manager.FontProperties(family='Arial',
                                         style='normal',
                                         size=11)
years = mdates.YearLocator()  # every year
months = mdates.MonthLocator()  # every month
bimonthly = mdates.MonthLocator(interval=2)
years_fmt = mdates.DateFormatter('%Y/%m')
# months_fmt = mdates.DateFormatter('%m')
# ----------------------------------------------------END



def plot_padist_panels(fname_cdf, energies, Lshells, fname_plot, nt = 2):
    plotshape = np.zeros((len(energies),len(Lshells)))

    plotenergies = np.array([energies]).T
    plotenergies = np.repeat(plotenergies, len(Lshells), axis = 1)
    plotLshells = np.array([Lshells])
    plotLshells = np.repeat(plotLshells, len(energies), axis = 0)

    cdf = pycdf.CDF(fname_cdf)

    dynamic = cdf.attrs[interpret_cdf.lab_dynamic][0]
    if dynamic > 0: #boolean is not preserved by CDF format (?)
        dynamic = True
    else:
        dynamic = False


    #which times to plot (for a dynamic run):
    ax_t = cdf[interpret_cdf.lab_axt]
    if dynamic:
        t_plot = np.linspace(ax_t[-1], ax_t[0], nt)
    else:
        t_plot = [ax_t[0]]
        nt = 1

    n = np.shape(plotshape)[0]
    m = np.shape(plotshape)[1]
    print("", "{}r x {}c plot".format(n, m))

    #fig, ax_array_all = plt.subplots(n, m, sharex=True)#, gridspec_kw={'width_ratios': [0.6]*m})
    fig, ax_array_all = plt.subplots(n, m+(m-1), sharex=True, sharey=False,
                                 gridspec_kw={'width_ratios': [1]+[0.075, 1]*(m-1)})
    cmap = matplotlib.cm.get_cmap('viridis')
    usegrid = True
    addlabel = True

    #arrange the axes:
    ax_array = []
    if n == 1:
        ax_array_all = [ax_array_all]
    for row in ax_array_all:
        ax_row = []
        if m == 1:
            ax_row.append(row)
        else:
            for count, col in enumerate(row):
                if (count % 2 == 0):
                    ax_row.append(col)
                else:
                    col.axis('off')
        ax_array.append(ax_row)


    for rowidx, ax_row in enumerate(ax_array):
        print("","starting row #{}".format(rowidx+1))
        for colidx, ax_col in enumerate(ax_row):
            print("","","starting col #{}".format(colidx+1))
            en = plotenergies[rowidx][colidx]
            L_ = plotLshells[rowidx][colidx]

            ax_target = ax_col

            #interpolate to the required time
            for time_plot in t_plot:

                #get PAD:
                plot_alpha, plot_f = interpret_tools.getpad_(cdf, L_, en, dynamic, time_plot)
                #plot_alpha, plot_f = getpad_(filesets_eachK, L_plot, en, dynamic, time_plot)
                plot_j = interpret_tools.f2j(en, plot_f)

                colour = interpret_tools.get_colour_from_time(time_plot, ax_t, cmap)
                # if not dynamic: #don't want scatter marks in the dynamic plots
                #     ax_target.scatter(plot_alpha, plot_j, color=colour,# linewidth=1.4,
                #                label=str(round(en,1))+" MeV", alpha=0.75,marker='x',s=15)
                ax_target.plot(plot_alpha, plot_j, color=colour, linewidth=1.4, label="{:.1f} MeV".format(en), alpha=1)


                miny = min(plot_j[plot_j>0])

            x1 = ax_target.get_xlim()[0]
            x2 = ax_target.get_xlim()[1]
            y1 = ax_target.get_ylim()[0]
            y2 = ax_target.get_ylim()[1]

            ax_target.set_ylim([0, y2])

            #x axis tidy up:
            ax_target.set_xlim([interpret_tools.get_lc(Lshells[-1]), 90])

            start, end = ax_target.get_xlim()
            da = 20
            xtickarray = np.arange(end, start - da, -da)
            ax_target.xaxis.set_ticks(xtickarray)

            ax_target.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter(useMathText=True, useOffset=False))
            ax_target.yaxis.get_offset_text()
            exponent = ax_target.yaxis.get_offset_text()
            exponent.set_fontsize(14)
            ax_target.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

            usegrid = True
            if usegrid:
                ax_target.grid(linestyle='--', color='grey')


            # axis ticks:
            ax_target.yaxis.set_ticks_position('default')
            ax_target.yaxis.get_offset_text().set_fontsize(xtickSize-2)
            #ax_target.xaxis.set_tick_params(direction='out', length=5, labelsize=xtickSize-2)  # , which='bottom'
            ax_target.yaxis.set_tick_params(direction='out', length=2, labelsize=ytickSize-2, right = False)  # , which='bottom'
            ax_target.tick_params(labelsize=xtickSize-2)
            ax_target.set_ylabel('')


            if addlabel:
                ax_row[colidx].text(0.05,0.85,"L = {:.2f}, {:.2f}MeV".format(L_, en), transform=ax_row[colidx].transAxes) 

            #ax_row[colidx].yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter(useMathText=True, useOffset=False))
        #ax_row[0].set_ylabel('f [km$^{-6}$s$^{3}$]', fontdict={'fontsize': yaxisSize})


    # Additional plot information:
    #aeq on x axes:
    for colidx, ax_col in enumerate(ax_array[-1]):
        ax_col.set_xlabel('$\\alpha_{eq} [^{\\circ}]$', fontdict={'fontsize': xaxisSize})


    #E labels:
    for rowidx, ax_row in enumerate(ax_array):
        en = plotenergies[rowidx][0]
        ax_row[0].text(-0.45, 0.5, "{} MeV".format(en) + '\n'+'j [cm$^{-2}$s$^{-1}$str$^{-1}$MeV$^{-1}$]', transform=ax_row[0].transAxes,
                va='center',ha='center',fontdict={'fontsize': yaxisSize}, rotation=90)
    #L labels:
    for colidx, ax_col in enumerate(ax_array[0]):
        L_ = plotLshells[0][colidx]
        ax_col.text(0.5, 1.22, "L = {}".format(L_), transform=ax_col.transAxes,
                va='center',ha='center',fontdict={'fontsize': yaxisSize}, rotation=0)


    if dynamic:
        #label each dynamic simulation:
        shadow = False
        textax = ax_array[-1][-1]
        for idx_t, time_plot in enumerate(t_plot):
            colour = interpret_tools.get_colour_from_time(time_plot, ax_t, cmap)
            time_dt = datetime.fromtimestamp(time_plot)
            if len(t_plot) > n*12:
                labelheight = n*(idx_t * 1.0/(len(t_plot)-1))
            else:
                labelheight = n*(idx_t * 1.0/(n*12))

            text = textax.text(1.1, labelheight, time_dt.strftime("%j/%Y"), transform=textax.transAxes, color=colour,
                va='center',fontdict={'fontsize': yaxisSize-4})

            if shadow:
                text.set_path_effects([path_effects.Stroke(linewidth=0.8, foreground='black'),
                                       path_effects.Normal()])


    fig.set_size_inches(m*2.1, n*3.2)
    fig.set_size_inches(m*4.1, n*3.2)
    #plt.savefig(fname_plot, dpi = 400)
    plt.savefig(fname_plot + ".pdf", format= 'pdf')
    print("","saved figure to", fname_plot)

    return 1


def plot_f_vs_L_panels(fname_cdf, z_, iKs, fname_plot, fixedenergy = False, nt = 2):
    z_ = np.array(z_)
    plotshape = np.zeros(np.shape(z_))

    cdf = pycdf.CDF(fname_cdf)

    dynamic = cdf.attrs[interpret_cdf.lab_dynamic][0]
    if dynamic > 0: #boolean is not preserved by CDF format (?)
        dynamic = True
    else:
        dynamic = False

    #which times to plot (for a dynamic run):
    ax_t = cdf[interpret_cdf.lab_axt]
    if dynamic:
        t_plot = np.linspace(ax_t[-1], ax_t[0], nt)
    else:
        t_plot = [ax_t[0]]
        nt = 1

    n = np.shape(plotshape)[0]
    m = np.shape(plotshape)[1]
    print("", "{}r x {}c plot".format(n, m))

    #set up plot:
    fig, ax_array_all = plt.subplots(n, m+(m-1), sharex=True, sharey=False,
                                 gridspec_kw={'width_ratios': [1]+[0.0, 1]*(m-1)})
    scale_axes_factor = 3.35
    scale_axes_factor = 5.3
    width = 0
    for axwidth in [1]+[0.0, 1]*(m-1): width += axwidth
    width = width * scale_axes_factor #scale factor
    fig.set_size_inches(width, scale_axes_factor*n/1.6)
    usegrid = True
    cmap = matplotlib.cm.get_cmap('viridis')
    plot_energy_bars = True
    #legendon = False


    #arrange the axes:
    ax_array = []
    if n == 1:
        ax_array_all = [ax_array_all]

    for row in ax_array_all:
        ax_row = []
        if m == 1:
            ax_row.append(row)
        else:
            for count, col in enumerate(row):
                if (count % 2 == 0):
                    ax_row.append(col)
                else:
                    col.axis('off')
        ax_array.append(ax_row)


    ax_t = cdf[interpret_cdf.lab_axt]
    ax_mu = cdf[interpret_cdf.lab_axmu]
    ax_K = cdf[interpret_cdf.lab_axK]
    ax_L = cdf[interpret_cdf.lab_axL]
    map_alpha = cdf[interpret_cdf.lab_map]


    ylim_max = 0
    for rowidx, ax_row in enumerate(ax_array):
        print("","starting row #{}".format(rowidx+1))
        for colidx, ax_col in enumerate(ax_row):
            print("","","starting col #{}".format(colidx+1))
            en = z_[rowidx][colidx]
            mu = en #could be mu or energy supplied
            ax_target = ax_col


            #------------------------------------------------------------------------------+
            # make plot                                                                    |
            #------------------------------------------------------------------------------+
            minys = []
            maxys = []
            minys_zoom = []
            maxys_zoom = []

            ax_target.set_yscale('log')

            # if zoombox:
            #     #zoom box:
            #     ax_target_zoom = zoomed_inset_axes(ax_target, 2, loc=4, # zoom = 6
            #         bbox_to_anchor = [1, 0.075], bbox_transform =ax_target.transAxes)
            #     #mark_inset(ax_target, ax_target_zoom, loc1=2, loc2=4, fc="none", ec="red", alpha=0.5, lw=1, zorder = 1)
            #     mark_inset(ax_target, ax_target_zoom, loc1=2, loc2=4, fc="none", ec="black", alpha=0.5, lw=1, ls="solid")


            #     x1_zoom = 1.19
            #     x2_zoom = 1.35

            for idx_K in iKs:
                if idx_K > len(ax_K) - 1:
                    continue
                K_now = ax_K[idx_K]
                sol_en = cdf[interpret_cdf.lab_en][0, :, idx_K, :]
                sol_f = cdf[interpret_cdf.lab_f][:, :, idx_K, :]

                idxLvf = 0
                while min(sol_en[:, idxLvf]) < 0:
                    idxLvf+=1
                    if idxLvf==len(ax_L):
                        break
                if idxLvf>=len(ax_L)-2:
                    print("","warning: cannot plot at K={:.2f}".format(K_now))
                    break

                if (fixedenergy):
                    for idxL in range(idxLvf, len(ax_L)):
                        if (en > sol_en[-1, idxL] or en < sol_en[0, idxL]):
                            print("Error: energy {}MeV is out of bounds".format(en))
                            sys.exit(1)


                for time_plot in t_plot: #do it backward so we can stick our K label to t0 solution plotted
                    colour = interpret_tools.get_colour_from_time(time_plot, ax_t, cmap)

                    if dynamic: #interpolate to the current t we need:
                        #get idx_t_0 and idx_t_1 surrounding the time we want to plot:
                        idx_t_0 = -1
                        idx_t_1 = idx_t_0
                        if time_plot < ax_t[0] or time_plot > ax_t[-1]:
                            print("Error: time_plot is out of range on K idx", idx_K)
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

                        sol_f1d_t_0 = []
                        sol_f1d_t_1 = []
                        if (not fixedenergy):
                            sol_en1d = [] #used for displaying energy bars when inspectmu == True
                            #get f at every L at the energy under investigation:
                            # t0
                            for idxL in range(idxLvf, len(ax_L)):
                                #sol_f1d_t_0.append(float(np.interp(np.log10(mu),ax_mu[:], sol_f[idx_t_0*len(ax_mu):(1+idx_t_0)*len(ax_mu), idxL])))
                                sol_f1d_t_0.append(np.interp(np.log10(mu), ax_mu[:], sol_f[idx_t_0, :, idxL]))
                                sol_en1d.append(np.interp(np.log10(mu), ax_mu[:], sol_en[:, idxL]))
                            

                            # t1
                            for idxL in range(idxLvf, len(ax_L)):
                                #sol_f1d_t_1.append(float(np.interp(np.log10(mu),ax_mu[:], sol_f[idx_t_1*len(ax_mu):(1+idx_t_1)*len(ax_mu), idxL])))  
                                #sol_en1d.append(float(np.interp(np.log10(mu), ax_mu[:], sol_en[:, idxL])))    
                                sol_f1d_t_1.append(np.interp(np.log10(mu), ax_mu[:], sol_f[idx_t_1, :, idxL]))

                        else:
                            #get f at every L at the energy under investigation:
                            # t0
                            for idxL in range(idxLvf, len(ax_L)):
                                #sol_f1d_t_0.append(float(np.interp(en, sol_en[:, idxL], sol_f[idx_t_0*len(ax_mu):(1+idx_t_0)*len(ax_mu), idxL])))
                                sol_f1d_t_0.append(np.interp(en, sol_en[:, idxL], sol_f[idx_t_0, :, idxL]))

                            # t1
                            for idxL in range(idxLvf, len(ax_L)):
                                #sol_f1d_t_1.append(float(np.interp(en, sol_en[:, idxL], sol_f[idx_t_1*len(ax_mu):(1+idx_t_1)*len(ax_mu), idxL])))
                                sol_f1d_t_1.append(np.interp(en, sol_en[:, idxL], sol_f[idx_t_1, :, idxL]))


                        sol_f1d = []
                        ax_t_surround = [ax_t[idx_t_0], ax_t[idx_t_1]]
                        for idxL in range(len(sol_f1d_t_0)): #interpolate f at each L from surrounding times
                            f_surround = [sol_f1d_t_0[idxL], sol_f1d_t_1[idxL]]
                            sol_f1d.append(np.interp(time_plot, ax_t_surround, f_surround))

                        
                    else:
                        #ask user whether to plot f vs. mu or energy:
                        sol_f1d = []
                        if (not fixedenergy):
                            sol_en1d = [] #used for displaying energy bars when inspectmu == True
                            for idxL in range(idxLvf,len(ax_L)):
                                #sol_f1d.append(float(np.interp(np.log10(mu),ax_mu[:],sol_f[:,idxL])))
                                sol_f1d.append(np.interp(np.log10(mu), ax_mu[:], sol_f[0, :, idxL]))
                                sol_en1d.append(np.interp(np.log10(mu), ax_mu[:], sol_en[:, idxL]))
                        else:
                            for idxL in range(idxLvf,len(ax_L)):
                                #sol_f1d.append(float(np.interp(en,sol_en[:,idxL],sol_f[:,idxL])))
                                sol_f1d.append(np.interp(en,sol_en[:,idxL],sol_f[0, :, idxL]))
                        

                    sol_f1d = np.array(sol_f1d)
                    sol_j1d = interpret_tools.f2j(en, sol_f1d)

                    #if not (labelarr):
                    if (not fixedenergy):
                        label = "{:.0f}MeV/G".format(mu)
                        sol_plot = sol_f1d
                    else:
                        label = "{:.0f}MeV".format(en)
                        sol_plot = sol_j1d


                    if not len(sol_plot[sol_plot>0]): continue

                    zorder = len(ax_K) - idx_K
                    ax_L_plot = ax_L[idxLvf:]
                    

                    ax_target.plot(ax_L_plot, sol_plot, color=colour, linewidth=0.8, label=label, alpha=1, zorder = zorder)
                    #save the y axis limits:
                    minys.append(np.min(sol_plot[sol_plot>0]))
                    maxys.append(np.max(sol_plot))

                    # if zoombox:
                    #     ax_target_zoom.plot(ax_L_plot,sol_plot, color=colour, linewidth=1, label=label, alpha=1,zorder = zorder)
                    #     #save the y axis limits:
                    #     minys_zoom.append(np.min(sol_plot[ax_L_plot>=x1_zoom]))
                    #     maxys_zoom.append(np.max(sol_plot[ax_L_plot<=x2_zoom]))


                #label K
                ax_target.text(ax_L[-1]+0.02,sol_plot[-1],"{:.3f}".format(K_now),rotation=0,
                    color='black', size=9,ha="left",va="center")#,transform=ax_col.transAxes)


            #------------------------------------------------------------------------------+
            # adjust axes, etc.                                                            |
            #------------------------------------------------------------------------------+


            #automatically detect limits:
            xmin_select = ax_target.get_xlim()[0]
            xmax_select = ax_target.get_xlim()[1]
            ymin_select = max(minys)
            ymax_select = max(maxys)
            

            major_ticks = np.power(10.0,np.linspace(-10,10,21))

            ax_target.set_yticks(major_ticks)
            ax_target.minorticks_off()

            ax_target.set_xlim([ax_L[0],ax_L[-1]])
            ax_target.set_ylim([ymin_select,ymax_select])

            # if zoombox:
            #     #adjust zoomed axes:
            #     ax_target_zoom.set_yscale('log')
            #     ax_target_zoom.set_xlim(x1_zoom, x2_zoom)
            #     ax_target_zoom.set_ylim(min(minys_zoom), max(maxys_zoom))

            #write K= 
            ax_target.text(1.02,1.035,"$K=$",rotation=0,
                transform=ax_target.transAxes,
                color='black', size=9,ha="left",va="center")#,transform=ax_col.transAxes)



            if usegrid:
                ax_target.grid(which='minor', alpha=0.2)
                ax_target.grid(which='major', alpha=0.5)
                # if zoombox:
                #     ax_target_zoom.grid(which='minor', alpha=0.2)
                #     ax_target_zoom.grid(which='major', alpha=0.5)

            #------------------------------------------------------------------------------+
            # print a set of lines indicating energy for given mu                          |
            #------------------------------------------------------------------------------+

            #enplot = np.array([0.5, 0.75, 1, 1.5, 2, 5, 7.5, 10, 15, 20, 30, 45, 65, 80, 100])
            enplot_L = np.array([1.2,1.4,1.6,1.8,2.0])

            if plot_energy_bars and not fixedenergy:
                sol_en = cdf[interpret_cdf.lab_en][0, :, 0, :]
                sol_f = cdf[interpret_cdf.lab_f][0, :, 0, :] #t0, K = 0
                idxLvf = 0
                while min(sol_en[:, idxLvf]) < 0:
                    idxLvf+=1
                    if idxLvf==len(ax_L):
                        break
                if idxLvf>=len(ax_L)-2:
                    print("","warning: cannot plot energy bars at K=0")
                else:
                    sol_f1d = []
                    sol_en1d = [] #used for displaying energy bars when inspectmu == True
                    for idxL in range(idxLvf,len(ax_L)):
                        sol_f1d.append(np.interp(np.log10(mu),ax_mu[:],sol_f[:, idxL]))
                        sol_en1d.append(np.interp(np.log10(mu), ax_mu[:], sol_en[:, idxL]))
                    sol_f1d = np.array(sol_f1d)
                    ax_en = np.array(sol_en1d)

                    ax_L_nonfill = ax_L[idxLvf:]#[ax_en>0]

                    enplot_E = np.interp(enplot_L,ax_L_nonfill,ax_en)
                    for idx in range(len(enplot_L)):
                        ysol = np.interp(enplot_L[idx], ax_L_nonfill, sol_f1d)
                        ysol_norm = (np.log10(ysol)-np.log10(ymin_select))/(np.log10(ymax_select)-np.log10(ymin_select))
                        ax_target.axvline(enplot_L[idx],ymin = ysol_norm,
                            ymax=1, linestyle=":",color='black',zorder=0,linewidth=1.1,alpha=0.85)

                        labelx = (enplot_L[idx]-xmin_select)/(xmax_select - xmin_select)
                        if idx == len(enplot_L)-1:
                            ha = "right"
                        else:
                            ha = "center"
                        ax_target.text(labelx,1.035,'{:.1f}'.format(enplot_E[idx]),rotation=0,
                            color='dimgrey', transform=ax_target.transAxes,size=9,ha=ha,va="center")




            # text:
            #ax_col.set_title("$\\mu$ = " + str(mu) + "MeV/G", fontdict={'fontsize': yaxisSize - 1}, rotation=0)
            if not fixedenergy:
                ax_target.text(0.035, 0.975,"$\\mu=$" + "{:.1f}MeV/G".format(mu), fontdict={'fontsize': labSize}, rotation=0,
                            horizontalalignment='left', verticalalignment='top',
                            transform=ax_col.transAxes)
            else:
                ax_target.text(0.035, 0.975,"$E=$" + "{:.1f}MeV".format(en), fontdict={'fontsize': labSize}, rotation=0,
                            horizontalalignment='left', verticalalignment='top',
                            transform=ax_col.transAxes)
            # axis ticks:
            ax_target.set_ylabel('')
            ax_target.set_xlabel('')
            ax_target.yaxis.set_ticks_position('both')
            ax_target.xaxis.set_tick_params(direction='out', length=5, labelsize=xtickSize)#, which='bottom'
            ax_target.yaxis.set_tick_params(direction='out', length=2, labelsize=ytickSize)#, which='bottom'
            ax_target.tick_params(labelsize=xtickSize)

            ylims = ax_col.get_ylim()

        if not fixedenergy:
            ax_row[0].set_ylabel('$f$ [km$^{-6}$s$^{3}$]', fontdict={'fontsize': yaxisSize-2}) #m_{0}^{3}
        else:
            ax_row[0].set_ylabel('$j$ [cm$^{-2}$s$^{-1}$str$^{-1}$MeV$^{-1}$]', fontdict={'fontsize': yaxisSize-2})


    for ax in ax_array[n-1]:
        ax.set_xlabel('$L$', fontdict={'fontsize': yaxisSize})

    if not fixedenergy:
        ax_array[0][0].text(-0.2,1.035,'$E$ (MeV)\nat $K=0$',rotation=0,
            color='dimgrey', transform=ax_array[0][0].transAxes,size=9,ha="center",va="center")
        ax.arrow(-0.122,1.035, 0.09, 0, transform=ax_array[0][0].transAxes, head_width=0.022, head_length=0.01,
            fc='dimgrey', ec='dimgrey', clip_on=False)


    #write K units:
    xlim=ax_array[-1][-1].get_xlim()
    ax_array[-1][-1].text(1.02,0.,"K units\nG$^{0.5}R_{E}$",rotation=0,
        transform=ax_array[-1][-1].transAxes,
        color='black', size=9,ha="left",va="bottom")#,transform=ax_col.transAxes)


    if dynamic:
        #label each dynamic simulation:
        shadow = False
        textax = ax_array[-1][-1]
        for idx_t, time_plot in enumerate(t_plot):
            colour = interpret_tools.get_colour_from_time(time_plot, ax_t, cmap)
            time_dt = datetime.fromtimestamp(time_plot)
            if len(t_plot) > n*12:
                labelheight = n*(idx_t * 1.0/(len(t_plot)-1))
            else:
                labelheight = n*(idx_t * 1.0/(n*12))
            text = textax.text(1.03, labelheight, time_dt.strftime("%j/%Y"), transform=textax.transAxes, color=colour,
                va='center',fontdict={'fontsize': yaxisSize-4})
            if shadow:
                text.set_path_effects([path_effects.Stroke(linewidth=0.8, foreground='black'),
                                       path_effects.Normal()])


    #plt.savefig(fname_plot, dpi = 400)
    plt.savefig(fname_plot + ".pdf", format= 'pdf')
    print("","saved figure to", fname_plot)

    return 1

def plot_f_vs_L_mu_panels(cdf, mus, iKs, fname_plot, nt): plot_f_vs_L_panels(cdf, mus, iKs, fname_plot, fixedenergy = False, nt = nt)
def plot_f_vs_L_E_panels(cdf, energies, iKs, fname_plot, nt): plot_f_vs_L_panels(cdf, energies, iKs, fname_plot, fixedenergy = True, nt = nt)


def plot_n_vs_L_panels(fname_cdf, energies, L_axis_plot, fname_plot, nt = 2):
    #reshape the energies array:
    plotenergies = np.array(energies)
    plotshape = np.zeros(np.shape(plotenergies))
    cdf = pycdf.CDF(fname_cdf)

    dynamic = cdf.attrs[interpret_cdf.lab_dynamic][0]
    if dynamic > 0: #boolean is not preserved by CDF format (?)
        dynamic = True
    else:
        dynamic = False


    #which times to plot (for a dynamic run):
    ax_t = cdf[interpret_cdf.lab_axt]
    if dynamic:
        t_plot = np.linspace(ax_t[-1], ax_t[0], nt)
    else:
        t_plot = [ax_t[0]]
        nt = 1

    n = np.shape(plotshape)[0]
    m = np.shape(plotshape)[1]
    print("", "{}r x {}c plot".format(n, m))
    sharey = False
    if n == 1: sharey = True
    fig, ax_array_all = plt.subplots(n, m, sharex=True, sharey=sharey)#, gridspec_kw={'width_ratios': [0.6]*m})
    cmap = matplotlib.cm.get_cmap('viridis')
    usegrid = True

    #arrange the axes:
    ax_array = []
    if n == 1:
        ax_array_all = [ax_array_all]
    for row in ax_array_all:
        ax_row = []
        if m == 1:
            ax_row.append(row)
        else:
            for count, col in enumerate(row):
                ax_row.append(col)
        ax_array.append(ax_row)

    ylim_max = 0
    #print(t_plot)
    for rowidx, ax_row in enumerate(ax_array):
        print("","starting row #{}".format(rowidx+1))
        for colidx, ax_col in enumerate(ax_row):
            en = plotenergies[rowidx][colidx]
            ax_target = ax_col
            print("","","starting col #{}".format(colidx+1))

            #interpolate to the required time
            for time_plot in t_plot:
                #print(time_plot)
                n_ = []
                ratio = []
                for L_ in L_axis_plot:
                    alpha_lc = interpret_tools.get_lc(L_)

                    alpha_eqpad_deg, f_ = interpret_tools.getpad_(cdf, L_, en, dynamic, time_plot)

                    j_ = interpret_tools.f2j(en,f_)

                    # ----- CORRECTION FOR MISSING LOSS CONE DATA POINT AT INTERPOLATED L ----- #
                    correct_missing_lc_datapoint = True
                    if correct_missing_lc_datapoint:
                        j_ = j_[1:]
                        alpha_eqpad_deg = alpha_eqpad_deg[1:]

                        j_ = np.array([0, 0] + list(j_))
                        alpha_eqpad_deg = np.array([0, alpha_lc] + list(alpha_eqpad_deg))
                    # ------------------------------------------------------------------------- #

                    j_tofit_sinnc = [alpha_eqpad_deg[2:], j_[2:]]

                    #get the x axis in radians:
                    alpha_eqpad_rad = np.radians(alpha_eqpad_deg)


                    p0 = [j_[-1], 4] #A, n
                    popt, pcov = opti.curve_fit(interpret_tools.f_sinn_simple, j_tofit_sinnc[0], j_tofit_sinnc[1], p0= p0,
                        bounds=((0, 0), (2*j_[-1], 1000)), maxfev=2000)


                    perr = np.sqrt(np.diag(pcov))
                    n_.append(popt[-1])
                    #a = f_sinn_simple(np.array([70.0,90.0]),*popt)
                    #ratio.append(a[1]/a[0])

                colour = interpret_tools.get_colour_from_time(time_plot, ax_t, cmap)
                ax_target.plot(L_axis_plot, n_, color=colour, linewidth=1.4,alpha=1)
                #ax_target.scatter(L_axis_plot, n_, color=colour,marker='o', alpha=1)
                #ax_target.plot(L_axis_plot, ratio, color=colour, linewidth=1.4,alpha=1,linestyle=":")
                #ax_target.scatter(L_axis_plot, ratio, color=colour,marker='o', alpha=1)

            ax_target.set_xlim([L_axis_plot[0],L_axis_plot[-1]])
            if usegrid: ax_target.grid(linestyle='--', color='grey')



            
            #format axis ticks:
            #ax_row[colidx].set_yscale('log')
            ax_target.yaxis.set_ticks_position('both')
            ax_target.xaxis.set_tick_params(direction='out', length=5, labelsize=xtickSize)  # , which='bottom'
            ax_target.yaxis.set_tick_params(direction='out', length=2, labelsize=ytickSize)  # , which='bottom'
            ax_target.tick_params(labelsize=xtickSize)
            ax_target.set_ylabel('')

            #keep a record of the max value:
            ylim = ax_target.get_ylim()
            if (ylim[1] > ylim_max):
                ylim_max = ylim[1]


            #if addlabel:
            #    ax_row[colidx].text(0.1,0.1,"L = {}, E = {}".format(Lplot, en), transform=ax_row[colidx].transAxes) #DELETE THIS 

            #ax_row[colidx].yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter(useMathText=True, useOffset=False))

        #ax_row[0].set_ylabel('f [km$^{-6}$s$^{3}$]', fontdict={'fontsize': yaxisSize})


    #add additional plot information:
    #aeq on x axes:
    for colidx, ax_col in enumerate(ax_array[-1]):
        ax_col.set_xlabel('$L$', fontdict={'fontsize': xaxisSize})

    plt.tight_layout(pad=5, w_pad=0, h_pad=0)
    

    for rowidx, ax_row in enumerate(ax_array):
        for colidx, ax_col in enumerate(ax_row):
            en = plotenergies[rowidx][colidx]

            ax_col.text(0.99, 0.89, "{} MeV".format(en), transform=ax_col.transAxes,
                    va='center',ha='right',fontdict={'fontsize': yaxisSize})
            ax_col.yaxis.set_major_locator(MultipleLocator(20))
            ax_col.yaxis.set_minor_locator(MultipleLocator(10))
            if ylim_max <20:
                ax_col.yaxis.set_major_locator(MultipleLocator(10))
            ax_col.grid(linestyle='--', color='grey',which='both',alpha=0.6,linewidth=1.3)
            ax_col.set_ylim([0,ylim_max])

            #set_size(ax_col, 2*(figsx/m),2*(figsy/n)*(ylim_max/60))


    if dynamic:
        #label each dynamic simulation:
        shadow = False
        textax = ax_array[-1][-1]
        for idx_t, time_plot in enumerate(t_plot):
            colour = interpret_tools.get_colour_from_time(time_plot, ax_t, cmap)
            time_dt = datetime.fromtimestamp(time_plot)
            if len(t_plot) > n*12:
                labelheight = n*(idx_t * 1.0/(len(t_plot)-1))
            else:
                labelheight = n*(idx_t * 1.0/(n*12))
            text = textax.text(1.005, labelheight, time_dt.strftime("%j/%Y"), transform=textax.transAxes, color=colour,
                va='center',fontdict={'fontsize': yaxisSize-4})
            if shadow:
                text.set_path_effects([path_effects.Stroke(linewidth=0.8, foreground='black'),
                                       path_effects.Normal()])
    

    
    figsx = m*4
    figsy = n*2
    fig.set_size_inches(figsx, figsy)
    fig.tight_layout() 
    
    
    #plt.savefig(fname_plot, dpi = 400)
    plt.savefig(fname_plot + ".pdf", format= 'pdf')
    print("","saved figure to", fname_plot)

    return 1