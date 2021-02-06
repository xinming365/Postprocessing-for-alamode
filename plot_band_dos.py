#!/usr/bin/env python
#
# plotband.py
#
# Simple script to visualize phonon dispersion relations
#
# Copyright (c) 2014 Terumasa Tadano
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#

import numpy as np
import optparse
import os
import matplotlib as mpl
import mpl_toolkits
from matplotlib.gridspec import GridSpec
from plotdos import get_natoms_and_symbols, get_x_minmax, get_y_minmax, change_xscale,sum_atom_projected_dos
try:
    mpl.use("Qt5agg")
except:
    pass
import matplotlib.pyplot as plt



# parser options
usage = "usage: %prog [options] file1.bands file2.bands ... "
parser = optparse.OptionParser(usage=usage)

parser.add_option("--nokey", action="store_false", dest="print_key", default=True,
                  help="don't print the key in the figure")
parser.add_option("-u", "--unit", action="store", type="string", dest="unitname", default="kayser",
                  help="print the band dispersion in units of UNIT. Available options are kayser, meV, and THz", metavar="UNIT")
parser.add_option("--emin", action="store", type="float", dest="emin",
                  help="minimum value of the energy axis")
parser.add_option("--emax", action="store", type="float", dest="emax",
                  help="maximum value of the energy axis")
parser.add_option("--normalize", action="store_true", dest="normalize_xaxis", default=False,
                  help="normalize the x axis to unity.")


# font styles
mpl.rc('font', **{'family': 'Times New Roman','weight':'light', 'sans-serif': ['Helvetica']})
mpl.rc('xtick', labelsize=12)
mpl.rc('ytick', labelsize=16)
mpl.rc('axes', labelsize=16)
mpl.rc('lines', linewidth=1.5)
mpl.rc('legend', fontsize='small')
mpl.rcParams.update({"mathtext.fontset":'stix',"font.serif": ['SimSun'],})
# line colors and styles
# color = ['b', 'g', 'r', 'm', 'k', 'c', 'y', 'r']
color = np.array([(21,9,242), (64,26,217), (106,43,191),(149,60,166),(191,77,140),(234,94,115)])/255
# color = np.array([(21,9,242), (64,26,217), (106,43,191),(149,60,166),(191,77,140),(234,94,115),(0,0,0)])/255
lsty = ['-', '-', '-', '-', '-', '-','--', '--', '--']
# LightGrey: (211,211,211)
# Gray: (128,128,128)

def get_kpath_and_kval(file_in):

    ftmp = open(file_in, 'r')
    kpath = ftmp.readline().rstrip('\n').split()
    kval = ftmp.readline().rstrip('\n').split()
    ftmp.close()

    if kpath[0] == '#' and kval[0] == '#':
        kval_float = [float(val) for val in kval[1:]]
        kpath_list = []
        for i in range(len(kpath[1:])):
            if kpath[i + 1] == 'G':
                kpath_list.append('$\Gamma$')
            else:
                kpath_list.append("$\mathrm{%s}$" % kpath[i + 1])

        return kpath_list, kval_float
    else:
        return [], []


def change_scale(array, str_scale):

    str_tmp = str_scale.lower()

    if str_tmp == 'kayser':
        print("Band structure will be shown in units of cm^{-1}")
        return array

    elif str_tmp == 'mev':
        print("Band structure will be shown in units of meV")
        kayser_to_mev = 0.0299792458 * 1.0e+12 * \
            6.62606896e-34 / 1.602176565e-19 * 1000

        for i in range(len(array)):
            for j in range(len(array[i])):
                for k in range(1, len(array[i][j])):
                    array[i][j][k] *= kayser_to_mev

        return array

    elif str_tmp == 'thz':
        print("Band structure will be shown in units of THz")
        kayser_to_thz = 0.0299792458

        for i in range(len(array)):
            for j in range(len(array[i])):
                for k in range(1, len(array[i][j])):
                    array[i][j][k] *= kayser_to_thz

        return array

    else:
        print("Unrecognizable option for --unit %s" % str_scale)
        print("Band structure will be shown in units of cm^{-1}")
        return array


def normalize_to_unity(array, array_axis):

    for i in range(len(array)):
        max_val = array[i][-1][0]

        factor_normalize = 1.0 / max_val

        for j in range(len(array[i])):
            array[i][j][0] *= factor_normalize

    max_val = array_axis[-1]
    factor_normalize = 1.0 / max_val

    for i in range(len(array_axis)):
        array_axis[i] *= factor_normalize

    return array, array_axis


def get_xy_minmax(array):

    xmin, xmax, ymin, ymax = [0, 0, 0, 0]

    for i in range(len(array)):
        xtmp = array[i][-1][0]
        xmax = max(xmax, xtmp)

    for i in range(len(array)):
        for j in range(len(array[i])):
            for k in range(1, len(array[i][j])):
                ytmp = array[i][j][k]
                ymin = min(ymin, ytmp)
                ymax = max(ymax, ytmp)

    return xmin, xmax, ymin, ymax


def gridspec_setup(data_merged, xtickslabels, xticksvars):

    xmaxs = []
    xmins = []

    xticks_grids = []
    xticklabels_grids = []
    xticklabels_tmp = []
    xticks_tmp = []

    for i in range(len(xtickslabels)):

        if i == 0:
            xmins.append(xticksvars[0])
        else:
            if xticksvars[i] == xticksvars[i-1]:
                xmaxs.append(xticksvars[i - 1])
                xmins.append(xticksvars[i])

                xticks_grids.append(xticks_tmp)
                xticklabels_grids.append(xticklabels_tmp)
                xticklabels_tmp = []
                xticks_tmp = []

        xticklabels_tmp.append(xtickslabels[i])
        xticks_tmp.append(xticksvars[i])

    xticks_grids.append(xticks_tmp)
    xticklabels_grids.append(xticklabels_tmp)
    xmaxs.append(xticksvars[-1])

    naxes = len(xticks_grids)
    nfiles = len(data_merged)

    data_all_axes = []

    for i in range(naxes):
        data_ax = []

        xmin_ax = xmins[i]
        xmax_ax = xmaxs[i]

        for j in range(nfiles):

            kval = np.array(data_merged[j][0:, 0])
            ix_xmin_arr = np.where(kval <= xmin_ax)
            ix_xmax_arr = np.where(kval >= xmax_ax)

            if len(ix_xmin_arr[0]) > 0:
                ix_xmin = int(ix_xmin_arr[0][-1])
            else:
                ix_xmin = 0

            if len(ix_xmax_arr[0]) > 0:
                ix_xmax = int(ix_xmax_arr[0][0])
            else:
                ix_xmax = -2

            data_ax.append(data_merged[j][ix_xmin:(ix_xmax+1), :])

        data_all_axes.append(data_ax)

    return naxes, xticks_grids, xticklabels_grids, xmins, xmaxs, data_all_axes


def preprocess_data(files, unitname, normalize_xaxis):

    xtickslabels, xticksvars = get_kpath_and_kval(files[0])

    data_merged = []

    for file in files:
        data_tmp = np.loadtxt(file, dtype=float)
        data_merged.append(data_tmp)

    data_merged = change_scale(data_merged, unitname)

    if normalize_xaxis:
        data_merged, xticksvars = normalize_to_unity(data_merged, xticksvars)

    xmin, xmax, ymin, ymax = get_xy_minmax(data_merged)

    if options.emin is None and options.emax is None:
        factor = 1.05
        ymin *= factor
        ymax *= factor
    else:
        if options.emin is not None:
            ymin = options.emin
        if options.emax is not None:
            ymax = options.emax
        if ymin > ymax:
            print("Warning: emin > emax")

    naxes, xticks_grids, xticklabels_grids, xmins, xmaxs, data_merged_grids \
        = gridspec_setup(data_merged, xtickslabels, xticksvars)

    return naxes, xticks_grids, xticklabels_grids, \
        xmins, xmaxs, ymin, ymax, data_merged_grids


def run_plot(nax, xticks_ax, xticklabels_ax, xmin_ax, xmax_ax, ymin, ymax, data_merged_ax):

    fig = plt.figure()

    width_ratios = []
    used_colors = []
    for xmin, xmax in zip(xmin_ax, xmax_ax):
        width_ratios.append(xmax - xmin)

    gs = GridSpec(nrows=1, ncols=nax, width_ratios=width_ratios)
    gs.update(wspace=0.1)


    for iax in range(nax):
        ax = plt.subplot(gs[iax])

        for i in range(len(data_merged_ax[iax])):
            if len(data_merged_ax[iax][i]) > 0:
                ax.plot(data_merged_ax[iax][i][0:, 0], data_merged_ax[iax][i][0:, 1],
                        linestyle=lsty[i], color=color[i], label=files[i])
                used_colors.append(color[i])

                for j in range(2, len(data_merged_ax[iax][i][0][0:])):
                    ax.plot(data_merged_ax[iax][i][0:, 0], data_merged_ax[iax][i][0:, j],
                            linestyle=lsty[i], color=color[i])

        # # ax_0 = plt.subplot(gs[0,:])
        # # ax_0.set_axis_off()
        # cmp = mpl.colors.ListedColormap(used_colors)
        # norm = mpl.colors.BoundaryNorm(boundaries=bn, ncolors=cmp.N)
        # #fig.subplots_adjust(top=0.9)
        #
        # p0 = ax.get_position().get_points().flatten()
        # # l, b, w = [0.15, 0.92, 0.7]
        # # cax = fig.add_axes([0.15, 0.9, 0.2, 0.05])
        # # cax = fig.add_axes([p0[0], 1, p0[2]-p0[0], 0.05])
        # # cb_ax = mpl_toolkits.axes_grid1.inset_locator.inset_axes(ax, loc=3)
        # cb =fig.colorbar(mappable=mpl.cm.ScalarMappable(norm=norm, cmap=cmp), use_gridspec=False,
        #                  ax=ax, orientation='horizontal',shrink=1, fraction=0.1, pad=0.1)

        if iax == 0:
            if options.unitname.lower() == "mev":
                ax.set_ylabel("Frequency (meV)", labelpad=20)
            elif options.unitname.lower() == "thz":
                ax.set_ylabel("Frequency (THz)", labelpad=20)
            else:
                ax.set_ylabel("Frequency (cm${}^{-1}$)", labelpad=10)

        else:
            ax.set_yticklabels([])
            ax.set_yticks([])

        plt.axis([xmin_ax[iax], xmax_ax[iax], ymin, ymax])
        ax.set_xticks(xticks_ax[iax])
        ax.set_xticklabels(xticklabels_ax[iax])
        ax.xaxis.grid(True, linestyle='-')

        if options.print_key and iax == 0:
            ax.legend(loc='best', prop={'size': 10})

    plt.show()

def plot_dos(ax, dos_file, text_name,text_color_name, emax, emin):
    color = ['k', 'b', 'g', 'r', 'm', 'c', 'y', 'r',
             'darkred', 'darkblue', 'darkgreen', 'darkmagenta']
    lsty = ['-', '-', '-', '-', '--', '--', '--', '--', '-', '-', '-', '-']
    dos_lw = 1

    usage = "usage: %prog [options] file1.dos file2.dos ... "
    parser = optparse.OptionParser(usage=usage)

    parser.add_option("--pdos", action="store_true", dest="print_pdos", default=False,
                      help="print atom-projected phonon DOS")
    parser.add_option("--nokey", action="store_false", dest="print_key", default=True,
                      help="don't print the key in the figure")
    parser.add_option("-u", "--unit", action="store", type="string", dest="unitname", default="kayser",
                      help="print the band dispersion in units of UNIT. Available options are kayser, meV, and THz",
                      metavar="UNIT")
    parser.add_option("--emin", action="store", type="float", dest="emin",
                      help="minimum value of the energy axis")
    parser.add_option("--emax", action="store", type="float", dest="emax",
                      help="maximum value of the energy axis")

    # fakeArgs = ['-u', 'Thz', dos_file,"--pdos","--emax",str(emax),"--emin",'0']
    # fakeArgs = ['-u', 'Thz', dos_file,"--pdos","--emax",str(emax), "--emin", str(emin)]
    fakeArgs = ['-u', 'mev', dos_file, "--pdos", "--emax", str(emax), "--emin", str(emin)]
    options, args = parser.parse_args(fakeArgs)
    files = args[0:]
    nfiles = len(files)

    if nfiles == 0:
        print("Usage: plotdos.py [options] file1.dos file2.dos ...")
        print("For details of available options, please type\n$ python plotdos.py -h")
        exit(1)
    else:
        print("Number of files = %d" % nfiles)

    energy_axis = []
    dos_merged = []

    for file in files:
        data_tmp = np.loadtxt(file, dtype=float)
        energy_axis.append(data_tmp[:, 0])
        dos_merged.append(data_tmp[:, 1:])

    energy_axis = change_xscale(energy_axis, options.unitname)
    # 交换了x和y的表示！此处出现函数不一致。
    ymin, ymax= get_x_minmax(energy_axis)
    xmin, xmax = get_y_minmax(dos_merged)

    counter_line = 0

    for i in range(len(dos_merged)):
        counter_line = counter_line % 12

        # ax.plot(energy_axis[i][:], dos_merged[i][:, 0],
        #          linestyle=lsty[counter_line], color=color[counter_line], label="File" + str(i + 1) + ".Total")
        ax.plot(dos_merged[i][:, 0],energy_axis[i][:],
                 linestyle=lsty[counter_line], color=color[counter_line], label="Total", linewidth=dos_lw)

        counter_line += 1

        if options.print_pdos:
            symbols, natoms = get_natoms_and_symbols(files[i])

            if len(dos_merged[i][0, 1:]) != np.sum(natoms):
                print(
                    "Error: Projected DOS is not contained in the %d-th file" % (i + 1))
                exit(1)
            else:
                pdos = sum_atom_projected_dos(dos_merged[i][:, 1:], natoms)

                for j in range(len(pdos[0, :])):

                    # ax.plot(energy_axis[i][:], pdos[:, j], linestyle=lsty[counter_line],
                    #          color=color[counter_line], label="File" + str(i + 1) + "." + symbols[j])
                    ax.plot(pdos[:, j], energy_axis[i][:], linestyle=lsty[counter_line],
                             color=color[counter_line], label=symbols[j], linewidth=dos_lw)

                    counter_line += 1

    # if options.unitname.lower() == "mev":
    #     ax.set_xlabel("Frequency (meV)", fontsize=16)
    # elif options.unitname.lower() == "thz":
    #     ax.set_xlabel("Frequency (THz)", fontsize=16)
    # else:
    #     ax.set_xlabel("Frequency (cm${}^{-1}$)", fontsize=16)


    # ax.set_ylabel("Phonon DOS", fontsize=16, labelpad=20)
    ax.set_xlabel("DOS", fontsize=12)
    ax.text(0.26,0.5, text_name, transform=ax.transAxes, color=text_color_name)

    if options.emin == None and options.emax == None:
        factor = 1.00
        ymin *= factor
        ymax *= factor
    else:
        if options.emin != None:
            ymin = options.emin
        if options.emax != None:
            ymax = options.emax

        if ymin > ymax:
            print("Warning: emin > emax")

    # ymax *= 1.05
    # plt.axis([xmin, xmax, ymin, ymax])
    ax.set(xlim=(xmin,xmax), ylim=(ymin, ymax))

    ax.tick_params(axis='y', direction='in', labelsize=10,right=True)
    # ax.tick_params(axis='x', direction='in', labelsize=10, top=True)
    # ax.minorticks_on()
    # ax.tick_params(axis='both', direction='in', right=True, top=True)
    # ax.tick_params(axis='y',direction='in',right=True,which='minor')

    # ax.set_xticks(fontsize=16)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.set_xticks([])

    if options.print_key:
        ax.legend(loc='upper right', framealpha=0.6, handlelength=1.3, prop={'size': 8})




def plot_band(ax, nax, xticks_ax, xticklabels_ax, xmin_ax, xmax_ax, ymin, ymax, data_merged_ax):
    width_ratios = []
    used_colors = []
    for xmin, xmax in zip(xmin_ax, xmax_ax):
        width_ratios.append(xmax - xmin)

    gs = GridSpec(nrows=1, ncols=nax, width_ratios=width_ratios)
    gs.update(wspace=0.1)


    for iax in range(nax):

        for i in range(len(data_merged_ax[iax])):
            if len(data_merged_ax[iax][i]) > 0:
                ax.plot(data_merged_ax[iax][i][0:, 0], data_merged_ax[iax][i][0:, 1],
                        linestyle=lsty[i], color=color[i])
                used_colors.append(color[i])

                for j in range(2, len(data_merged_ax[iax][i][0][0:])):
                    ax.plot(data_merged_ax[iax][i][0:, 0], data_merged_ax[iax][i][0:, j],
                            linestyle=lsty[i], color=color[i])

        # # ax_0 = plt.subplot(gs[0,:])
        # # ax_0.set_axis_off()
        # cmp = mpl.colors.ListedColormap(used_colors)
        # norm = mpl.colors.BoundaryNorm(boundaries=bn, ncolors=cmp.N)
        # #fig.subplots_adjust(top=0.9)
        #
        # p0 = ax.get_position().get_points().flatten()
        # # l, b, w = [0.15, 0.92, 0.7]
        # # cax = fig.add_axes([0.15, 0.9, 0.2, 0.05])
        # # cax = fig.add_axes([p0[0], 1, p0[2]-p0[0], 0.05])
        # # cb_ax = mpl_toolkits.axes_grid1.inset_locator.inset_axes(ax, loc=3)
        # cb =fig.colorbar(mappable=mpl.cm.ScalarMappable(norm=norm, cmap=cmp), use_gridspec=False,
        #                  ax=ax, orientation='horizontal',shrink=1, fraction=0.1, pad=0.1)

        if iax == 0:
            if options.unitname.lower() == "mev":
                ax.set_ylabel("Frequency (meV)", labelpad=10)
            elif options.unitname.lower() == "thz":
                ax.set_ylabel("Frequency (THz)", labelpad=10, fontsize=12)
            else:
                ax.set_ylabel("Frequency (cm${}^{-1}$)", labelpad=10)

        else:
            ax.set_yticklabels([])
            ax.set_yticks([])

        plt.axis([xmin_ax[iax], xmax_ax[iax], ymin, ymax])
        ax.set_xticks(xticks_ax[iax])
        ax.set_xticklabels(xticklabels_ax[iax])
        ax.xaxis.grid(True, linestyle='-')
        ax.tick_params(axis='y', direction='in', labelsize=12)
        # ax.minorticks_on()
        ax.tick_params(axis='both', direction='in', right=True,top=True)
        # ax.tick_params(axis='y',direction='in',right=True,which='minor')

        # ax.set_yticklabels(labels = ax.get_yticklabels(), fontsize=8)

        # if options.print_key and iax == 0:
        #     ax.legend(loc='best', prop={'size': 8})



def plot_band_dos(fig_name):
    fig=plt.figure(figsize=(6,4))
    gs = GridSpec(nrows=2, ncols=3, width_ratios=[3,1,1], height_ratios=[0.04,0.96])
    gs.update(wspace=0.1, hspace=0.05)
    ax0 = plt.subplot(gs[0,0])
    # ax_0.set_axis_off()
    cmp = mpl.colors.ListedColormap(colors = color)
    norm = mpl.colors.BoundaryNorm(boundaries=bn, ncolors=cmp.N)
    # fig.subplots_adjust(top=0.9)

    # l, b, w = [0.15, 0.92, 0.7]
    # cax = fig.add_axes([0.15, 0.9, 0.2, 0.05])
    # cax = fig.add_axes([p0[0], 1, p0[2]-p0[0], 0.05])
    # cb_ax = mpl_toolkits.axes_grid1.inset_locator.inset_axes(ax, loc=3)
    cb = mpl.colorbar.ColorbarBase(ax0, cmap=cmp, norm=norm, orientation='horizontal', ticklocation='top',extend='neither')

    cb.ax.tick_params(labelsize=10, pad=-0.5)



    ax1 = plt.subplot(gs[1,0])
    nax, xticks_ax, xticklabels_ax, xmin_ax, xmax_ax, ymin, ymax, \
        data_merged_ax = preprocess_data(files, options.unitname, options.normalize_xaxis)
    plot_band(ax1, nax, xticks_ax, xticklabels_ax,
             xmin_ax, xmax_ax, ymin, ymax, data_merged_ax)

#color = np.array([(21,9,242), (64,26,217), (106,43,191),(149,60,166),(191,77,140),(234,94,115)])/255
    ax2 = plt.subplot(gs[1,1])
    plot_dos(ax2, dos_file_1['dos_path'], dos_file_1['dos_temp'], np.array((21,9,242))/255, ymax, ymin)
    ax3 = plt.subplot(gs[1,2])
    plot_dos(ax3, dos_file_2['dos_path'], dos_file_2['dos_temp'],np.array((191,77,140))/255, ymax, ymin)

    # plt.tight_layout()
    if fig_name:
        plt.savefig(fig_name, dpi=400)
    plt.show()

def plot_band_dos_2(fig_name):
    fig=plt.figure(figsize=(6,4))
    gs = GridSpec(nrows=2, ncols=3, width_ratios=[3,1,1], height_ratios=[0.04,0.96])
    gs.update(wspace=0.1, hspace=0.05)
    ax0 = plt.subplot(gs[0,0])
    # ax_0.set_axis_off()
    color = np.array(
        [(21, 9, 242), (64, 26, 217), (106, 43, 191), (149, 60, 166), (191, 77, 140), (234, 94, 115)]) / 255
    cmp = mpl.colors.ListedColormap(colors = color)
    norm = mpl.colors.BoundaryNorm(boundaries=bn, ncolors=cmp.N)
    # fig.subplots_adjust(top=0.9)

    # l, b, w = [0.15, 0.92, 0.7]
    # cax = fig.add_axes([0.15, 0.9, 0.2, 0.05])
    # cax = fig.add_axes([p0[0], 1, p0[2]-p0[0], 0.05])
    # cb_ax = mpl_toolkits.axes_grid1.inset_locator.inset_axes(ax, loc=3)
    cb = mpl.colorbar.ColorbarBase(ax0, cmap=cmp, norm=norm, orientation='horizontal', ticklocation='top',extend='neither')
    cb.ax.tick_params(labelsize=10, pad=-0.5)



    ax1 = plt.subplot(gs[1,0])
    nax, xticks_ax, xticklabels_ax, xmin_ax, xmax_ax, ymin, ymax, \
        data_merged_ax = preprocess_data(files, options.unitname, options.normalize_xaxis)
    plot_band(ax1, nax, xticks_ax, xticklabels_ax,
             xmin_ax, xmax_ax, ymin, ymax, data_merged_ax)

#color = np.array([(21,9,242), (64,26,217), (106,43,191),(149,60,166),(191,77,140),(234,94,115)])/255
    ax2 = plt.subplot(gs[1,1])
    plot_dos(ax2, dos_file_1['dos_path'], dos_file_1['dos_temp'], np.array((21,9,242))/255, ymax, ymin)
    ax3 = plt.subplot(gs[1,2])
    plot_dos(ax3, dos_file_2['dos_path'], dos_file_2['dos_temp'],np.array((191,77,140))/255, ymax, ymin)

    # plt.tight_layout()
    if fig_name:
        plt.savefig(fig_name, dpi=600)
    plt.show()


if __name__ == '__main__':
    '''
    Simple script for visualizing phonon dispersion relations.
    Usage:
    $ python plot_band.py [options] file1.bands file2.bands ...
    For details of available options, please type
    $ python plot_band.py -h
    '''
    """
    Test arguments:
    ./BrCsSn/scfph/300/pc_disp.bands
    ./BrCsSn/scfph/400/pc_disp.bands
    """
    # fakeArgs = ['-u', 'Thz',  './BrCsSn/scfph/kapa_100/pc_disp.bands','./BrCsSn/scfph/kapa_150/pc_disp.bands'
    #             ,'./BrCsSn/scfph/kapa_200/pc_disp.bands','./BrCsSn/scfph/kapa_250/pc_disp.bands',
    #            './BrCsSn/scfph/kapa_300/pc_disp.bands', './BrCsSn/scfph/kapa_350/pc_disp.bands']

    # fakeArgs = ['-u', 'Thz']
    fakeArgs = ['-u', 'mev']
    start_temp=150
    stop_temp =450
    temp_range = range(start_temp, stop_temp, 50)

    ########---------------BrCsSn------------------------------#############
    dos_file_1 = {'dos_path': './BrCsSn/scfph/kapa_100/pc_dos.dos', 'dos_temp':'T=100K'}
    dos_file_2 = {'dos_path': './BrCsSn/scfph/kapa_300/pc_dos.dos', 'dos_temp': 'T=300K'}
    base_dir = './BrCsSn/scfph/'
    fig_name = './BrCsSn_bands_mev.png'

    # ########-----------------------BrCdCs-----------------------------############
    # dos_file_1 = {'dos_path': './Br3Cd1Cs1/scfph2/kapa_100/pc_dos.dos', 'dos_temp':'T=100K'}
    # dos_file_2 = {'dos_path': './Br3Cd1Cs1/scfph2/kapa_300/pc_dos.dos', 'dos_temp': 'T=300K'}
    # base_dir = './Br3Cd1Cs1/scfph2/'
    # fig_name = './BrCdCs_bands_mev.eps'

    # #########------------------------BrCaCs----------------------------############
    # dos_file_1 = {'dos_path': './Br3Ca1Cs1/scfph/kapa_100/pc_dos.dos', 'dos_temp':'T=100K'}
    # ##### dos_file_1 = {'dos_path': './Br3Ca1Cs1/phonons/pc_dos.dos', 'dos_temp':'T=100K'}
    # dos_file_2 = {'dos_path': './Br3Ca1Cs1/scfph/kapa_300/pc_dos.dos', 'dos_temp': 'T=300K'}
    # base_dir = './Br3Ca1Cs1/scfph/'
    # fig_name = './BrCaCs_bands_mev.png'


    file_dir = [base_dir+'kapa_{}'.format(i)+'/pc_disp.bands' for i in temp_range]
    for i in file_dir:
        fakeArgs.append(i)
    # fakeArgs.append('./BrCsSn/phonons/pc_disp.bands')
    print(fakeArgs)
    print(type(temp_range))
    bn=[x for x in temp_range]
    bn.insert(0, start_temp - 50)

    # fakeArgs.append('./Br3Ca1Cs1/phonons/pc_disp.bands')


    print(bn)

    options, args = parser.parse_args(fakeArgs)
    files = args[0:]
    nfiles = len(files)

    if nfiles == 0:
        print("Usage: plotband.py [options] file1.bands file2.bands ...")
        print("For details of available options, please type\n$ python plotband.py -h")
        exit(1)
    else:
        print("Number of files = %d" % nfiles)
    plot_band_dos(fig_name)



