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
import matplotlib as mpl
from matplotlib import colors
import mpl_toolkits
from matplotlib.gridspec import GridSpec
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
mpl.rc('font', **{'family': 'Times New Roman', 'sans-serif': ['Helvetica']})
mpl.rc('xtick', labelsize=12)
mpl.rc('ytick', labelsize=16)
mpl.rc('axes', labelsize=16)
mpl.rc('lines', linewidth=1.5)
mpl.rc('legend', fontsize='small')
mpl.rcParams.update({"mathtext.fontset":'stix',"font.serif": ['SimSun'],})
# line colors and styles
color = np.array([(144,53,235), (255,26,26), (82,82,82)])/255
# 紫色，红色。灰色
lsty = ['-', '-', '-', '-', '--', '--', '--', '--']
bn = [0, 300, 400]

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


def read_pr(file):
    f = open(file)
    natoms = 5  # number of atoms in the unit cell!!!!!!!!!!!!!!!!!!!!
    nlines = 15 #  number of lines in the phonon spectra!!!!!!!!!!!!
    nk=250 # number of k points !!!!!!!!!!!!!!!!!!!
    i = -1
    notations = f.readline()
    line = f.readline()
    atom1 = np.zeros(shape=(nk, nlines))
    for line in f.readlines():
        list_from_line = line.strip().split()
        if list_from_line==[]:
            pass
        elif list_from_line[0] =='#':
            i = i + 1
            j1 = 0
        else:
            pr = float(list_from_line[2])
            atom1[i,j1]=pr
            j1 = j1 + 1
    f.close()
    return atom1


def read_apr(file):
    f = open(file)
    natoms = 5  # number of atoms in the unit cell!!!!!!!!!!!!!!!!!!!!
    nlines = 15 #  number of lines in the phonon spectra!!!!!!!!!!!!
    nk=750 # number of k points !!!!!!!!!!!!!!!!!!!
    i = -1
    notations = f.readline()
    line = f.readline()
    atom1 = np.zeros(shape=(nk, nlines))
    atom2 = np.zeros(shape=(nk, nlines))
    atom3 = np.zeros(shape=(nk, nlines))
    atom4 = np.zeros(shape=(nk, nlines))
    atom5 = np.zeros(shape=(nk, nlines))


    for line in f.readlines():
        list_from_line = line.strip().split()
        if list_from_line==[]:
            pass
        elif list_from_line[0] =='#':
            i = i + 1
            j1,j2,j3,j4,j5 = [0,0,0,0,0]
        else:
            atom_i = list_from_line[2]
            apr = float(list_from_line[3])
            if atom_i=='1':
                atom1[i,j1]=apr
                j1 = j1 + 1
            elif atom_i=='2':
                atom2[i,j2]=apr
                j2 = j2 + 1
            elif atom_i=='3':
                atom3[i,j3]=apr
                j3 = j3 + 1
            elif atom_i=='4':
                atom4[i,j4]=apr
                j4 = j4 + 1
            elif atom_i=='5':
                atom5[i,j5]=apr
                j5 = j5 + 1
            else:
                print('wrong information!')
    f.close()
    return atom1, atom2, atom3, atom4, atom5


def plot_band_minus(fig_name, atom, cb_string):
    # fig = plt.figure(figsize=(6, 4))
    # gs = GridSpec(nrows=1, ncols=3, width_ratios=[1, 1, 1])
    # gs.update(wspace=0.1)
    # print(files)
    # ax_i = plt.subplot(gs[0, i])
    font = {'family': 'sans-serif',
            'weight': 'normal',
            'size': 16,
            }
    nax, xticks_ax, xticklabels_ax, xmin_ax, xmax_ax, ymin, ymax, \
    data_merged_ax = preprocess_data(files, options.unitname, options.normalize_xaxis)
    fig = plt.figure()

    width_ratios = []
    used_colors = []
    for xmin, xmax in zip(xmin_ax, xmax_ax):
        width_ratios.append(xmax - xmin)

    gs = GridSpec(nrows=1, ncols=nax, width_ratios=width_ratios)
    gs.update(wspace=0.1)
    colors_list = np.array([(207,224,94),(155,216,103), (115,204,113), (80,192,112),(57,175,131),(39,155,136), (27,131,139),(23,99,135),(28,23,110)])/255
    cmap= colors.ListedColormap(name ='own2', colors=colors_list)
    norm = colors.Normalize(vmin=0, vmax=1)
    for iax in range(nax):
        ax = plt.subplot(gs[iax])
        # ax.plot(np.arange(start=xmin, stop=xmax,step=(xmax-xmin)/100), np.zeros(shape=(100,)), linestyle='-', color='grey')
        ax.axhline(y=0, xmin=xmin, xmax=xmax, linestyle='--', color='grey')

        for i in range(len(data_merged_ax[iax])):
            if len(data_merged_ax[iax][i]) > 0:
                ax.scatter(data_merged_ax[iax][i][0:, 0], data_merged_ax[iax][i][0:, 1],
                        marker='s', c=atom[:, i], s=1, alpha=0.8, cmap=cmap, norm=norm)
                used_colors.append(color[i])

                for j in range(2, len(data_merged_ax[iax][i][0][0:])):
                    ax.scatter(data_merged_ax[iax][i][0:, 0], data_merged_ax[iax][i][0:, j],
                            marker='s', alpha =0.8,  c=atom[:,j-1], s=1, cmap=cmap, norm=norm)

        if iax == 0:
            if options.unitname.lower() == "mev":
                ax.set_ylabel("Frequency (meV) ", labelpad=20)
            elif options.unitname.lower() == "thz":
                ax.set_ylabel("Frequency (THz)", labelpad=20)
            else:
                ax.set_ylabel("Frequency (cm${}^{-1}$)", labelpad=10)

        else:
            ax.set_yticklabels([])
            ax.set_yticks([])
        sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
        cb = plt.colorbar(sm, drawedges=False, shrink=0.9, pad=0.03)
        tick_len = 0
        tick_size=12
        cb.ax.tick_params(labelsize=tick_size, length=tick_len)
        cb_text_string = cb_string
        cb.ax.text(x=0, y=1.02, s=cb_text_string, fontsize=12, fontweight='normal')
        cb.outline.set_visible(False)
        plt.axis([xmin_ax[iax], xmax_ax[iax], ymin, ymax])
        ax.set_xticks(xticks_ax[iax])
        ax.set_xticklabels(xticklabels_ax[iax], fontdict=font)

        ax.xaxis.grid(True, linestyle='-')
        ax.tick_params(axis='y', direction='in')
        # ax.minorticks_on()
        ax.tick_params(axis='both', direction='in', right=True,top=True)

    plt.tight_layout()
    if fig_name:
        plt.savefig(fig_name, dpi=400)
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

    # #################### BrCaCS#########################################
    fakeArgs = ['-u', 'meV', './Br3Ca1Cs1/phonons/pc_disp.bands']
    atom1, atom2, atom3, atom4, atom5 = read_apr('./Br3Ca1Cs1/phonons/pc_disp.apr')
    fig_name = 'phonons_minus_BrCaCs_Br.eps'

    ####################BrCdCs#####################################
    # fakeArgs = ['-u', 'meV', './Br3Cd1Cs1/phonons/pc_disp.bands']
    # atom1, atom2, atom3, atom4, atom5 = read_apr('./Br3Cd1Cs1/phonons/pc_disp.apr')
    # fig_name = 'phonons_minus_BrCdCs_Cd.png'

    # ####################BrCsSn#####################################
    # fakeArgs = ['-u', 'meV', './BrCsSn/phonons/pc_disp.bands']
    # atom1, atom2, atom3, atom4, atom5 = read_apr('./BrCsSn/phonons/pc_disp.apr')
    # fig_name = 'phonons_minus_BrCsSn_Br.eps'
    #

    atom =atom1+atom2+atom3
    # atom = atom4
    # atom=atom5
    cb_string='Br'

    options, args = parser.parse_args(fakeArgs)
    files = args[0:]
    nfiles = len(files)

    if nfiles == 0:
        print("Usage: plotband.py [options] file1.bands file2.bands ...")
        print("For details of available options, please type\n$ python plotband.py -h")
        exit(1)
    else:
        print("Number of files = %d" % nfiles)

    plot_band_minus(fig_name, atom, cb_string)

