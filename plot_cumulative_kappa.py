import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import FormatStrFormatter

mpl.rc('xtick', labelsize=14)
mpl.rc('ytick', labelsize=14)
mpl.rc('axes', labelsize=14)
mpl.rc('font', **{'family': 'Times New Roman', 'sans-serif': ['Helvetica']})
mpl.rcParams.update({"mathtext.fontset":'stix',"font.serif": ['SimSun'],})

def read_cumulative_kappa(file_name):
    f = open(file_name)
    frequency=[]
    kl_spec=[]
    kl_cumulative_spec=[]
    for line in f.readlines():
        list_from_line = line.strip().split()
        frequency.append(float(list_from_line[0]))
        kl_spec.append(float(list_from_line[1]))
        kl_cumulative_spec.append(float(list_from_line[4]))
    return np.array(frequency), np.array(kl_spec), np.array(kl_cumulative_spec)

def kayser_to_mev(array):
    kayser_to_mev = 0.0299792458 * 1.0e+12 * \
                    6.62606896e-34 / 1.602176565e-19 * 1000

    for i in range(len(array)):
        array[i] *= kayser_to_mev
    return array

def plot_cumulative_kappa(data_name_1, data_name_2, fig_name, text_name):
    fig, ax1 = plt.subplots()
    fig.set_size_inches(8,4)
    frequency, kl, kl_cumulative = read_cumulative_kappa(data_name_1)
    frequency2, kl2, kl_cumulative2 = read_cumulative_kappa(data_name_2)
    color1='mediumslateblue'
    xmax = 23 # 32 for BrCaCs/BrCdCs
    # ymax1=3
    # ymax2 = 1.2  # 1.5 for BrCaCs; 0.8 for BrCdCs

    # The SCP + BTE method
    ######################BrCaCS###################
    # ymax1=5
    # ymax2 = 2.5
    # x1_text, y1_text, x2_text, y2_text = [0.7, 0.4, 0.74, 0.8]

    ######################BrCdCS###################
    ymax1=3
    ymax2= 1.0
    x1_text, y1_text, x2_text, y2_text= [0.7, 0.45, 0.74, 0.8]

    ######################BrCSSn###################
    ymax1=4.5
    ymax2= 2.0
    x1_text, y1_text, x2_text, y2_text= [0.6, 0.45, 0.64, 0.8]

    x1=kayser_to_mev(frequency)
    x2 = kayser_to_mev(frequency2)
    ax1.plot(x1, kl*100, linestyle='-', c='blueviolet')
    ax1.fill(x1, kl*100, c='blueviolet', alpha=0.25)
    ax1.plot(x2, kl2*100, linestyle='-', c='sandybrown')
    ax1.fill(x2, kl2*100, c='sandybrown', alpha=0.25)
    ax1.set_xlabel('Frequency (meV)')
    ax1.set_ylabel('$\kappa_L(\Omega)$ [10$^{-2}$ W/mK/cm$^{-1}$')
    ax1.set_ylim(bottom=0, top=ymax1)
    ax1.set_xlim(left=0, right=xmax)
    ax1.tick_params(axis='both', which='both', direction='in', top=True, right=True)

    ax2 = ax1.twinx()
    ax2.plot(x1, kl_cumulative, linestyle='--', c='blueviolet')
    ax2.plot(x2, kl_cumulative2, linestyle ='--', c='sandybrown')

    ax2.set_ylim(bottom=0, top= ymax2)
    ax2.set_ylabel('Cumulative $\kappa_L $ (W/mK)')
    ax2.tick_params(axis='both', which='both', direction='in', top=True, right=True)
    plt.text(x= 0.1, y=0.87, s=text_name, fontsize = 16, transform=ax2.transAxes )

    # The SCP + BTE method
    plt.text(x= x1_text, y=y1_text, s='300K', fontsize = 16, color= 'blueviolet', transform=ax2.transAxes )
    plt.text(x= x2_text, y=y2_text, s='100K', fontsize = 16, color='sandybrown', transform=ax2.transAxes )

    plt.tight_layout()
    if fig_name:
        plt.savefig(fig_name, dpi=400)
    plt.show()

def plot_figure5():
    fig = plt.figure(figsize=(8, 9))
    gs = GridSpec(nrows=3, ncols=1)
    gs.update(wspace=0.3)

    data_name_1_a = './Br3Ca1Cs1/scfph/kapa_300/kappa_spec_culm.dat'
    data_name_2_a = './Br3Ca1Cs1/scfph/kapa_100/kappa_spec_culm.dat'
    ax1 = plt.subplot(gs[0])
    frequency1_a, kl1_a, kl_cumulative1_a = read_cumulative_kappa(data_name_1_a)
    frequency2_a, kl2_a, kl_cumulative2_a = read_cumulative_kappa(data_name_2_a)

    data_name_1_b = './BrCsSn/scfph/kapa_300/kappa_spec_culm.dat'
    data_name_2_b = './BrCsSn/phonons/kappa_spec_culm.dat'
    ax2 = plt.subplot(gs[1])
    frequency1_b, kl1_b, kl_cumulative1_b = read_cumulative_kappa(data_name_1_b)
    frequency2_b, kl2_b, kl_cumulative2_b = read_cumulative_kappa(data_name_2_b)

    data_name_1_c = './Br3Cd1Cs1/scfph2/kapa_300/kappa_spec_culm.dat'
    data_name_2_c = './Br3Cd1Cs1/phonons/kappa_spec_culm.dat'
    ax3 = plt.subplot(gs[2])
    frequency1_c, kl1_c, kl_cumulative1_c = read_cumulative_kappa(data_name_1_c)
    frequency2_c, kl2_c, kl_cumulative2_c = read_cumulative_kappa(data_name_2_c)



    plot_single_figure(ax1, frequency1=frequency1_a , kl=kl1_a, frequency2=frequency2_a,
                       kl2=kl2_a, xmax=32, ymax1=5, ymax2=2.5, kl_cumulative1=kl_cumulative1_a,
                       kl_cumulative2=kl_cumulative2_a,
                       x1_text=0.7, y1_text=0.4, x2_text=0.74, y2_text=0.8, s1='300 K', s2='100 K',
                       text_name='CsCaBr3 (SCP+BTE)')

    plot_single_figure(ax2, frequency1=frequency1_b , kl=kl1_b, frequency2=frequency2_b,
                       kl2=kl2_b, xmax=23, ymax1=3, ymax2=1.2, kl_cumulative1=kl_cumulative1_b,
                       kl_cumulative2=kl_cumulative2_b,
                       x1_text=0.6, y1_text=0.78, x2_text=0.6, y2_text=0.5, s1='SCP + BTE', s2='HP + BTE',
                       text_name='CsSnBr3 (300 K)')

    plot_single_figure(ax3, frequency1=frequency1_c , kl=kl1_c, frequency2=frequency2_c,
                       kl2=kl2_c, xmax=23, ymax1=3, ymax2=0.8, kl_cumulative1=kl_cumulative1_c,
                       kl_cumulative2=kl_cumulative2_c,
                       x1_text=0.7, y1_text=0.78, x2_text=0.6, y2_text=0.2, s1='SCP + BTE', s2='HP + BTE',
                       text_name='CsCdBr3 (300 K)', xlabel='Frequency (meV)')

    plt.tight_layout()
    fig_name='./figure5.eps'
    if fig_name:
        plt.savefig(fig_name, dpi=500)
    plt.show()


def plot_single_figure(ax1, frequency1, kl, frequency2, kl2, xmax, ymax1, ymax2, kl_cumulative1, kl_cumulative2,
                       x1_text, y1_text, x2_text, y2_text, s1, s2, text_name, xlabel=None):
    """
    :param ax1: The axes class
    :param frequency1: x axis
    :param kl: the lattice thermal conductivity
    :param frequency2: x axis of the second file for comparison
    :param kl2: the lattice thermal conductivity for comparison
    :param xmax: set the maximum of x-axis
    :param ymax1: set the maximum of the value of the kappa_L
    :param ymax2: set the maximum of the value of the cumulative kappa_L
    :param kl_cumulative1:
    :param kl_cumulative2:
    :param x1_text: text positon of the first file in x axis (temperature or methods)
    :param y1_text: text positon of the first file in y axis
    :param x2_text: text positon of the second file in x axis
    :param y2_text: text positon of the second file in y axis
    :param s1: string. the text name of the first file
    :param s2: string. the text name of the second file
    :param text_name: string. the figure notation.

    :return:
    """
    x1 = kayser_to_mev(frequency1)
    x2 = kayser_to_mev(frequency2)
    fs = 14
    ax1.plot(x1, kl * 100, linestyle='-', c='blueviolet')
    ax1.fill(x1, kl * 100, c='blueviolet', alpha=0.25)
    ax1.plot(x2, kl2 * 100, linestyle='-', c='sandybrown')
    ax1.fill(x2, kl2 * 100, c='sandybrown', alpha=0.25)
    if xlabel:
        ax1.set_xlabel(xlabel)
    ax1.set_ylabel('$\kappa_L(\Omega)$ [10$^{-2}$ W/mK/cm$^{-1}$')
    ax1.set_ylim(bottom=0, top=ymax1)
    ax1.set_xlim(left=0, right=xmax)
    ax1.tick_params(axis='both', which='both', direction='in', top=True, right=True)

    ax2 = ax1.twinx()
    ax2.plot(x1, kl_cumulative1, linestyle='--', c='blueviolet')
    ax2.plot(x2, kl_cumulative2, linestyle='--', c='sandybrown')

    ax2.set_ylim(bottom=0, top=ymax2)
    ax2.set_ylabel('Cumulative $\kappa_L $ (W/mK)')
    xmajor_formatter = FormatStrFormatter('%.1f')
    ax2.yaxis.set_major_formatter(xmajor_formatter)
    ax2.tick_params(axis='both', which='both', direction='in', top=True, right=True)
    ax2.text(x=0.1, y=0.87, s=text_name, fontsize=fs, transform=ax2.transAxes)

    # The SCP + BTE method
    ax2.text(x=x1_text, y=y1_text, s=s1, fontsize=fs, color='blueviolet', transform=ax2.transAxes)
    ax2.text(x=x2_text, y=y2_text, s=s2, fontsize=fs, color='sandybrown', transform=ax2.transAxes)


def plot_materials(mat):
    if mat=='BrCaCs':
        data_name_1 = './Br3Ca1Cs1/scfph/kapa_300/kappa_spec_culm.dat'
        data_name_2 = './Br3Ca1Cs1/phonons/kappa_spec_culm.dat'
        fig_name = 'kappa_spec_BrCaCs.png'
        text_name = 'CsCaBr3 (300K)'
        plot_cumulative_kappa(data_name_1, data_name_2, fig_name, text_name)
    elif mat=='BrCdCs':
        data_name_1 = './Br3Cd1Cs1/scfph2/kapa_300/kappa_spec_culm.dat'
        data_name_2 = './Br3Cd1Cs1/phonons/kappa_spec_culm.dat'
        fig_name = 'kappa_spec_BrCdCs.png'
        text_name = 'CsCdBr3 (300K)'
        plot_cumulative_kappa(data_name_1, data_name_2, fig_name, text_name)
    elif mat=='BrCsSn':
        data_name_1 = './BrCsSn/scfph/kapa_300/kappa_spec_culm.dat'
        data_name_2 = './BrCsSn/phonons/kappa_spec_culm.dat'
        fig_name = 'kappa_spec_BrCsSn.png'
        text_name = 'CsSnBr3 (300K)'
        plot_cumulative_kappa(data_name_1, data_name_2, fig_name, text_name)

if __name__=='__main__':

    # data_name_1 = './Br3Ca1Cs1/scfph/kapa_300/kappa_spec_culm.dat'
    # data_name_2 = './Br3Ca1Cs1/scfph/kapa_100/kappa_spec_culm.dat'
    # fig_name = 'kappa_spec_BrCaCs_temp.png'
    # text_name = 'CsCaBr3 (SCP+BTE)'
    # plot_cumulative_kappa(data_name_1, data_name_2, fig_name, text_name)

    # data_name_1 = './Br3Cd1Cs1/scfph2/kapa_300/kappa_spec_culm.dat'
    # data_name_2 = './Br3Cd1Cs1/scfph2/kapa_100/kappa_spec_culm.dat'
    # fig_name = 'kappa_spec_BrCdCs_temp.png'
    # text_name = 'CsCdBr3 (SCP+BTE)'
    # plot_cumulative_kappa(data_name_1, data_name_2, fig_name, text_name)


    # data_name_1 = './BrCsSn/scfph/kapa_300/kappa_spec_culm.dat'
    # data_name_2 = './BrCsSn/scfph/kapa_100/kappa_spec_culm.dat'
    # fig_name = 'kappa_spec_BrCsSn_temp.png'
    # text_name = 'CsSnBr3 (SCP+BTE)'
    # plot_cumulative_kappa(data_name_1, data_name_2, fig_name, text_name)
    plot_figure5()