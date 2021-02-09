import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import os

mpl.rc('font', **{'family': 'Times New Roman','weight':'normal', 'sans-serif': ['Helvetica']})
mpl.rc('xtick', labelsize=12)
mpl.rc('ytick', labelsize=12)
mpl.rc('axes', labelsize=12)
mpl.rc('lines', linewidth=1.5)
mpl.rc('legend', fontsize='small')
mpl.rcParams.update({"mathtext.fontset":'stix',"font.serif": ['SimSun'],})



kayser_to_mev = 0.0299792458 * 1.0e+12 * \
                6.62606896e-34 / 1.602176565e-19 * 1000


def plot_sr_temp(file_name_low_T, file_name_high_T, file_name_max_T):
    tau_low = np.loadtxt(file_name_low_T)
    tau_high = np.loadtxt(file_name_high_T)
    tau_max = np.loadtxt(file_name_max_T)

    # fig = plt.figure(figsize=(8,8), dpi=300)
    start_index = 3
    markersize= 7
    labelsize = 16

    plt.scatter(tau_low[start_index:,2] * kayser_to_mev, 1/tau_low[start_index:,3], c='lightskyblue', s=markersize, marker='o')
    plt.scatter(tau_high[start_index:,2] * kayser_to_mev, 1/tau_high[start_index:,3], c='#ff474c', s=markersize, marker='o')
    plt.scatter(tau_max[start_index:,2] * kayser_to_mev, 1/tau_max[start_index:, 3], c='purple', s=markersize, marker='o')
    x=np.arange(0, 30, 0.1)
    plt.plot(x, x, color='k', linestyle='-', linewidth=1.5)
    plt.yscale('log')
    plt.tick_params(axis='both', which='both', direction='in', right=True, top=True)
    plt.grid(b=True, which='major', axis='both', linestyle='--')

    plt.annotate(r'$\Gamma$ (100 K)', color='lightskyblue', xy=(16, 0.12), xytext=(10.1, 0.031),
                 arrowprops=dict(facecolor='lightskyblue', width=1, headwidth=4.5,
                                 edgecolor='lightskyblue'),ha='center', va='center', fontsize=labelsize)

    plt.annotate(r'$\Gamma$ (300 K)', color='#ff474c', xy=(19.3, 0.57), xytext=(27.3, 0.058),
                 arrowprops=dict(facecolor='#ff474c', width=1, headwidth=4.5,
                                 edgecolor='#ff474c'),ha='right', va='top', fontsize=labelsize)

    plt.annotate(r'$\Gamma$ (600 K)', color='purple', xy=(18.9, 2.63), xytext=(13.5, 5.5),
                 arrowprops=dict(facecolor='purple', width=1, headwidth=4.5,
                                 edgecolor='purple'),ha='center', va='center', fontsize=labelsize)

    # plt.axis(xmin=0, xmax=30, ymin=0, ymax=100)
    plt.xlim((0, 30))

    plt.xlabel('Energy (meV)', fontsize=labelsize)
    plt.ylabel('Scattering rates (ps$^{-1}$)', fontsize=labelsize)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    # plt.gca().set_aspect(1.0 / plt.gca().get_data_ratio(), adjustable='box')
    # plt.tight_layout()
    plt.savefig('./pictures/scattering_rate_1.png', dpi=400)
    plt.show()


def plot_sr_methods(file_path_1, file_path_2, fig_name, **kwargs):

    """

    :param file_path_1: str.
    :param file_path_2: str.
    :param fig_name: str.
    :param kwargs:
    xmax: the max value of the x-axis
    ymax: the max value of the y-axis

    :return:
    """
    allowed_kwargs={'xmin', 'xmax', 'ymin', 'ymax'}
    xmax= kwargs.get('xmax')
    ymin = kwargs.get('ymin')
    tau_1 = np.loadtxt(file_path_1)
    tau_2 = np.loadtxt(file_path_2)

    # fig = plt.figure(figsize=(8,8), dpi=300)
    start_index = 3
    markersize= 7
    labelsize = 16

    plt.scatter(tau_1[start_index:,2] * kayser_to_mev, 1/tau_1[start_index:,3], c='lightskyblue', s=markersize, marker='o')
    plt.scatter(tau_2[start_index:,2] * kayser_to_mev, 1/tau_2[start_index:,3], c='#ff474c', s=markersize, marker='o')
    x=np.arange(0, 30, 0.1)
    plt.plot(x, 2*x, color='k', linestyle='-', linewidth=1.5)
    plt.yscale('log')
    plt.tick_params(axis='both', which='both', direction='in', right=True, top=True)
    plt.grid(b=True, which='major', axis='both', linestyle='--')

    plt.annotate(r'HA + 3ph', color='lightskyblue', xy=(5.4, 14), xytext=(11.5, 24),
                 arrowprops=dict(facecolor='lightskyblue', width=1, headwidth=4.5,
                                 edgecolor='lightskyblue'),ha='center', va='center', fontsize=labelsize)

    plt.annotate(r'SCP + 3ph', color='#ff474c', xy=(10.3, 0.57), xytext=(16.3, 0.058),
                 arrowprops=dict(facecolor='#ff474c', width=1, headwidth=4.5,
                                 edgecolor='#ff474c'),ha='center', va='top', fontsize=labelsize)


    # plt.axis(xmin=0, xmax=30, ymin=0, ymax=100)
    plt.xlim((0, xmax))
    plt.ylim(bottom=ymin)

    plt.xlabel('Energy (meV)', fontsize=labelsize)
    plt.ylabel('Scattering rates (ps$^{-1}$)', fontsize=labelsize)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    # plt.gca().set_aspect(1.0 / plt.gca().get_data_ratio(), adjustable='box')
    # plt.tight_layout()
    fname = os.path.join('./pictures', fig_name)
    plt.savefig(fname, dpi=400)
    plt.show()

if __name__=='__main__':
    # tau_Br3CaCs_100 = './Br3Ca1Cs1/scfph/kapa_300/tau100K_iso.dat'
    # tau_Br3CaCs_300 = './Br3Ca1Cs1/scfph/kapa_300/tau300K_iso.dat'
    # tau_Br3CaCs_600 = './Br3Ca1Cs1/scfph/kapa_300/tau600K_iso.dat'
    # plot_sr_temp(tau_Br3CaCs_100, tau_Br3CaCs_300, tau_Br3CaCs_600)

    # tau_1 = './Br3Cd1Cs1/phonons/tau300K_iso.dat'
    # tau_2 = './Br3Cd1Cs1/scfph2/kapa_300/tau300K_iso.dat'
    # plot_sr_methods(tau_1, tau_2, fig_name='scattering_rate_BrCdCs.png', xmax=25, ymin=0.01)


    tau_1 = './BrCsSn/phonons/tau300K_iso.dat'
    tau_2 = './BrCsSn/scfph/kapa_300/tau300K_iso.dat'
    plot_sr_methods(tau_1, tau_2, fig_name='scattering_rate_BrCsSn.png', xmax=20, ymin=0.01)

    # tau_1 = './Br3Cd1Cs1/phonons/tau300K_iso.dat'
    # tau_2 = './BrCsSn/scfph/kapa_300/tau300K_iso.dat'
    # plot_sr_methods(tau_1, tau_2, fig_name='scattering_rate_BrCsSn.png', xmax=20, ymin=0.01)





