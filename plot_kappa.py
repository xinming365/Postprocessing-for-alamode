import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


def plot_kappa(file1, file2, file3, file4, file5, file6):
    f1 = np.loadtxt(file1)
    f2 = np.loadtxt(file2)
    f3 = np.loadtxt(file3)
    f4 = np.loadtxt(file4)
    f5 = np.loadtxt(file5)
    f6 = np.loadtxt(file6)


    plt.plot(f1[:,0], f1[:,1], linestyle='--', marker='s', color='gray', markerfacecolor='white',
             label="CsCaBr$_3$, HP+BTE")
    plt.plot(f2[:, 0], f2[:, 1], linestyle='--', marker='s', color='gray', label="CsCaBr$_3$, SCP+BTE")

    plt.plot(f3[:,0], f3[:,1], linestyle='--', marker='o', color='darkviolet', markerfacecolor='white',
             label = "CsCdBr$_3$, HP+BTE")
    plt.plot(f4[:, 0], f4[:, 1], linestyle='--', marker='o', color='darkviolet', label = "CsCdBr$_3$, SCP+BTE")

    plt.plot(f5[:,0], f5[:,1], linestyle='--', marker='^', color='red', markerfacecolor='white',
             label = "CsSnBr$_3$, HP+BTE")
    plt.plot(f6[:, 0], f6[:, 1], linestyle='--', marker='^', color='red', label="CsSnBr$_3$, SCP+BTE")

    plt.legend()
    labelsize=13
    ticksize = 13
    plt.xlabel('Temperature (K)', fontsize=labelsize)
    plt.ylabel('$\kappa_L$ (W/mK)', fontsize=labelsize)
    plt.tick_params(axis='both', which='both', direction='in', right=True, top=True)
    plt.xticks(fontsize=ticksize)
    plt.yticks(fontsize=ticksize)
    plt.savefig('./pictures/kappa.png', dpi=400)
    plt.show()


def fit_curve(xdata, ydata, c, data_label):
    def func(x, a, b, c):
        return a * np.power(x, -b) + c

    popt, pcov = curve_fit(func, xdata=xdata, ydata=ydata)

    plt.plot(xdata, ydata, linestyle='--', marker='s', color=c, label=data_label)
    plt.plot(xdata, func(xdata, *popt),linestyle='-', color = 'r',
             label='fit: %5.3f * $T^{- %5.3f }$ + %5.3f'%tuple(popt))
    print('{}: '.format(data_label))
    print('The fitted temperature dependence is : $\kappa_L$ = %5.3f * $T^{- %5.3f }$ + %5.3f'%tuple(popt))
    # plt.plot(xdata, func(xdata, *popt), 'r-', label = 'fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))


def fit_file(file1, file2, file3):
    f1 = np.loadtxt(file1)
    f2 = np.loadtxt(file2)
    f3 = np.loadtxt(file3)
    fit_curve(f1[:,0], f1[:,1], c='gray', data_label='CsCaBr$_3$')
    fit_curve(f2[:, 0], f2[:, 1], c='darkviolet', data_label='CsCdBr$_3$')
    fit_curve(f3[:, 0], f3[:, 1], c='pink', data_label='CsSnBr$_3$')
    plt.xlabel('Temperature (K)')
    plt.ylabel('$\kappa_L$')
    plt.legend()
    # plt.savefig('./pictures/kappa_fit.png', dpi=400)
    plt.savefig('./pictures/kappa_fit_ha.png', dpi=400)
    plt.show()


if __name__ == '__main__':
    file1 = './Br3Ca1Cs1/phonons/pc_rta.kl'
    file2 = './Br3Ca1Cs1/scfph/kl.temp'

    file3 = './Br3Cd1Cs1/phonons/pc_rta.kl'
    file4 = './Br3Cd1Cs1/scfph2/kl.temp'

    file5 = './BrCsSn/phonons/pc_rta.kl'
    file6 = './BrCsSn/scfph/kl.temp'

    # plot_kappa(file1,file2, file3, file4, file5, file6)
    # fit_file(file2, file4, file6)
    fit_file(file1, file3, file5)