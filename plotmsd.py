import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

color = ['k', 'b', 'g', 'r', 'm', 'c', 'y', 'r',
         'darkred', 'darkblue', 'darkgreen', 'darkmagenta']
lsty = ['-', '-', '-', '-', '--', '--', '--', '--', '-', '-', '-', '-']


def get_msd(file_path):
    msd_data = []
    with open(file_path, 'r') as f:
        _ = f.readline().strip().split()
        notation = f.readline().strip().split()
        for x in f.readlines():
            line = x.strip().split()
            if line[0] != '#':
                data = [float(i) for i in line]
                msd_data.append(data)
    return np.array(msd_data)


def plot_msd(ha_msd, ah_msd, fig_name):
    # 在msd文件中，1-9列为Br原子的数据。10-12列的数据为Ca，Cd的数据，13-15列为Cs的数据。
    # plt.figure(figsize=(6,6))
    plt.plot(ha_msd[:, 0], ha_msd[:, 2], linestyle='-',
             color='orange', label="HA-$U_{11}$($U_{22}$)Br")
    plt.plot(ah_msd[:, 0], ah_msd[:, 2], linestyle='--', color='orange', marker='o', markerfacecolor='w',
             label="SCP-$U_{11}$($U_{22}$)Br")

    plt.plot(ha_msd[:, 0], ha_msd[:, 1], linestyle='-',
             color='gray', label="HA-$U_{33}$Br")
    plt.plot(ah_msd[:, 0], ah_msd[:, 1], linestyle='--', color='gray', marker='o', markerfacecolor='w',
             label="SCP-$U_{33}$Br")

    if fig_name.split('.')[0]=='Br3CsSn_msd':
        #### for the BrCsSn system, the atom positions are Br Cs Sn, so it's a little different from other two materials.
        plt.plot(ha_msd[:, 0], ha_msd[:, 11], linestyle='-',
                 color='cornflowerblue', label="HA-$U_{ii}$Cs")
        plt.plot(ah_msd[:, 0], ah_msd[:, 11], linestyle='--', color='cornflowerblue', marker='s', markerfacecolor='w',
                 label="SCP-$U_{ii}$Cs")

        plt.plot(ha_msd[:, 0], ha_msd[:, 13], linestyle='-',
                 color='gold', label="HA-$U_{ii}$Sn")
        plt.plot(ah_msd[:, 0], ah_msd[:, 13], linestyle='--', color='gold', marker='s', markerfacecolor='w',
                 label="SCP-$U_{ii}$Sn")

        plt.plot(300, 0.147, marker='s', color='orange')
        plt.plot(300, 0.025, marker='s', color='gray')
    else:
        plt.plot(ha_msd[:, 0], ha_msd[:, 13], linestyle='-',
                 color='cornflowerblue', label="HA-$U_{ii}$Cs")
        plt.plot(ah_msd[:, 0], ah_msd[:, 13], linestyle='--', color='cornflowerblue', marker='s', markerfacecolor='w',
                 label="SCP-$U_{ii}$Cs")

    plt.xlabel("Temperature (K)", fontsize=14)
    plt.ylabel("Mean Square Displacement ($\AA^{2}$)", fontsize=14, labelpad=20)
    # xmin=0
    # plt.axis([xmin, xmax, ymin, ymax])
    plt.xlim(left=0, right=1000)
    plt.minorticks_on()
    plt.tick_params(axis='both', which='both', direction='in', top=True, right=True)

    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(loc='upper left', prop={'size': 12})

    plt.tight_layout()
    plt.gca().set_aspect(1.0 / plt.gca().get_data_ratio(), adjustable='box')
    # plt.savefig(fig_name, dpi=600)
    plt.show()

# ############------Br3Ca1Cs1-------#####################
# file_path_1 = './Br3Ca1Cs1/msd/pc_dos.msd'
# file_path_2 = './Br3Ca1Cs1/msd/pc_scfph.scph_msd'
# ha_msd = get_msd(file_path_1)
# ah_msd = get_msd(file_path_2)
# fig_name = 'Br3Ca1Cs1_msd.png'
# plot_msd(ha_msd,ah_msd, fig_name)

# #############------Br3Cs1Sn-------#####################
file_path_1 = './BrCsSn/msd/pc_dos.msd'
file_path_2 = './BrCsSn/msd/pc_scfph.scph_msd'
ha_msd = get_msd(file_path_1)
ah_msd = get_msd(file_path_2)
fig_name = 'Br3CsSn_msd.eps'
plot_msd(ha_msd,ah_msd, fig_name)

############------Br3CdCs-------#####################
# file_path_1 = './Br3Cd1Cs1/msd/pc_dos.msd'
# file_path_2 = './Br3Cd1Cs1/msd/pc_scfph.scph_msd'
# ha_msd = get_msd(file_path_1)
# ah_msd = get_msd(file_path_2)
# fig_name = 'Br3Cd1Cs1_msd.png'
# plot_msd(ha_msd,ah_msd, fig_name)
