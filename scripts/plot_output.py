#!/usr/bin/env python3
"""
Script to read in the ionosphere output and plot the result.
"""
import numpy as np
import matplotlib.pyplot as plt


def read_output_file(filename):
    """Read the IE output file"""
    with open(filename, 'r') as f:

        names = [x.strip() for x in f.readline().split(',')[2:]]
        print(names)
        nTheta = int(f.readline().split()[-1])
        nPhi = int(f.readline().split()[-1])
        nVars = len(names) - 2

        Buffer = dict()
        for variable in names:
            Buffer[variable] = np.zeros((nTheta, nPhi))

        for line in f:

            linesplit = line.split()
            for ivar, variable in enumerate(names):

                i = int(linesplit[0])
                j = int(linesplit[1])

                Buffer[variable][i - 1, j - 1] = float(linesplit[2 + ivar])

    return Buffer


def plot_output(Buffer):
    """Plot the output"""
    fig, ax = plt.subplots(3, 1, figsize=(8.5, 11), gridspec_kw=dict(hspace=0.4))

    for i, vname in enumerate(['Epsilon', 'Jr', 'Potential']):
        mp = ax[i].contourf(Buffer['Phi'], Buffer['Theta'], Buffer[vname])
        ax[i].set_title(vname)
        ax[i].set_xlabel('Phi')
        ax[i].set_ylabel('Theta')
        ax[i].set_aspect('equal')

        size = 0.02
        bbox = ax[i].get_position()
        cax = ax[i].figure.add_axes([bbox.x1 + 0.01, bbox.y0, 0.02, bbox.y1 - bbox.y0])
        cbar = ax[i].figure.colorbar(mp, cax=cax)
        # cbar.set_label(self.cbarlabel, rotation=270, labelpad=1.0)

    # fig.tight_layout()
    plt.show()

def plot_residual():
    """
    """
    residual_files = ['output/runlog_91x180_gmres',
                      'output/runlog_91x180_bicgstab',
                      'output/runlog_181x360_gmres',
                      'output/runlog_181x360_bicgstab']

    all_iter = list()
    all_residual = list()

    for rfile in residual_files:
        iterations = list()
        residual = list()
        with open(rfile, 'r') as f:
            for line in f:
                if 'matvecs' in line:
                    iterations.append(int(line.split()[0]))
                    residual.append(float(line.split()[-1]))

        all_iter.append(iterations)
        all_residual.append(residual)

    fig = plt.figure()
    ax = plt.subplot(111)
    for iter, res in zip(all_iter, all_residual):
        ax.semilogy(iter, res)

    ax.set_xlabel('Iterations')
    ax.set_ylabel('Residual (relative)')

    ax.legend([
        'GMRES (91x180)',
        'BICGSTAB (91x180)',
        'GMRES (181x360)',
        'BICGSTAB (181x360)'])

    plt.show()


# ========================================================================
if __name__ == '__main__':

    filename = 'output/output.dat'
    Buffer = read_output_file(filename)
    plot_output(Buffer)
    # plot_residual()
