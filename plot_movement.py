# !/usr/bin/python
# -*- coding: utf-8 -*-

"""

Versions:
- 0.1

Todo:
    * Improve log messages

*GNU Terry Pratchett*
"""

from astropy.io import fits
from astropy.table import Table
from pandas import concat, read_csv, Series
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


class PlotFalseMovement:
    def __init__(self):
        """
        """
        self.cats_d = {}
        self.pms = [0.1, 0.3]
        # self.pms = [0.1]
        self.mags = ['20-21', '21-22']
        self.mode = 'false'

        self.create_catalogue_dict()
        self.plot_movement_to_pdf()

    def create_catalogue_dict(self):
        """

        :return: catalogue
        """
        for pm_ in self.pms:
            dir_ = 'false_positives/{}_{}.csv'.format(self.mode, pm_)
            self.cats_d[pm_] = read_csv(dir_, index_col=0)

    def plot_movement_to_pdf(self):
        """
        """
        pm_ = 0.1
        for mag_ in self.mags:
            pdf_name = '{}_{}_{}.pdf'.format(self.mode, pm_, mag_)
            # Creates pdf
            with PdfPages(pdf_name) as pdf:
                cat = self.cats_d[pm_]
                cat = cat[cat['MAG_AUTO'].isin([mag_])]

                unique_sources = list(set(cat['SOURCE'].tolist()))
                for source_ in unique_sources:
                    fig = plt.figure(figsize=(16.53, 11.69), dpi=100)
                    ax = fig.add_subplot(1, 1, 1)
                    title_ = 'source: {}\npm: {} - mag: {}'.format(source_,
                                                                   pm_,
                                                                   mag_)
                    ax.set_title(title_)

                    source_df = cat[cat['SOURCE'].isin([source_])]

                    dithers = source_df['DITHER'].tolist()

                    alpha_list = source_df['ALPHA_J2000'].tolist()
                    alpha_list_seconds = []
                    for alpha_ in alpha_list:
                        alpha_list_seconds.append(alpha_ * 3600)
                    # Both in seconds
                    max_alpha = max(alpha_list_seconds)
                    min_alpha = min(alpha_list_seconds)

                    delta_list = source_df['DELTA_J2000'].tolist()
                    delta_list_seconds = []
                    for delta_ in delta_list:
                        delta_list_seconds.append(delta_ * 3600)
                    # Both in seconds
                    max_delta = max(delta_list_seconds)
                    min_delta = min(delta_list_seconds)

                    # Ticks alpha creation
                    x_ticks = {}
                    x_ticks_label = {}
                    # steps_alpha = int((max_alpha - min_alpha) / 0.01)
                    x_ticks['major_t'] = []
                    x_ticks_label['major_t'] = []
                    x_ticks['major_t'].append(min_alpha - 0.01)
                    x_first_label = min_alpha - 0.01 - int(min_alpha - 0.01)
                    x_first_label = float("{0:.3f}".format(x_first_label))
                    x_ticks_label['major_t'].append(x_first_label)
                    x_ticks['major_t'].append(min_alpha)
                    x_second_label = min_alpha - int(min_alpha)
                    x_second_label = float("{0:.3f}".format(x_second_label))
                    x_ticks_label['major_t'].append(x_second_label)

                    idx_alpha = 1
                    alpha_step = min_alpha
                    while alpha_step < max_alpha:
                        alpha_step = min_alpha + idx_alpha * 0.01
                        x_ticks['major_t'].append(alpha_step)
                        x_step_label = alpha_step - int(alpha_step)
                        x_step_label = float("{0:.3f}".format(x_step_label))
                        x_ticks_label['major_t'].append(x_step_label)
                        idx_alpha += 1

                    x_ticks['minor_t'] = []
                    x_ticks['minor_t'].append(min_alpha - 0.005)
                    x_ticks['minor_t'].append(min_alpha)

                    idx_alpha = 1
                    alpha_step = min_alpha
                    while alpha_step < max_alpha:
                        alpha_step = min_alpha + idx_alpha * 0.005
                        x_ticks['minor_t'].append(alpha_step)
                        idx_alpha += 1

                    # Ticks delta creation
                    y_ticks = {}
                    y_ticks_label = {}
                    # steps_alpha = int((max_alpha - min_alpha) / 0.01)
                    y_ticks['major_t'] = []
                    y_ticks_label['major_t'] = []
                    y_ticks['major_t'].append(min_delta - 0.01)
                    y_first_label = min_delta - 0.01 - int(min_delta - 0.01)
                    y_first_label = float("{0:.3f}".format(y_first_label))
                    y_ticks_label['major_t'].append(y_first_label)
                    y_ticks['major_t'].append(min_delta)
                    y_second_label = min_delta - int(min_delta)
                    y_second_label = float("{0:.3f}".format(y_second_label))
                    y_ticks_label['major_t'].append(y_second_label)

                    idx_delta = 1
                    delta_step = min_delta
                    while delta_step < max_delta:
                        delta_step = min_delta + idx_delta * 0.01
                        y_ticks['major_t'].append(delta_step)
                        y_step_label = delta_step - int(delta_step)
                        y_step_label = float("{0:.3f}".format(y_step_label))
                        y_ticks_label['major_t'].append(y_step_label)
                        idx_delta += 1

                    y_ticks['minor_t'] = []
                    y_ticks['minor_t'].append(min_delta - 0.005)
                    y_ticks['minor_t'].append(min_delta)
                    idx_delta = 1
                    delta_step = min_delta
                    while delta_step < max_delta:
                        delta_step = min_delta + idx_delta * 0.005
                        y_ticks['minor_t'].append(delta_step)
                        idx_delta += 1

                    colors = ['r', 'b', 'y', 'g']
                    for idx in range(0, len(alpha_list_seconds), 1):
                        ax.scatter(alpha_list_seconds[idx],
                                   delta_list_seconds[idx],
                                   c=colors[idx],
                                   label='dither_{}'.format(dithers[idx]), s=36)

                    for idx, txt in enumerate(dithers):
                        ax.annotate(txt, (alpha_list_seconds[idx],
                                          delta_list_seconds[idx]),
                                    xytext=(alpha_list_seconds[idx] + 0.0002,
                                            delta_list_seconds[idx] + 0.0002))

                    # Test
                    if len(alpha_list_seconds) == 3:

                        ax.plot(alpha_list_seconds[:2],
                                delta_list_seconds[:2])
                        ax.plot(alpha_list_seconds[1:3],
                                delta_list_seconds[1:3])
                        ax.plot([alpha_list_seconds[0],
                                 alpha_list_seconds[2]],
                                [delta_list_seconds[0],
                                 delta_list_seconds[2]], linestyle='--')
                    elif len(alpha_list_seconds) == 4:
                        ax.plot(alpha_list_seconds[:2],
                                delta_list_seconds[:2])
                        ax.plot(alpha_list_seconds[1:3],
                                delta_list_seconds[1:3])
                        ax.plot(alpha_list_seconds[2:4],
                                delta_list_seconds[2:4])
                        ax.plot([alpha_list_seconds[0],
                                 alpha_list_seconds[3]],
                                [delta_list_seconds[0],
                                 delta_list_seconds[3]], linestyle='--')

                    ax.set_xticks(x_ticks['major_t'], minor=False)
                    ax.set_xticklabels(x_ticks_label['major_t'])
                    ax.set_xticks(x_ticks['minor_t'], minor=True)

                    ax.set_yticks(y_ticks['major_t'], minor=False)
                    ax.set_yticklabels(y_ticks_label['major_t'])
                    ax.set_yticks(y_ticks['minor_t'], minor=True)

                    # Formats grids
                    ax.grid(b=True, which='major',
                            linestyle='-', linewidth=1)
                    ax.grid(b=True, which='minor',
                            linestyle='--', linewidth=1)

                    ax.legend()

                    pdf.savefig()
                    plt.close(fig)

        pm_ = 0.3
        for mag_ in self.mags:
            pdf_name = '{}_{}_{}.pdf'.format(self.mode, pm_, mag_)
            # Creates pdf
            with PdfPages(pdf_name) as pdf:
                cat = self.cats_d[pm_]
                cat = cat[cat['MAG_AUTO'].isin([mag_])]

                unique_sources = list(set(cat['SOURCE'].tolist()))
                for source_ in unique_sources:
                    fig = plt.figure(figsize=(16.53, 11.69), dpi=100)
                    ax = fig.add_subplot(1, 1, 1)
                    title_ = 'source: {}\npm: {} - mag: {}'.format(source_,
                                                                   pm_,
                                                                   mag_)
                    ax.set_title(title_)

                    source_df = cat[cat['SOURCE'].isin([source_])]

                    dithers = source_df['DITHER'].tolist()

                    alpha_list = source_df['ALPHA_J2000'].tolist()
                    alpha_list_seconds = []
                    for alpha_ in alpha_list:
                        alpha_list_seconds.append(alpha_ * 3600)
                    # Both in seconds
                    max_alpha = max(alpha_list_seconds)
                    min_alpha = min(alpha_list_seconds)

                    delta_list = source_df['DELTA_J2000'].tolist()
                    delta_list_seconds = []
                    for delta_ in delta_list:
                        delta_list_seconds.append(delta_ * 3600)
                    # Both in seconds
                    max_delta = max(delta_list_seconds)
                    min_delta = min(delta_list_seconds)

                    # Ticks alpha creation
                    x_ticks = {}
                    x_ticks_label = {}
                    # steps_alpha = int((max_alpha - min_alpha) / 0.01)
                    x_ticks['major_t'] = []
                    x_ticks_label['major_t'] = []
                    x_ticks['major_t'].append(min_alpha - 0.03)
                    x_first_label = min_alpha - 0.03 - int(min_alpha - 0.03)
                    x_first_label = float("{0:.3f}".format(x_first_label))
                    x_ticks_label['major_t'].append(x_first_label)
                    x_ticks['major_t'].append(min_alpha)
                    x_second_label = min_alpha - int(min_alpha)
                    x_second_label = float("{0:.3f}".format(x_second_label))
                    x_ticks_label['major_t'].append(x_second_label)

                    idx_alpha = 1
                    alpha_step = min_alpha
                    while alpha_step < max_alpha:
                        alpha_step = min_alpha + idx_alpha * 0.03
                        x_ticks['major_t'].append(alpha_step)
                        x_step_label = alpha_step - int(alpha_step)
                        x_step_label = float("{0:.3f}".format(x_step_label))
                        x_ticks_label['major_t'].append(x_step_label)
                        idx_alpha += 1

                    x_ticks['minor_t'] = []
                    x_ticks['minor_t'].append(min_alpha - 0.01)
                    x_ticks['minor_t'].append(min_alpha)

                    idx_alpha = 1
                    alpha_step = min_alpha
                    while alpha_step < max_alpha:
                        alpha_step = min_alpha + idx_alpha * 0.01
                        x_ticks['minor_t'].append(alpha_step)
                        idx_alpha += 1

                    # Ticks delta creation
                    y_ticks = {}
                    y_ticks_label = {}
                    # steps_alpha = int((max_alpha - min_alpha) / 0.01)
                    y_ticks['major_t'] = []
                    y_ticks_label['major_t'] = []
                    y_ticks['major_t'].append(min_delta - 0.03)
                    y_first_label = min_delta - 0.03 - int(min_delta - 0.03)
                    y_first_label = float("{0:.3f}".format(y_first_label))
                    y_ticks_label['major_t'].append(y_first_label)
                    y_ticks['major_t'].append(min_delta)
                    y_second_label = min_delta - int(min_delta)
                    y_second_label = float("{0:.3f}".format(y_second_label))
                    y_ticks_label['major_t'].append(y_second_label)

                    idx_delta = 1
                    delta_step = min_delta
                    while delta_step < max_delta:
                        delta_step = min_delta + idx_delta * 0.03
                        y_ticks['major_t'].append(delta_step)
                        y_step_label = delta_step - int(delta_step)
                        y_step_label = float("{0:.3f}".format(y_step_label))
                        y_ticks_label['major_t'].append(y_step_label)
                        idx_delta += 1

                    y_ticks['minor_t'] = []
                    y_ticks['minor_t'].append(min_delta - 0.01)
                    y_ticks['minor_t'].append(min_delta)
                    idx_delta = 1
                    delta_step = min_delta
                    while delta_step < max_delta:
                        delta_step = min_delta + idx_delta * 0.01
                        y_ticks['minor_t'].append(delta_step)
                        idx_delta += 1

                    colors = ['r', 'b', 'y', 'g']
                    for idx in range(0, len(alpha_list_seconds), 1):
                        ax.scatter(alpha_list_seconds[idx],
                                   delta_list_seconds[idx],
                                   c=colors[idx],
                                   label='dither_{}'.format(dithers[idx]), s=36)

                    for idx, txt in enumerate(dithers):
                        ax.annotate(txt, (alpha_list_seconds[idx],
                                          delta_list_seconds[idx]),
                                    xytext=(alpha_list_seconds[idx] + 0.0002,
                                            delta_list_seconds[idx] + 0.0002))

                    # Test
                    if len(alpha_list_seconds) == 3:
                        ax.plot(alpha_list_seconds[:2],
                                delta_list_seconds[:2])
                        ax.plot(alpha_list_seconds[1:3],
                                delta_list_seconds[1:3])
                        ax.plot([alpha_list_seconds[0],
                                 alpha_list_seconds[2]],
                                [delta_list_seconds[0],
                                 delta_list_seconds[2]], linestyle='--')
                    elif len(alpha_list_seconds) == 4:
                        ax.plot(alpha_list_seconds[:2],
                                delta_list_seconds[:2])
                        ax.plot(alpha_list_seconds[1:3],
                                delta_list_seconds[1:3])
                        ax.plot(alpha_list_seconds[2:4],
                                delta_list_seconds[2:4])
                        ax.plot([alpha_list_seconds[0],
                                 alpha_list_seconds[3]],
                                [delta_list_seconds[0],
                                 delta_list_seconds[3]], linestyle='--')

                    ax.set_xticks(x_ticks['major_t'], minor=False)
                    ax.set_xticklabels(x_ticks_label['major_t'])
                    ax.set_xticks(x_ticks['minor_t'], minor=True)

                    ax.set_yticks(y_ticks['major_t'], minor=False)
                    ax.set_yticklabels(y_ticks_label['major_t'])
                    ax.set_yticks(y_ticks['minor_t'], minor=True)

                    # Formats grids
                    ax.grid(b=True, which='major',
                            linestyle='-', linewidth=1)
                    ax.grid(b=True, which='minor',
                            linestyle='--', linewidth=1)

                    ax.legend()

                    pdf.savefig()
                    plt.close(fig)

    #     for source_ in list(set(df['SOURCE'].tolist())):
    #         source_df = df[df['SOURCE'].isin([source_])]
    #         source_df.to_csv('{}.csv'.format(source_))


if __name__ == "__main__":
    PlotFalseMovement()
