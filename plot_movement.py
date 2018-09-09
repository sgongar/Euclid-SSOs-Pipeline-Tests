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

from misc import extract_settings_elvis, check_source, setting_logger
from misc import get_norm_mag, get_norm_speed
from misc_cats import get_dither


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
        self.prfs_d = extract_settings_elvis()
        self.mags = ['20-21', '21-22', '22-23', '23-24', '24-25', '25-26']

        self.create_catalogue_dict()
        self.plot_movement_to_pdf()

    def create_catalogue_dict(self):
        """

        :return: catalogue
        """
        for pm_ in self.prfs_d['pms']:
            dir_ = 'false_positives/false_{}.csv'.format(pm_)
            self.cats_d[pm_] = read_csv(dir_, index_col=0)

    def plot_movement_to_pdf(self):
        """
        """
        for pm_ in self.prfs_d['pms']:
            for mag_ in self.mags:
                pdf_name = '{}_{}.pdf'.format(pm_, mag_)
                # Creates pdf
                with PdfPages(pdf_name) as pdf:
                    cat = self.cats_d[pm_]
                    cat = cat[cat['MAG_AUTO'].isin([mag_])]

                    unique_sources = list(set(cat['SOURCE'].tolist()))
                    for source_ in unique_sources:
                        fig = plt.figure(figsize=(16.53, 11.69), dpi=100)
                        ax = fig.add_subplot(1, 1, 1)

                        source_df = cat[cat['SOURCE'].isin([source_])]
                        ax.plot(source_df['ALPHA_J2000'].tolist(),
                                source_df['DELTA_J2000'].tolist(), 'bs')

                        pdf.savefig()

    #     for source_ in list(set(df['SOURCE'].tolist())):
    #         source_df = df[df['SOURCE'].isin([source_])]
    #         source_df.to_csv('{}.csv'.format(source_))


if __name__ == "__main__":
    PlotFalseMovement()
