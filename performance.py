#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Python script for performance

Versions:
- 0.2 Now supports confidence intervals
- 0.3 Sextractor performance's added
- 0.4 Proper motion performance's added
- 0.5 Proper motion now creates an stats file
- 0.6 Plotting methods out to plots.py file

Todo:
    * Check useful methods

"""
from os import makedirs, path

from astropy.io import fits
from astropy.table import Table
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from numpy import arange, array, median
from pandas import concat, read_csv, DataFrame

from misc import extract_settings_luca, check_source
from misc import speeds_range
from regions_creation_luca_ssos import Create_regions


__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2017"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.6"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


class StatsPerformance:
    def __init__(self):
        """

        """
        self.prfs_d = extract_settings_luca()

    def error(self, logger, mag, sex_cf):
        """

        :param logger:
        :param mag:
        :param sex_cf:
        :return:
        """
        self.mag = mag
        errors_a_stars = []
        errors_b_stars = []
        errors_a_galaxies = []
        errors_b_galaxies = []
        errors_a_ssos = []
        errors_b_ssos = []
        keys = ['ALPHA_J2000', 'DELTA_J2000']

        # Input sources
        input_ssos_d = {}
        for d in range(1, 5, 1):
            cat_loc = '{}/{}/Catalogs'.format(self.prfs_d['fits_dir'], mag)
            cat_name = '{}/Cat_{}_d{}'.format(cat_loc, self.mag, d)
            input_ssos_d[d] = '{}.dat'.format(cat_name)
        input_ssos_d = Create_regions(input_ssos_d).check_luca(self.mag, False, True)

        # Creates a DataFrame from an input dictionary
        input_ssos_l = []
        for key_ in input_ssos_d.keys():
            input_ssos_l.append(input_ssos_d[key_])

        i_ssos_df = concat(input_ssos_l, axis=0)
        i_ssos_df = i_ssos_df.reset_index(drop=True)

        # Input stars
        input_stars_d = {}
        for d in range(1, 5, 1):
            cat_loc = '{}/{}/Catalogs'.format(self.prfs_d['fits_dir'], mag)
            cat_name = '{}/Cat_20-21_d{}'.format(cat_loc, d)
            input_stars_d[d] = '{}.dat'.format(cat_name)
        input_stars_d = Create_regions(input_stars_d).check_luca(mag, False, True)

        # Creates a DataFrame from an input dictionary
        input_stars_l = []
        for key_ in input_stars_d.keys():
            input_stars_l.append(input_stars_d[key_])

        i_stars_df = concat(input_stars_l, axis=0)
        i_stars_df = i_stars_df.reset_index(drop=True)

        # Input stars
        input_galaxies_d = {}
        for d in range(1, 5, 1):
            cat_loc = '{}/{}/Catalogs'.format(self.prfs_d['fits_dir'], mag)
            cat_name = '{}/Cat_20-21_d{}'.format(cat_loc, d)
            input_galaxies_d[d] = '{}.dat'.format(cat_name)
        input_galaxies_d = Create_regions(input_galaxies_d).check_luca(mag, False,
                                                                       True)

        # Creates a DataFrame from an input dictionary
        input_galaxies_l = []
        for key_ in input_galaxies_d.keys():
            input_galaxies_l.append(input_galaxies_d[key_])

        i_galaxies_df = concat(input_galaxies_l, axis=0)
        i_galaxies_df = i_galaxies_df.reset_index(drop=True)

        # Creates an input dictionary with all input sources
        logger.debug('checking performance for {}'.format(sex_cf))

        opts = [[0, 0], [0, 1], [0, 2], [1, 0], [1, 1],
                [1, 2], [2, 0], [2, 1], [2, 2]]
        for opt in opts:
            for dither in range(1, 5, 1):
                # Lists definition
                errors_a_stars = []
                errors_b_stars = []
                errors_a_galaxies = []
                errors_b_galaxies = []
                errors_a_ssos = []
                errors_b_ssos = []

                cat_n = 'mag_{}_CCD_x{}_y{}_d{}.cat'.format(mag, opt[0],
                                                            opt[1], dither)
                cat_o_n = '{}/{}/CCDs/{}/{}'.format(self.prfs_d['fits_dir'],
                                                    mag, sex_cf, cat_n)

                hdu_list = fits.open(cat_o_n)
                o_cat = Table(hdu_list[2].data).to_pandas()

                i_stars_d_df = i_stars_df[i_stars_df['dither_values'].isin([dither])]
                i_galaxies_d_df = i_galaxies_df[i_galaxies_df['dither_values'].isin([dither])]
                i_ssos_d_df = i_ssos_df[i_ssos_df['dither_values'].isin([dither])]
                for idx, row in o_cat.iterrows():
                    i_alpha = row['ALPHA_J2000']
                    i_delta = row['DELTA_J2000']
                    o_df = check_source(i_stars_d_df, i_alpha, i_delta, keys)
                    if o_df.empty is not True:
                        errors_a_stars.append(row['ERRA_IMAGE'])
                        errors_b_stars.append(row['ERRB_IMAGE'])
                    o_df = check_source(i_galaxies_d_df, i_alpha, i_delta, keys)
                    if o_df.empty is not True:
                        errors_a_galaxies.append(row['ERRA_IMAGE'])
                        errors_b_galaxies.append(row['ERRB_IMAGE'])
                    o_df = check_source(i_ssos_d_df, i_alpha, i_delta, keys)
                    if o_df.empty is not True:
                        errors_a_ssos.append(row['ERRA_IMAGE'])
                        errors_b_ssos.append(row['ERRB_IMAGE'])

        mean_a_stars = array(errors_a_stars).mean()
        mean_b_stars = array(errors_b_stars).mean()

        mean_a_galaxies = array(errors_a_galaxies).mean()
        mean_b_galaxies = array(errors_b_galaxies).mean()

        mean_a_ssos = array(errors_a_ssos).mean()
        mean_b_ssos = array(errors_b_ssos).mean()

        std_a_stars = array(errors_a_stars).std()
        std_b_stars = array(errors_b_stars).std()

        std_a_galaxies = array(errors_a_galaxies).std()
        std_b_galaxies = array(errors_b_galaxies).std()

        std_a_ssos = array(errors_a_ssos).std()
        std_b_ssos = array(errors_b_ssos).std()

        stats_d = {'conf': sex_cf,
                   'm_a_stars': mean_a_stars, 'm_b_stars': mean_b_stars,
                   'm_a_gals': mean_a_galaxies, 'm_b_gals': mean_b_galaxies,
                   'm_a_ssos': mean_a_ssos, 'm_b_ssos': mean_b_ssos,
                   's_a_stars': std_a_stars, 's_b_stars': std_b_stars,
                   's_a_gals': std_a_galaxies, 's_b_gals': std_b_galaxies,
                   's_a_ssos': std_a_ssos, 's_b_ssos': std_b_ssos}

        return stats_d


class PMPerformance:

    def __init__(self):
        """

        """
        self.bypassfilter = True
        self.plot = False
        self.prfs_d = extract_settings_luca()

        pass

    def check(self, logger, mag, scmp_cf, sex_cf, confidence_):
        """

        :param logger:
        :param mag:
        :param scmp_cf:
        :param sex_cf:
        :param confidence_:
        :return:
        """
        input_pm = []
        output_pm = []

        # Creates an input dictionary with all input sources
        logger.debug('checking performance for {} and {}'.format(scmp_cf,
                                                                 sex_cf))
        input_d = {}
        for d in range(1, 5, 1):
            input_d[d] = '{}/Cat_20-21_d{}.dat'.format(self.prfs_d['input_ref'],
                                                       d)
        input_d = Create_regions(input_d).check_luca(mag, True, True)

        # Creates a DataFrame from an input dictionary
        input_l = []
        for key_ in input_d.keys():
            input_l.append(input_d[key_])

        i_df = concat(input_l, axis=0)
        # Look for < 3 coincidences
        i_df = concat(g for _, g in i_df.groupby('source')
                      if len(g) >= 3)
        i_df = i_df.reset_index(drop=True)

        i_df.to_csv('input.csv')

        # Open particular file!
        filt_n = 'filt_{}_{}_4.csv'.format(scmp_cf, mag)
        filter_o_n = '{}/{}/{}/{}'.format(self.prfs_d['filter_dir'],
                                          sex_cf, scmp_cf, filt_n)

        # Cross with filtered data - Opens datafile
        o_cat = read_csv('{}'.format(filter_o_n), index_col=0)

        # Gets unique sources from input data
        unique_sources = list(set(i_df['source'].tolist()))

        # Loops over input data
        for idx_source, source_ in enumerate(unique_sources):
            # Gets associated data in input catalog
            cat = i_df[i_df['source'].isin([source_])]
            # cat has, at least, three or more rows, one for each dither

            # Creates lists for each source
            boolean_l = []
            tmp_catalog = []
            tmp_source = []
            # tmp lists
            tmp_input_pm = []
            tmp_output_pm = []

            # Iterate over each detection of each source
            for i, row in enumerate(cat.itertuples(), 1):
                catalog_n = row.catalog
                i_pm = row.pm_values
                i_alpha = row.alpha_j2000
                i_delta = row.delta_j2000

                # Checks if there is a source closed to input one
                # o_cat contains data from output (filtered) catalog
                o_df = check_source(catalog_n, o_cat, i_alpha, i_delta)

                # If there is one saves data from input data
                if o_df.empty is not True and o_df['PM'].size == 1:
                    pm_mask = self.pm_filter(o_df, i_pm, confidence_)
                    if pm_mask:
                        if o_df['SOURCE_NUMBER'].size != 1:
                            boolean_l.append('False')
                        else:
                            boolean_l.append('True')
                        tmp_catalog.append(catalog_n)
                        tmp_source.append(o_df['SOURCE_NUMBER'].iloc[0])
                        tmp_input_pm.append(i_pm)
                        o_pm = o_df['PM'].iloc[0]
                        tmp_output_pm.append(o_pm)
                else:
                    boolean_l.append('False')

            # check
            if len(tmp_input_pm) >= 3:
                input_pm.append(i_pm)
                output_pm.append(o_pm)

        pm_d = {}
        for key_ in list(set(input_pm)):
            pm_d[key_] = []

        for idx, o_pm in enumerate(output_pm):
            i_pm = input_pm[idx]
            pm_d[i_pm].append(o_pm)

        stats_d = {'pm': [], 'mean': [], 'median': [], 'std': [],
                   'max': [], 'min': [], 'detected': []}

        # mean and std determination
        for key_ in pm_d.keys():
            stats_d['pm'].append(key_)
            stats_d['mean'].append(array(pm_d[key_]).mean())
            stats_d['median'].append(median(array(pm_d[key_])))
            stats_d['std'].append(array(pm_d[key_]).std())
            stats_d['max'].append(max(pm_d[key_]))
            stats_d['min'].append(min(pm_d[key_]))
            stats_d['detected'].append(len(pm_d[key_]))
            # print(key_, len(pm_d[key_]))
            # print(array(pm_d[key_]).mean())
            # print(median(array(pm_d[key_])))
            # print(array(pm_d[key_]).std())
            # print(max(pm_d[key_]))
            # print(min(pm_d[key_]))
            # print(" ")

        if self.plot:
            if not self.plot_comp(sex_cf, scmp_cf, input_pm, output_pm,
                                  confidence_):
                raise Exception

        return stats_d

    def pm_filter(self, o_df, pm, confidence_):
        """

        :param o_df:
        :param pm:
        :param prfs_d:
        :param confidence_:
        :param bypassfilter:
        :return:
        """
        pm_ranges = speeds_range(self.prfs_d, confidence_)
        pm_range = pm_ranges[pm]

        if pm_range[0] < float(o_df['PM']) < pm_range[1]:
            return True
        elif self.bypassfilter:
            return True
        else:
            return False

    def plot_comp(self, sex_cf, scmp_cf, input_pm, output_pm, confidence_):
        """ fixme support for multiple mags

        :param sex_cf:
        :param scmp_cf:
        :param input_pm:
        :param output_pm:
        :param confidence_:
        :return:
        """
        sextractor_folder = '{}/pm/{}'.format(self.prfs_d['plots_dir'], sex_cf)
        if not path.exists(sextractor_folder):
            makedirs(sextractor_folder)

        with PdfPages('{}/pm/{}/{}_{}_cmp.pdf'.format(self.prfs_d['plots_dir'],
                                                      sex_cf, scmp_cf,
                                                      confidence_)) as pdf:
            fig = plt.figure(figsize=(16.53, 11.69), dpi=100)
            ax = fig.add_subplot(1, 1, 1)

            input_pm_cutted = filter(lambda a: a < 4, input_pm)
            output_pm_cutted = filter(lambda a: a < 4, output_pm)
            ax.plot(input_pm_cutted, output_pm_cutted, 'bs', markersize=1)
            ax.plot([0, 4], [0, 4])

            ax.set_xlabel('pm input [0 - 4]')
            ax.set_ylabel('pm output [0 - 4]')

            # x-scale
            x_major_ticks = arange(0, 4, 0.5)
            x_minor_ticks = arange(0, 4, 0.1)
            ax.set_xticks(x_major_ticks, minor=False)
            ax.set_xticks(x_minor_ticks, minor=True)

            y_major_ticks = arange(0, 4, 0.5)
            y_minor_ticks = arange(0, 4, 0.1)
            ax.set_yticks(y_major_ticks, minor=False)
            ax.set_yticks(y_minor_ticks, minor=True)

            ax.grid(b=True, which='major', linestyle='-', linewidth=2)
            ax.grid(b=True, which='minor', linestyle='--', linewidth=1)
            pdf.savefig()  # saves current figure
            plt.clf()  # clear current figure

            # [8-10]
            fig = plt.figure(figsize=(16.53, 11.69), dpi=100)
            ax = fig.add_subplot(1, 1, 1)

            input_pm_cutted = filter(lambda a: a > 4, input_pm)
            output_pm_cutted = filter(lambda a: a > 4, output_pm)
            ax.plot(input_pm_cutted, output_pm_cutted, 'bs', markersize=1)
            ax.plot([9.5, 10.5], [9.5, 10.5])

            ax.set_xlabel('pm input [9.5 - 10.5]')
            ax.set_ylabel('pm output [9.5 - 10.5]')

            # x-scale
            x_major_ticks = arange(9.5, 10.5, 0.5)
            x_minor_ticks = arange(9.5, 10.5, 0.05)
            ax.set_xticks(x_major_ticks, minor=False)
            ax.set_xticks(x_minor_ticks, minor=True)

            y_major_ticks = arange(9.5, 10.5, 0.5)
            y_minor_ticks = arange(9.5, 10.5, 0.05)
            ax.set_yticks(y_major_ticks, minor=False)
            ax.set_yticks(y_minor_ticks, minor=True)

            ax.grid(b=True, which='major', ls='-', lw=2)
            ax.grid(b=True, which='minor', ls='--', lw=1)
            pdf.savefig()  # saves current figure
            plt.clf()  # clear current figure

        return True
