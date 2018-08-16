# !/usr/bin/python
# -*- coding: utf-8 -*-

"""Python script for performance

Versions:
- 0.1

Todo:
    * Improve log messages - create logger instance
    * Multi-threading support
    * Get number of input sources for each situation - automatically

*GNU Terry Pratchett*
"""
from sys import argv

from matplotlib import cm
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from numpy import array, isnan, nan, nanmean
from pandas import concat, DataFrame, read_csv, Series

from misc import extract_settings_luca
from regions import Create_regions

__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def compute_factors(stats_d, tmp_d):
    """
    N_meas: number of all detected sources(including false detections)
    N_se: number of simulated sources recovered by source extraction
    N_true: number of simulated input sources
    f_dr: detection rate f_pur: purity
    f_com: completeness

    f_dr = N_meas / N_true = (N_se + N_false) / N_true
    f_pur = N_se / N_meas = N_se / (N_se + N_false)
    f_com = f_dr * f_pur = N_se / N_true

    :param stats_d:
    :param tmp_d:
    :return:
    """
    prfs_d = extract_settings_luca()

    for mag_ in prfs_d['mags']:
        for pm_ in prfs_d['pms']:
            n_meas = tmp_d[mag_][pm_]['right'] + tmp_d[mag_][pm_]['false']
            # stats_d[mag_][pm_]['n_meas'].append(n_meas)
            stats_d[mag_][pm_]['n_meas'] = n_meas
            n_false = tmp_d[mag_][pm_]['false']
            # stats_d[mag_][pm_]['n_false'].append(n_false)
            stats_d[mag_][pm_]['n_false'] = n_false
            n_se = tmp_d[mag_][pm_]['right']
            # stats_d[mag_][pm_]['n_se'].append(n_se)
            stats_d[mag_][pm_]['n_se'] = n_se
            n_true = tmp_d[mag_][pm_]['total']
            # stats_d[mag_][pm_]['n_true'].append(n_true)
            stats_d[mag_][pm_]['n_true'] = n_true
            # Factors computation
            try:
                f_dr = float(n_meas) / float(n_true)
                f_dr = float("{0:.2f}".format(f_dr))
                stats_d[mag_][pm_]['f_dr'] = f_dr
            except ZeroDivisionError:
                stats_d[mag_][pm_]['f_dr'] = nan
            try:
                f_pur = float(n_se) / float(n_meas)
                f_pur = float("{0:.2f}".format(f_pur))
                stats_d[mag_][pm_]['f_pur'] = f_pur
            except ZeroDivisionError:
                stats_d[mag_][pm_]['f_pur'] = nan
            try:
                f_com = float(n_se) / float(n_true)
                f_com = float("{0:.2f}".format(f_com))
                stats_d[mag_][pm_]['f_com'] = f_com
            except ZeroDivisionError:
                stats_d[mag_][pm_]['f_com'] = nan

    return stats_d


def redo_check_d():
    """ Creates a dictionary

    :return:
    """
    check_d = {'detections': 0}

    return check_d


def redo_stats_d():
    """ Creates a dictionary

    :return: tmp_d
    """
    prfs_d = extract_settings_luca()
    stats_d = {}
    for mag_ in prfs_d['mags']:
        stats_d[mag_] = {}
        for pm_ in prfs_d['pms']:
            stats_d[mag_][pm_] = {'n_meas': [], 'n_false': [], 'n_se': [],
                                  'n_true': [], 'f_dr': [], 'f_pur': [],
                                  'f_com': []}

    return stats_d


def redo_tmp_d():
    """ Creates a dictionary
    TODO - Automatic number!

    :return: tmp_d
    """
    prfs_d = extract_settings_luca()
    total_d = {'20-21': {0: 0, 0.001: 21, 0.003: 22, 0.01: 14, 0.03: 17,
                         0.1: 19, 0.3: 22, 1: 23, 3: 19, 10: 18, 30: 21},
               '21-22': {0: 0, 0.001: 22, 0.003: 21, 0.01: 19, 0.03: 24,
                         0.1: 16, 0.3: 23, 1: 18, 3: 24, 10: 24, 30: 19},
               '22-23': {0: 0, 0.001: 18, 0.003: 24, 0.01: 21, 0.03: 22,
                         0.1: 16, 0.3: 22, 1: 20, 3: 20, 10: 21, 30: 26},
               '23-24': {0: 0, 0.001: 21, 0.003: 19, 0.01: 25, 0.03: 23,
                         0.1: 22, 0.3: 19, 1: 20, 3: 17, 10: 20, 30: 21},
               '24-25': {0: 0, 0.001: 22, 0.003: 19, 0.01: 20, 0.03: 14,
                         0.1: 22, 0.3: 17, 1: 21, 3: 20, 10: 21, 30: 18},
               '25-26': {0: 0, 0.001: 14, 0.003: 23, 0.01: 19, 0.03: 27,
                         0.1: 25, 0.3: 29, 1: 23, 3: 21, 10: 16, 30: 22},
               '26-27': {0: 0, 0.001: 21, 0.003: 20, 0.01: 22, 0.03: 22,
                         0.1: 22, 0.3: 20, 1: 18, 3: 25, 10: 23, 30: 14}}

    tmp_d = {}
    for mag_ in prfs_d['mags']:
        tmp_d[mag_] = {}
        tmp_d[mag_][0] = {'right': 0, 'false': 0, 'total': total_d[mag_][0]}
        for pm_ in prfs_d['pms']:
            tmp_d[mag_][pm_] = {'right': 0, 'false': 0,
                                'total': total_d[mag_][pm_]}

    return tmp_d


def check_source(catalog_n, o_cat, i_alpha, i_delta):
    """

    :param catalog_n:
    :param o_cat:
    :param i_alpha:
    :param i_delta:
    :return:
    """
    prfs_d = extract_settings_luca()

    o_df = o_cat[o_cat['catalog'].isin([catalog_n])]
    o_df = o_df[o_df['alpha_j2000'] + prfs_d['tolerance'] > i_alpha]
    o_df = o_df[i_alpha > o_df['alpha_j2000'] - prfs_d['tolerance']]
    o_df = o_df[o_df['delta_j2000'] + prfs_d['tolerance'] > i_delta]
    o_df = o_df[i_delta > o_df['delta_j2000'] - prfs_d['tolerance']]

    return o_df


def speeds_range(prfs_d, confidence):
    """ given a confidence value returns a dict with speeds

    @param prfs_d:
    @param confidence:

    @return speeds_dict:
    """
    speeds_dict = {}
    for pm_ in prfs_d['pms']:
        speeds_dict[pm_] = [pm_ - pm_ * confidence / 100,
                            pm_ + pm_ * confidence / 100]

    return speeds_dict


def get_norm_speed(o_pm):
    """

    :return:
    """
    prfs_d = extract_settings_luca()

    speeds_d = speeds_range(prfs_d, 50)

    pm_norm = 0
    for key_ in speeds_d.keys():
        low = speeds_d[key_][0]
        high = speeds_d[key_][1]
        if low < o_pm < high:
            pm_norm = key_

    return pm_norm


def creates_input_df(input_dict, mag, save):
    """ Creates an input dataframe from an input dictionary.

    :return: input dataframe
    """
    occurrences = 3

    input_list = []
    for key_ in input_dict.keys():
        input_list.append(input_dict[key_])

    input_df = concat(input_list, axis=0)
    # Look for >= 3 coincidences
    input_df = concat(g for _, g in input_df.groupby('source')
                      if len(g) >= occurrences)
    input_df = input_df.reset_index(drop=True)

    if save:
        print('Saves input catalog')
        input_df.to_csv('tmp/inputs_{}.csv'.format(mag))

    return input_df


def gets_filtered_catalog(scmp_cf, sex_cf, mag, filter_p_number):
    """

    :return: filtered_cat, filtered catalog
    """
    prfs_d = extract_settings_luca()
    filter_n = 'filt_{}_{}_{}.csv'.format(scmp_cf, mag, filter_p_number)
    filter_o_n = '{}/{}/{}/{}/{}'.format(prfs_d['filtered'], mag, sex_cf,
                                         scmp_cf, filter_n)

    print('Opens filtered catalog {}'.format(filter_n))
    # Cross with filtered data - Opens datafile
    filtered_cat = read_csv('{}'.format(filter_o_n), index_col=0)

    return filtered_cat


class ScampOutputFactorsSSOs:

    def __init__(self):
        """

        """
        self.filter_p_number = 9
        self.prfs_d = extract_settings_luca()

        self.scmp_cf = '10_1.1_0.5_0.04'
        self.sex_cf = '30_1.5_1.5_0.01_4'
        self.fig_size = [16.53, 11.69]
        self.dpi = 100

        self.save = True
        self.tmp_d = redo_tmp_d()

        i_mag_l = []
        i_pm_l = []
        n_se_l = []
        n_false_l = []
        n_true_l = []
        n_meas_l = []
        f_pur_l = []
        f_dr_l = []
        f_com_l = []

        for mag_ in self.prfs_d['mags']:
            self.mag = mag_
            stats_d = self.check_pm_distribution()

            for pm_ in self.prfs_d['pms']:
                i_mag_l.append(self.mag)
                i_pm_l.append(pm_)
                n_se_l.append(stats_d[mag_][pm_]['n_se'])
                n_false_l.append(stats_d[mag_][pm_]['n_false'])
                n_true_l.append(stats_d[mag_][pm_]['n_true'])
                n_meas_l.append(stats_d[mag_][pm_]['n_meas'])
                f_pur_l.append(stats_d[mag_][pm_]['f_pur'])
                f_dr_l.append(stats_d[mag_][pm_]['f_pur'])
                f_com_l.append(stats_d[mag_][pm_]['f_com'])

        i_mag_s = Series(i_mag_l, name='mag_bin')
        i_pm_s = Series(i_pm_l, name='i_pm')
        n_se_s = Series(n_se_l, name='n_se')
        n_false_s = Series(n_false_l, name='n_false')
        n_true_s = Series(n_true_l, name='n_true')
        n_meas_s = Series(n_meas_l, name='n_meas')
        f_dr_s = Series(f_dr_l, name='f_dr')
        f_pur_s = Series(f_pur_l, name='f_pur')
        f_com_s = Series(f_com_l, name='f_com')

        self.out_df = concat([i_mag_s, i_pm_s, n_se_s, n_false_s, n_true_s,
                              n_meas_s, f_dr_s, f_pur_s, f_com_s], axis=1)
        self.out_df.to_csv('total.csv')

        self.plot()

    def creates_input_dict(self):
        """ Creates an input dictionary. Each key contains SSOs' information
        for each dither.

        :return: input_dict
        """
        input_dict = {}
        # Loops over the four dithers
        for dither in range(1, 5, 1):
            cat_location = '{}/{}/Catalogs'.format(self.prfs_d['fits_dir'],
                                                   self.mag)
            cat_name = '{}/Cat_{}_d{}'.format(cat_location, self.mag, dither)
            input_dict[dither] = '{}.dat'.format(cat_name)
        input_ssos = Create_regions(input_dict).check_ssos(self.mag, True)

        return input_ssos

    def check_pm_distribution(self):
        """

        :return:
        """
        # Creates an input dictionary with all input sources
        print('Magnitude bin: {}'.format(self.mag))
        print('Scamp configuration: {}'.format(self.scmp_cf))
        print('Sextractor configuration: {}'.format(self.sex_cf))

        input_ssos, input_stars, input_galaxies = self.creates_input_dict()
        input_ssos_df = creates_input_df(input_ssos, self.save, self.mag)

        # Open particular file!
        filter_cat = gets_filtered_catalog(self.scmp_cf, self.sex_cf, self.mag,
                                           self.filter_p_number)
        # Filter by mag!
        filter_cat = filter_cat[filter_cat['MAG_ISO'] < float(self.mag[3:5])]
        filter_cat = filter_cat[filter_cat['MAG_ISO'] > float(self.mag[0:2])]

        # Gets unique sources from input data
        # unique_sources = list(set(input_ssos_df['source'].tolist()))
        unique_sources = list(set(filter_cat['SOURCE_NUMBER'].tolist()))
        unique_sources = filter(lambda a: a > 0, unique_sources)

        sources_n = len(unique_sources)
        self.tmp_d[self.mag]['total'] = sources_n

        print('Input sources to be analysed {}'.format(sources_n))
        # Loops over input data (Luca's catalog)
        for idx_source, source_ in enumerate(unique_sources):
            # Gets associated data in input catalog
            source_df = filter_cat[filter_cat['SOURCE_NUMBER'].isin([source_])]

            if int(source_df['SOURCE_NUMBER'].size) >= 3:
                pm_right = False
                pm_false = False
                check_d = redo_check_d()  # Creates a dictionary
                # Iterate over each detection of each source
                for i, row in enumerate(source_df.itertuples(), 1):
                    catalog_n = row.CATALOG_NUMBER
                    o_alpha = row.ALPHA_J2000
                    o_delta = row.DELTA_J2000
                    o_pm = row.PM

                    # Checks if there is a source closed to input one
                    o_df = check_source(catalog_n, input_ssos_df,
                                        o_alpha, o_delta)

                    if o_df.empty is not True and o_df['pm_values'].size == 1:
                        check_d['detections'] += 1
                        pm_right = float(o_df['pm_values'])
                    else:
                        pm_false = get_norm_speed(o_pm)
                        pass

                if check_d['detections'] >= 3:
                    self.tmp_d[self.mag][pm_right]['right'] += 1
                else:
                    self.tmp_d[self.mag][pm_false]['false'] += 1

        stats_d = redo_stats_d()
        stats_d = compute_factors(stats_d, self.tmp_d)

        return stats_d

    def plot(self):
        """

        :return:
        """
        tmp_mag = []
        tmp_pm = []
        tmp_f_com = []
        tmp_f_pur = []

        with PdfPages('completeness.pdf') as pdf:
            fig = plt.figure(figsize=self.fig_size, dpi=self.dpi)
            ax = fig.add_subplot(1, 1, 1)

            for mag_ in self.prfs_d['mags']:
                for pm_ in self.prfs_d['pms']:
                    # print(mag_, pm_)
                    out_df = self.out_df[self.out_df['mag_bin'].isin([mag_])]
                    out_df = out_df[out_df['i_pm'].isin([pm_])]
                    f_com = float(out_df['f_com'])

                    # print('f_com {}'.format(f_com))
                    p_mag = 0
                    if mag_ == '20-21':
                        p_mag = 20.5
                    elif mag_ == '21-22':
                        p_mag = 21.5
                    elif mag_ == '22-23':
                        p_mag = 22.5
                    elif mag_ == '23-24':
                        p_mag = 23.5
                    elif mag_ == '24-25':
                        p_mag = 24.5
                    elif mag_ == '25-26':
                        p_mag = 25.5
                    elif mag_ == '26-27':
                        p_mag = 26.5

                    tmp_mag.append(p_mag)
                    tmp_pm.append(pm_)
                    tmp_f_com.append(f_com)

                    plt.scatter(p_mag, pm_, c=f_com, cmap=cm.copper,
                                vmin=0., vmax=1.)

            ax.set_yscale('log')

            plt.colorbar()
            plt.grid(True)
            pdf.savefig()

        with PdfPages('purity.pdf') as pdf:
            fig = plt.figure(figsize=(16.53, 11.69), dpi=100)
            ax = fig.add_subplot(1, 1, 1)

            for mag_ in self.prfs_d['mags']:
                for pm_ in self.prfs_d['pms']:
                    # print(mag_, pm_)
                    out_df = self.out_df[self.out_df['mag_bin'].isin([mag_])]
                    out_df = out_df[out_df['i_pm'].isin([pm_])]
                    f_pur = float(out_df['f_pur'])

                    p_mag = 0
                    if mag_ == '20-21':
                        p_mag = 20.5
                    elif mag_ == '21-22':
                        p_mag = 21.5
                    elif mag_ == '22-23':
                        p_mag = 22.5
                    elif mag_ == '23-24':
                        p_mag = 23.5
                    elif mag_ == '24-25':
                        p_mag = 24.5
                    elif mag_ == '25-26':
                        p_mag = 25.5
                    elif mag_ == '26-27':
                        p_mag = 26.5

                    tmp_f_pur.append(f_pur)

                    plt.scatter(p_mag, pm_, c=f_pur, cmap=cm.copper,
                                vmin=0., vmax=1.)

            ax.set_yscale('log')

            plt.colorbar()
            plt.grid(True)
            pdf.savefig()


class ScampOutputValuesSSOs:

    def __init__(self):
        """

        """
        self.filter_p_number = 3  # We are going extract all values
        self.prfs_d = extract_settings_luca()

        self.scmp_cf = '10_1.1_0.5_0.04'
        self.sex_cf = '30_1.5_1.5_0.01_4'

        self.save = True
        self.tmp_d = redo_tmp_d()
        self.input_sources_d = {}

        for mag_ in self.prfs_d['mags']:
            self.mag = mag_
            self.input_sources_d[self.mag] = {'SSO': [], 'star': [],
                                              'galaxy': []}
            self.check_filtered_sources()
            self.create_stats()

    def creates_input_dict(self):
        """ Creates an input dictionary. Each key contains SSOs' information
        for each dither.

        :return: input_dict
        """
        input_dict = {}
        # Loops over the four dithers
        for dither in range(1, 5, 1):
            cat_location = '{}/{}/Catalogs'.format(self.prfs_d['fits_dir'],
                                                   self.mag)
            cat_name = '{}/Cat_{}_d{}'.format(cat_location, self.mag, dither)
            input_dict[dither] = '{}.dat'.format(cat_name)
        input_ssos = Create_regions(input_dict).check_ssos(self.mag, True)

        input_dict = {}
        # Loops over the four dithers
        for dither in range(1, 5, 1):
            cat_location = '{}/{}/Catalogs'.format(self.prfs_d['fits_dir'],
                                                   self.mag)
            cat_name = '{}/Cat_{}_d{}'.format(cat_location, self.mag, dither)
            input_dict[dither] = '{}.dat'.format(cat_name)
        input_stars = Create_regions(input_dict).check_stars(self.mag, True)

        input_dict = {}
        # Loops over the four dithers
        for dither in range(1, 5, 1):
            cat_location = '{}/{}/Catalogs'.format(self.prfs_d['fits_dir'],
                                                   self.mag)
            cat_name = '{}/Cat_{}_d{}'.format(cat_location, self.mag, dither)
            input_dict[dither] = '{}.dat'.format(cat_name)
        input_galaxies = Create_regions(input_dict).check_galaxies(self.mag,
                                                                   True)
        return input_ssos, input_stars, input_galaxies

    def check_filtered_sources(self):
        """

        :return:
        """
        # Creates an input dictionary with all input sources
        print('Magnitude bin: {}'.format(self.mag))
        print('Scamp configuration: {}'.format(self.scmp_cf))
        print('Sextractor configuration: {}'.format(self.sex_cf))

        input_ssos, input_stars, input_galaxies = self.creates_input_dict()
        input_ssos_df = creates_input_df(input_ssos, self.save, self.mag)
        input_stars_df = creates_input_df(input_stars, self.save, self.mag)
        input_galaxies_df = creates_input_df(input_galaxies,
                                             self.save, self.mag)

        # Open particular file!
        filter_cat = gets_filtered_catalog(self.scmp_cf, self.sex_cf, self.mag,
                                           self.filter_p_number)
        # Filter by mag!
        filter_cat = filter_cat[filter_cat['MAG_ISO'] < float(self.mag[3:5])]
        filter_cat = filter_cat[filter_cat['MAG_ISO'] > float(self.mag[0:2])]

        # Gets unique sources from input data
        # unique_sources = list(set(input_ssos_df['source'].tolist()))
        unique_sources = list(set(filter_cat['SOURCE_NUMBER'].tolist()))
        unique_sources = filter(lambda a: a > 0, unique_sources)

        sources_n = len(unique_sources)
        self.tmp_d[self.mag]['total'] = sources_n

        print('Input sources to be analysed {}'.format(sources_n))
        # Loops over input data (Luca's catalog)
        for idx_source, source_ in enumerate(unique_sources):
            # Gets associated data in input catalog
            source_df = filter_cat[filter_cat['SOURCE_NUMBER'].isin([source_])]

            if int(source_df['SOURCE_NUMBER'].size) >= 1:  # 1 - for everything
                # pm_right = False
                # pm_false = False
                # check_d = redo_check_d()  # Creates a dictionary
                # Iterate over each detection of each source
                for i, row in enumerate(source_df.itertuples(), 1):
                    catalog_n = row.CATALOG_NUMBER
                    o_alpha = row.ALPHA_J2000
                    o_delta = row.DELTA_J2000

                    # Checks if there is a source closed to input one
                    sso_df = check_source(catalog_n, input_ssos_df,
                                          o_alpha, o_delta)

                    star_df = check_source(catalog_n, input_stars_df,
                                           o_alpha, o_delta)

                    galaxy_df = check_source(catalog_n, input_galaxies_df,
                                             o_alpha, o_delta)

                    if sso_df.empty is not True:
                        self.input_sources_d[self.mag]['SSO'].append(source_)
                    elif star_df.empty is not True:
                        self.input_sources_d[self.mag]['star'].append(source_)
                    elif galaxy_df.empty is not True:
                        self.input_sources_d[self.mag]['galaxy'].append(source_)
                    else:
                        pass

    def create_stats(self):
        """

        todo - improve code!
        :return:
        """
        pms = [0, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0]

        # SSO stuff
        sso_d = {'catalog_n': [], 'median_a_image': [],
                 'median_erra_image': [], 'median_b_image': [],
                 'median_errb_image': [], 'median_class_star': [],
                 'ellipticity': [], 'median_mag_iso': [],
                 'median_magerr_iso': [], 'output_pm': [],
                 'median_flux_iso': []}
        for idx in range(0, len(self.prfs_d['pms']) + 1, 1):
            sso_d['median_mag_iso'].append([])
            sso_d['median_magerr_iso'].append([])
            sso_d['catalog_n'].append([])
            sso_d['median_a_image'].append([])
            sso_d['median_erra_image'].append([])
            sso_d['median_b_image'].append([])
            sso_d['median_errb_image'].append([])
            sso_d['median_class_star'].append([])
            sso_d['ellipticity'].append([])
            sso_d['output_pm'].append([])
            sso_d['median_flux_iso'].append([])

        o_df = gets_filtered_catalog(self.scmp_cf, self.sex_cf, self.mag,
                                     self.filter_p_number)
        ssos = self.input_sources_d[self.mag]['SSO']
        ssos_df = o_df[o_df['SOURCE_NUMBER'].isin(ssos)]

        # Gets unique sources from input data
        unique_sources = list(set(ssos_df['SOURCE_NUMBER'].tolist()))
        unique_sources = array(unique_sources)
        unique_sources = unique_sources[~isnan(unique_sources)]

        unique_sources = [value for value in unique_sources if value != nan]
        sources_n = len(unique_sources)

        print('Fixed SSOs to be analysed {}'.format(sources_n))
        # Loops over input data (Luca's catalog)
        for idx_source, source_ in enumerate(unique_sources):
            print('Star stats {} - Total {}'.format(idx_source, sources_n))
            # Gets associated data in input catalog
            cat_df = ssos_df[ssos_df['SOURCE_NUMBER'].isin([source_])]

            o_pm = cat_df['PM'].iloc[0]
            o_pm_norm = get_norm_speed(o_pm)

            idx = pms.index(o_pm_norm)

            median_mag_iso = cat_df['MEDIAN_MAG_ISO'].iloc[0]
            sso_d['median_mag_iso'][idx].append(median_mag_iso)
            median_magerr_iso = cat_df['MEDIAN_MAGERR_ISO'].iloc[0]
            sso_d['median_magerr_iso'][idx].append(median_magerr_iso)
            median_a_image_ = cat_df['MEDIAN_A_IMAGE'].iloc[0]
            sso_d['median_a_image'][idx].append(median_a_image_)
            median_erra_image_ = cat_df['MEDIAN_ERRA_IMAGE'].iloc[0]
            sso_d['median_erra_image'][idx].append(median_erra_image_)
            median_b_image_ = cat_df['MEDIAN_B_IMAGE'].iloc[0]
            sso_d['median_b_image'][idx].append(median_b_image_)
            median_errb_image_ = cat_df['MEDIAN_ERRB_IMAGE'].iloc[0]
            sso_d['median_errb_image'][idx].append(median_errb_image_)
            median_class_star_ = cat_df['MEDIAN_CLASS_STAR'].iloc[0]
            sso_d['median_class_star'][idx].append(median_class_star_)
            ellipticity_ = cat_df['ELLIPTICITY'].iloc[0]
            sso_d['ellipticity'][idx].append(ellipticity_)
            output_pm_ = o_pm  # really needed?
            sso_d['output_pm'][idx].append(output_pm_)
            median_flux_iso_ = cat_df['MEDIAN_FLUX_ISO'].iloc[0]
            sso_d['median_flux_iso'][idx].append(median_flux_iso_)

        keys = ['median_mag_iso', 'median_magerr_iso', 'median_a_image',
                'median_erra_image', 'median_b_image', 'median_errb_image',
                'median_flux_iso', 'median_class_star', 'ellipticity',
                'output_pm']

        for sso_key in keys:
            df2_keys = [0, 0.001, 0.003, 0.01, 0.03,
                        0.1, 0.3, 1, 3, 10, 30]
            sso_d2 = {0: [], 0.001: [], 0.003: [], 0.01: [], 0.03: [],
                      0.1: [], 0.3: [], 1: [], 3: [], 10: [], 30: []}

            sso_df = DataFrame(sso_d[sso_key])
            for column_ in sso_df.columns:
                for idx, value_ in enumerate(sso_df[column_]):
                    sso_d2[df2_keys[idx]].append(value_)

            for key_ in sso_d2.keys():
                sso_d2[key_].append('m {}'.format(nanmean(sso_d2[key_])))

            sso_df2 = DataFrame(sso_d2)
            sso_df2_filename = 'f_{}_{}_{}.csv'.format(self.mag, sso_key,
                                                       self.filter_p_number)
            sso_df2_dir = 'output/stats_ssos'
            sso_df2.to_csv('{}/{}'.format(sso_df2_dir, sso_df2_filename))

        # Stars stuff
        star_d = {'catalog_n': [], 'median_a_image': [],
                  'median_erra_image': [], 'median_b_image': [],
                  'median_errb_image': [], 'median_class_star': [],
                  'ellipticity': [], 'median_mag_iso': [],
                  'median_magerr_iso': [], 'output_pm': [],
                  'median_flux_iso': []}
        for idx in range(0, len(self.prfs_d['pms']) + 1, 1):
            star_d['median_mag_iso'].append([])
            star_d['median_magerr_iso'].append([])
            star_d['catalog_n'].append([])
            star_d['median_a_image'].append([])
            star_d['median_erra_image'].append([])
            star_d['median_b_image'].append([])
            star_d['median_errb_image'].append([])
            star_d['median_class_star'].append([])
            star_d['ellipticity'].append([])
            star_d['output_pm'].append([])
            star_d['median_flux_iso'].append([])

        o_df = gets_filtered_catalog(self.scmp_cf, self.sex_cf, self.mag,
                                     self.filter_p_number)
        stars = self.input_sources_d[self.mag]['star']
        stars_df = o_df[o_df['SOURCE_NUMBER'].isin(stars)]

        # Gets unique sources from input data
        unique_sources = list(set(stars_df['SOURCE_NUMBER'].tolist()))
        unique_sources = array(unique_sources)
        unique_sources = unique_sources[~isnan(unique_sources)]

        unique_sources = [value for value in unique_sources if value != nan]
        sources_n = len(unique_sources)

        print('Fixed stars to be analysed {}'.format(sources_n))
        # Loops over input data (Luca's catalog)
        for idx_source, source_ in enumerate(unique_sources):
            print('Star stats {} - Total {}'.format(idx_source, sources_n))
            # Gets associated data in input catalog
            cat_df = stars_df[stars_df['SOURCE_NUMBER'].isin([source_])]

            o_pm = cat_df['PM'].iloc[0]
            o_pm_norm = get_norm_speed(o_pm)

            idx = pms.index(o_pm_norm)

            median_mag_iso = cat_df['MEDIAN_MAG_ISO'].iloc[0]
            star_d['median_mag_iso'][idx].append(median_mag_iso)
            median_magerr_iso = cat_df['MEDIAN_MAGERR_ISO'].iloc[0]
            star_d['median_magerr_iso'][idx].append(median_magerr_iso)
            median_a_image_ = cat_df['MEDIAN_A_IMAGE'].iloc[0]
            star_d['median_a_image'][idx].append(median_a_image_)
            median_erra_image_ = cat_df['MEDIAN_ERRA_IMAGE'].iloc[0]
            star_d['median_erra_image'][idx].append(median_erra_image_)
            median_b_image_ = cat_df['MEDIAN_B_IMAGE'].iloc[0]
            star_d['median_b_image'][idx].append(median_b_image_)
            median_errb_image_ = cat_df['MEDIAN_ERRB_IMAGE'].iloc[0]
            star_d['median_errb_image'][idx].append(median_errb_image_)
            median_class_star_ = cat_df['MEDIAN_CLASS_STAR'].iloc[0]
            star_d['median_class_star'][idx].append(median_class_star_)
            ellipticity_ = cat_df['ELLIPTICITY'].iloc[0]
            star_d['ellipticity'][idx].append(ellipticity_)
            output_pm_ = o_pm  # really needed?
            star_d['output_pm'][idx].append(output_pm_)
            median_flux_iso_ = cat_df['MEDIAN_FLUX_ISO'].iloc[0]
            star_d['median_flux_iso'][idx].append(median_flux_iso_)

        keys = ['median_mag_iso', 'median_magerr_iso', 'median_a_image',
                'median_erra_image', 'median_b_image', 'median_errb_image',
                'median_flux_iso', 'median_class_star', 'ellipticity',
                'output_pm']

        for star_key in keys:
            df2_keys = [0, 0.001, 0.003, 0.01, 0.03,
                        0.1, 0.3, 1, 3, 10, 30]
            star_d2 = {0: [], 0.001: [], 0.003: [], 0.01: [], 0.03: [],
                       0.1: [], 0.3: [], 1: [], 3: [], 10: [], 30: []}

            star_df = DataFrame(star_d[star_key])
            for column_ in star_df.columns:
                for idx, value_ in enumerate(star_df[column_]):
                    star_d2[df2_keys[idx]].append(value_)

            for key_ in star_d2.keys():
                star_d2[key_].append('m {}'.format(nanmean(star_d2[key_])))

            star_df2 = DataFrame(star_d2)
            star_df2_filename = 'f_{}_{}_{}.csv'.format(self.mag, star_key,
                                                        self.filter_p_number)
            star_df2_dir = 'output/stats_stars'
            star_df2.to_csv('{}/{}'.format(star_df2_dir, star_df2_filename))

        # Galaxies stuff
        galaxy_d = {'catalog_n': [], 'median_a_image': [],
                    'median_erra_image': [], 'median_b_image': [],
                    'median_errb_image': [], 'median_class_star': [],
                    'ellipticity': [], 'median_mag_iso': [],
                    'median_magerr_iso': [], 'output_pm': [],
                    'median_flux_iso': []}
        for idx in range(0, len(self.prfs_d['pms']) + 1, 1):
            galaxy_d['median_mag_iso'].append([])
            galaxy_d['median_magerr_iso'].append([])
            galaxy_d['catalog_n'].append([])
            galaxy_d['median_a_image'].append([])
            galaxy_d['median_erra_image'].append([])
            galaxy_d['median_b_image'].append([])
            galaxy_d['median_errb_image'].append([])
            galaxy_d['median_class_star'].append([])
            galaxy_d['ellipticity'].append([])
            galaxy_d['output_pm'].append([])
            galaxy_d['median_flux_iso'].append([])

        o_df = gets_filtered_catalog(self.scmp_cf, self.sex_cf, self.mag,
                                     self.filter_p_number)
        galaxies = self.input_sources_d[self.mag]['galaxy']
        galaxies_df = o_df[o_df['SOURCE_NUMBER'].isin(galaxies)]

        # Gets unique sources from input data
        unique_sources = list(set(galaxies_df['SOURCE_NUMBER'].tolist()))
        unique_sources = array(unique_sources)
        unique_sources = unique_sources[~isnan(unique_sources)]

        unique_sources = [value for value in unique_sources if value != nan]
        sources_n = len(unique_sources)

        # Organiza los objetos por PM
        print('Fixed galaxies to be analysed {}'.format(sources_n))
        # Loops over input data (Luca's catalog)
        for idx_source, source_ in enumerate(unique_sources):
            print('Galaxy stats {} - Total {}'.format(idx_source, sources_n))
            # Gets associated data in input catalog
            cat_df = galaxies_df[galaxies_df['SOURCE_NUMBER'].isin([source_])]

            o_pm = cat_df['PM'].iloc[0]
            o_pm_norm = get_norm_speed(o_pm)

            idx = pms.index(o_pm_norm)

            median_mag_iso = cat_df['MEDIAN_MAG_ISO'].iloc[0]
            galaxy_d['median_mag_iso'][idx].append(median_mag_iso)
            median_magerr_iso = cat_df['MEDIAN_MAGERR_ISO'].iloc[0]
            galaxy_d['median_magerr_iso'][idx].append(median_magerr_iso)
            median_a_image_ = cat_df['MEDIAN_A_IMAGE'].iloc[0]
            galaxy_d['median_a_image'][idx].append(median_a_image_)
            median_erra_image_ = cat_df['MEDIAN_ERRA_IMAGE'].iloc[0]
            galaxy_d['median_erra_image'][idx].append(median_erra_image_)
            median_b_image_ = cat_df['MEDIAN_B_IMAGE'].iloc[0]
            galaxy_d['median_b_image'][idx].append(median_b_image_)
            median_errb_image_ = cat_df['MEDIAN_ERRB_IMAGE'].iloc[0]
            galaxy_d['median_errb_image'][idx].append(median_errb_image_)
            median_class_star_ = cat_df['MEDIAN_CLASS_STAR'].iloc[0]
            galaxy_d['median_class_star'][idx].append(median_class_star_)
            ellipticity_ = cat_df['ELLIPTICITY'].iloc[0]
            galaxy_d['ellipticity'][idx].append(ellipticity_)
            output_pm_ = o_pm  # really needed?
            galaxy_d['output_pm'][idx].append(output_pm_)
            median_flux_iso_ = cat_df['MEDIAN_FLUX_ISO'].iloc[0]
            galaxy_d['median_flux_iso'][idx].append(median_flux_iso_)

        keys = ['median_mag_iso', 'median_magerr_iso', 'median_a_image',
                'median_erra_image', 'median_b_image', 'median_errb_image',
                'median_flux_iso', 'median_class_star', 'ellipticity',
                'output_pm']

        for galaxy_key in keys:
            df2_keys = [0, 0.001, 0.003, 0.01, 0.03,
                        0.1, 0.3, 1, 3, 10, 30]
            galaxy_d2 = {0: [], 0.001: [], 0.003: [], 0.01: [], 0.03: [],
                         0.1: [], 0.3: [], 1: [], 3: [], 10: [], 30: []}

            galaxy_df = DataFrame(galaxy_d[galaxy_key])
            for column_ in galaxy_df.columns:
                for idx, value_ in enumerate(galaxy_df[column_]):
                    galaxy_d2[df2_keys[idx]].append(value_)

            for key_ in galaxy_d2.keys():
                galaxy_d2[key_].append('m {}'.format(nanmean(galaxy_d2[key_])))

            galaxy_df2 = DataFrame(galaxy_d2)
            galaxy_df2_filename = 'f_{}_{}_{}.csv'.format(self.mag, galaxy_key,
                                                          self.filter_p_number)
            galaxy_df2_dir = 'output/stats_galaxies'
            galaxy_df2.to_csv('{}/{}'.format(galaxy_df2_dir,
                                             galaxy_df2_filename))


if __name__ == "__main__":
    """ todo - improve description!
    
    """
    if argv[1] == '-factor':
        ScampOutputFactorsSSOs()
    elif argv[1] == '-value':
        ScampOutputValuesSSOs()
