# !/usr/bin/python
# -*- coding: utf-8 -*-

""" Gets
   - factors from magnitude bins

Versions:
- 0.1

Todo:
    * Improve log messages
    * Get out check_source

*GNU Terry Pratchett*
"""

from misc import extract_settings_elvis, check_source, setting_logger
from misc import get_norm_speed
from misc_cats import gets_data

from numpy import nan
from pandas import concat, DataFrame, read_csv, Series


__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def compute_factors(factors_d, stats_df):
    """
    N_meas: number of all detected sources(including false detections)
    N_se: number of simulated sources recovered by source extraction
    N_true: number of simulated input sources
    f_dr: detection rate f_pur: purity
    f_com: completeness

    f_dr = N_meas / N_true = (N_se + N_false) / N_true
    f_pur = N_se / N_meas = N_se / (N_se + N_false)
    f_com = f_dr * f_pur = N_se / N_true

    :param factors_d:
    :param stats_df:
    :return:
    """
    pm_list = [0.1, 0.3, 1.0, 3.0, 10.0]
    idxs_pm = [[4, 5, 6], [7, 8, 9], [10, 11, 12],
               [13, 14, 15], [16, 17, 18]]
    for idx_mag, data in enumerate(stats_df.itertuples(), 1):
        # Loop over the different magnitudes
        if 6 < idx_mag < 13:  # Only gets values from magnitudes [20, 26]
            mag_ = data[0]  # Gets magnitude value
            for idx_pm, pm_values in enumerate(idxs_pm): # each sublist is a different pm
                pm_ = pm_list[idx_pm]
                n_se = data[pm_values[0]]
                factors_d[mag_][pm_]['n_se'] = n_se
                n_false = data[pm_values[1]]
                factors_d[mag_][pm_]['n_false'] = n_false
                n_meas = n_se + n_false
                factors_d[mag_][pm_]['n_meas'] = n_meas
                n_true = data[pm_values[2]]
                factors_d[mag_][pm_]['n_true'] = n_true

                try:
                    f_dr = float(n_meas) / float(n_true)
                    f_dr = float("{0:.2f}".format(f_dr))
                    factors_d[mag_][pm_]['f_dr'] = f_dr
                except ZeroDivisionError:
                    factors_d[mag_][pm_]['f_dr'] = nan
                try:
                    f_pur = float(n_se) / float(n_meas)
                    f_pur = float("{0:.2f}".format(f_pur))
                    factors_d[mag_][pm_]['f_pur'] = f_pur
                except ZeroDivisionError:
                    factors_d[mag_][pm_]['f_pur'] = nan
                try:
                    f_com = float(n_se) / float(n_true)
                    f_com = float("{0:.2f}".format(f_com))
                    factors_d[mag_][pm_]['f_com'] = f_com
                except ZeroDivisionError:
                    factors_d[mag_][pm_]['f_com'] = nan
        else:
            pass

    return factors_d


def save_factors(factors_d):
    """

    :param factors_d:
    :return:
    """
    tmp_d = {'mag': [], 'pm': [], 'n_se': [],
             'n_false': [], 'n_meas': [], 'n_true': [],
             'f_pur': [], 'f_dr': [], 'f_com': []}
    pm_list = [0.1, 0.3, 1.0, 3.0, 10.0]
    mags = ['20-21', '21-22', '22-23', '23-24',
            '24-25', '25-26', '26-27']
    for mag_ in mags:
        mag_df = factors_d[mag_]

        for pm_ in pm_list:
            pm_df = mag_df[pm_]
            tmp_d['mag'].append(mag_)
            tmp_d['pm'].append(pm_)
            tmp_d['n_se'].append(pm_df['n_se'])
            tmp_d['n_false'].append(pm_df['n_false'])
            tmp_d['n_meas'].append(pm_df['n_meas'])
            tmp_d['n_true'].append(pm_df['n_true'])
            tmp_d['f_pur'].append(pm_df['f_pur'])
            tmp_d['f_dr'].append(pm_df['f_dr'])
            tmp_d['f_com'].append(pm_df['f_com'])

    tmp_d['n_scmp'] = [24, 23, 20, 17, 11, 25, 30, 21, 18, 15, 19, 37,
                       20, 20, 14, 17, 26, 24, 22, 12, 22, 19, 20, 16,
                       13, 27, 23, 16, 21, 2, 0, 0, 0, 0, 0]

    # For test reasons
    # for key_ in tmp_d.keys():
    #     print(key_, len(tmp_d[key_]))

    tmp_df = DataFrame(tmp_d, columns=['mag', 'pm', 'n_se',
                                       'n_false', 'n_meas', 'n_scmp', 'n_true',
                                       'f_pur', 'f_dr', 'f_com'])
    tmp_df.to_csv('stats.csv')

    return True


def get_dither(catalog_n):
    """ Gets the number of the dither for a catalogue number.

    :param catalog_n: The catalogue with the unknown dither.
    :return: dither_n, the number of the dither.
    """
    dither_n = 0  # In order to have a clear code
    cats_dict = {}

    for i in range(1, 5, 1):
        cats_dict[i] = []

    for cat_ccd in range(0, 144, 4):
        for cat_dither in range(1, 5, 1):
            cats_dict[cat_dither].append(cat_dither + cat_ccd)

    for dither_ in cats_dict.keys():
        if catalog_n in cats_dict[dither_]:
            dither_n = dither_

    return dither_n



def get_norm_mag(o_mag):
    """

    :param o_mag:
    :return: mag_bin
    """
    mags = [[11, 12], [12, 13], [13, 14], [14, 15], [15, 16], [16, 17],
            [17, 18], [18, 19], [19, 20], [20, 21], [21, 22], [22, 23],
            [23, 24], [24, 25], [25, 26], [26, 27], [27, 28], [28, 29]]

    mag_bin = ''
    for mag_ in mags:
        if mag_[0] < o_mag < mag_[1]:
            mag_bin = '{}-{}'.format(mag_[0], mag_[1])

    return mag_bin


def splits_by_mag_bin():
    """

    :return: total_d
    """
    prfs_d = extract_settings_elvis()
    total_d = {}
    not_rejected_by_scamp = True
    all_of_them = False
    if not_rejected_by_scamp:
        total_d = {'14-15': {0: 0, 0.01: 0, 0.03: 0, 0.1: 0, 0.3: 0,
                             1.0: 0, 3.0: 0, 10.0: 0, 30.0: 0},
                   '15-16': {0: 0, 0.01: 0, 0.03: 0, 0.1: 0, 0.3: 0,
                             1.0: 0, 3.0: 0, 10.0: 0, 30.0: 0},
                   '16-17': {0: 0, 0.01: 0, 0.03: 0, 0.1: 0, 0.3: 0,
                             1.0: 0, 3.0: 0, 10.0: 0, 30.0: 0},
                   '17-18': {0: 0, 0.01: 0, 0.03: 0, 0.1: 0, 0.3: 0,
                             1.0: 0, 3.0: 0, 10.0: 0, 30.0: 0},
                   '18-19': {0: 0, 0.01: 0, 0.03: 0, 0.1: 0, 0.3: 0,
                             1.0: 0, 3.0: 0, 10.0: 0, 30.0: 0},
                   '19-20': {0: 0, 0.01: 0, 0.03: 0, 0.1: 0, 0.3: 0,
                             1.0: 0, 3.0: 0, 10.0: 0, 30.0: 0},
                   '20-21': {0: 0, 0.01: 5, 0.03: 21, 0.1: 24, 0.3: 23,
                             1.0: 20, 3.0: 17, 10.0: 11, 30.0: 17},
                   '21-22': {0: 0, 0.01: 3, 0.03: 22, 0.1: 25, 0.3: 30,
                             1.0: 21, 3.0: 18, 10.0: 15, 30.0: 15},
                   '22-23': {0: 0, 0.01: 3, 0.03: 26, 0.1: 19, 0.3: 37,
                             1.0: 20, 3.0: 20, 10.0: 14, 30.0: 10},
                   '23-24': {0: 0, 0.01: 7, 0.03: 29, 0.1: 17, 0.3: 26,
                             1.0: 24, 3.0: 22, 10.0: 12, 30.0: 15},
                   '24-25': {0: 0, 0.01: 9, 0.03: 23, 0.1: 22, 0.3: 19,
                             1.0: 20, 3.0: 16, 10.0: 13, 30.0: 19},
                   '25-26': {0: 0, 0.01: 10, 0.03: 21, 0.1: 27, 0.3: 23,
                             1.0: 16, 3.0: 21, 10.0: 2, 30.0: 16},
                   '26-27': {0: 0, 0.01: 0, 0.03: 0, 0.1: 0, 0.3: 0,
                             1.0: 0, 3.0: 0, 10.0: 0, 30.0: 0},
                   '27-28': {0: 0, 0.01: 0, 0.03: 0, 0.1: 0, 0.3: 0,
                             1.0: 0, 3.0: 0, 10.0: 0, 30.0: 0}}

    elif all_of_them:
        # opens input_catalog
        cat_ssos = read_csv('cats/cat_clean_ssos.csv', index_col=0)

        mags = [[20, 21], [21, 22], [22, 23], [23, 24],
                [24, 25], [25, 26], [26, 27], [27, 28]]
        pms = prfs_d['pms']
        pms_range = {0.01: [0.005, 0.015], 0.03: [0.015, 0.05],
                     0.1: [0.05, 0.15], 0.3: [0.15, 0.5],
                     1.0: [0.5, 1.5], 3.0: [1.5, 5],
                     10.0: [5.0, 15.0], 30.0: [15.0, 50]}

        total_d['14-15'] = {}
        total_d['14-15'][0] = 0
        for pm_ in pms:
            total_d['14-15'][pm_] = 0

        total_d['15-16'] = {}
        total_d['15-16'][0] = 0
        for pm_ in pms:
            total_d['15-16'][pm_] = 0

        total_d['16-17'] = {}
        total_d['16-17'][0] = 0
        for pm_ in pms:
            total_d['16-17'][pm_] = 0

        total_d['17-18'] = {}
        total_d['17-18'][0] = 0
        for pm_ in pms:
            total_d['17-18'][pm_] = 0

        total_d['18-19'] = {}
        total_d['18-19'][0] = 0
        for pm_ in pms:
            total_d['18-19'][pm_] = 0

        total_d['19-20'] = {}
        total_d['19-20'][0] = 0
        for pm_ in pms:
            total_d['19-20'][pm_] = 0

        for mag_ in mags:
            mag_bin_cat = cat_ssos[cat_ssos['ABMAG'] < mag_[1]]
            mag_bin_cat = mag_bin_cat[mag_bin_cat['ABMAG'] > mag_[0]]
            total_d['{}-{}'.format(mag_[0], mag_[1])] = {}
            total_d['{}-{}'.format(mag_[0], mag_[1])][0] = 0
            for pm_ in pms:
                # gets proper motion bin
                pm_bin_cat = mag_bin_cat[mag_bin_cat['VEL'] < pms_range[pm_][1]]
                pm_bin_cat = pm_bin_cat[pm_bin_cat['VEL'] > pms_range[pm_][0]]
                try:
                    pm_bin_cat = concat(g for _, g in pm_bin_cat.groupby('SOURCE')
                                        if len(g) >= 3)
                    total_sources = len(list(set(pm_bin_cat['SOURCE'].tolist())))
                    total_d['{}-{}'.format(mag_[0], mag_[1])][pm_] = total_sources
                except ValueError:
                    total_d['{}-{}'.format(mag_[0], mag_[1])][pm_] = 0

        print(total_d)

    return total_d


def redo_data_d():
    """ Creates a dictionary
    TODO - Automatic number!

    :return: tmp_d
    """
    prfs_d = extract_settings_elvis()
    total_d = splits_by_mag_bin()

    data_d = {}
    mags = [[14, 15], [15, 16], [16, 17], [17, 18], [18, 19], [19, 20],
            [20, 21], [21, 22], [22, 23], [23, 24], [24, 25], [25, 26],
            [26, 27], [27, 28]]
    pms = prfs_d['pms']
    for mag_ in mags:
        mag_bin = '{}-{}'.format(mag_[0], mag_[1])
        data_d[mag_bin] = {}
        data_d[mag_bin][0] = {'right': 0, 'false': 0,
                              'total': total_d[mag_bin][0]}
        for pm_ in pms:
            data_d[mag_bin][pm_] = {'right': 0, 'false': 0,
                                    'total': total_d[mag_bin][pm_]}

    return data_d


def redo_factors_d():
    """ Creates a dictionary
    TODO - Automatic number!

    :return: factors_d
    """
    prfs_d = extract_settings_elvis()

    factors_d = {}
    mags = [[20, 21], [21, 22], [22, 23], [23, 24],
            [24, 25], [25, 26], [26, 27]]
    pms = prfs_d['pms']
    for mag_ in mags:
        mag_bin = '{}-{}'.format(mag_[0], mag_[1])
        factors_d[mag_bin] = {}
        for pm_ in pms:
            factors_d[mag_bin][pm_] = {'f_pur': 0.0, 'f_dr': 0.0, 'f_com': 0.0,
                                       'n_se': 0, 'n_false': 0, 'n_meas': 0,
                                       'n_true': 0}

    return factors_d


class FactorsScampPerformance:

    def __init__(self):
        """ Me da los valores de salida de todos las estrellas  y galaxias
        presentes en filt 3 obtenidos o no

        """
        self.filter_p_number = 3  # First one with enough data for statistics
        self.prfs_d = extract_settings_elvis()
        self.data_d = redo_data_d()
        factors_d = redo_factors_d()

        self.save = True

        logger_name = 'scamp_performance'  # Set as desired
        self.logger = setting_logger(self.prfs_d, logger_name)

        filt_cat = self.gets_filtered_catalog()  # Gets data from filtered
        input_df = gets_data()  # Gets data from catalogs
        stats_df = self.extract_stats(filt_cat, input_df)  # Splits due type
        factors_d = compute_factors(factors_d, stats_df)
        save_factors(factors_d)

    def gets_filtered_catalog(self):
        """

        :return: filtered_cat, filtered catalog
        """
        filter_n = 'filt__{}.csv'.format(self.filter_p_number)
        filter_o_n = '{}/{}'.format(self.prfs_d['filtered'], filter_n)

        print('Opens filtered catalogue {}'.format(filter_o_n))
        filtered_cat = read_csv('{}'.format(filter_o_n), index_col=0)

        return filtered_cat

    def extract_stats(self, filt_cat, input_df):
        """

        :param filt_cat:
        :param input_df:
        :return:
        """
        # Unique sources (?)
        prfs_d = extract_settings_elvis()
        unique_sources = list(set(filt_cat['SOURCE_NUMBER'].tolist()))

        # False positives dictionary
        false_positives = {1: {'RA': [], 'DEC': []}, 2: {'RA': [], 'DEC': []},
                           3: {'RA': [], 'DEC': []}, 4: {'RA': [], 'DEC': []}}

        yey = 0
        nop = 0

        print('Creating new catalogues from filtered catalogue due type')
        print('Unique sources: {}'.format(len(unique_sources)))
        for idx_source_, source_ in enumerate(unique_sources):
            source_df = filt_cat[filt_cat['SOURCE_NUMBER'].isin([source_])]
            # Takes the first value of MAG Series
            o_mag_bin = get_norm_mag(source_df['MEDIAN_MAG_ISO'].iloc[0])
            # Takes the first value of PM Series
            o_pm_norm = get_norm_speed(source_df['PM'].iloc[0])

            source_d = {'source': [], 'pm': [], 'mag': []}
            right_detections = 0
            for i, row in enumerate(source_df.itertuples(), 1):
                alpha = row.ALPHA_J2000
                delta = row.DELTA_J2000
                dither_n = get_dither(int(row.CATALOG_NUMBER))
                # Checks object type
                keys = ['RA', 'DEC']  # Catalogue version 2
                test_sso = check_source(input_df[dither_n]['SSOs'],
                                        alpha, delta, keys)
                if test_sso.empty is not True:
                    right_detections += 1
                    source_d['source'].append(row.SOURCE_NUMBER)
                    pm_norm = get_norm_speed(float(test_sso['VEL'].iloc[0]))
                    source_d['pm'].append(pm_norm)  # order by input not output
                    source_d['mag'].append(test_sso['ABMAG'].iloc[0])
                else:
                    false_positives[dither_n]['RA'].append(alpha)
                    false_positives[dither_n]['DEC'].append(delta)

            if right_detections >= 4:
                i_mag_bin = get_norm_mag(source_d['mag'][0])
                i_pm_norm = get_norm_speed(source_d['pm'][0])
                if i_mag_bin == '24-25' and i_pm_norm == 1.0:
                    yey += 1
                    print('yey {}'.format(yey))
                self.data_d[i_mag_bin][i_pm_norm]['right'] += 1
            else:
                if o_mag_bin == '24-25' and o_pm_norm == 1.0:
                    nop += 1
                    print('nop {}'.format(nop))
                self.data_d[o_mag_bin][o_pm_norm]['false'] += 1

        for dither_ in false_positives.keys():
            # print(false_positives[dither_]['RA'])
            alpha_list = false_positives[dither_]['RA']
            alpha_serie = Series(alpha_list, name='ALPHA_J2000')
            delta_list = false_positives[dither_]['DEC']
            delta_serie = Series(delta_list, name='DELTA_J2000')

            output = concat([alpha_serie, delta_serie], axis=1)
            cat_name = 'false_positives/regions/false_{}.reg'.format(dither_)
            output.to_csv(cat_name, index=False, header=False, sep=" ")

        mags = [[14, 15], [15, 16], [16, 17], [17, 18], [18, 19],
                [19, 20], [20, 21], [21, 22], [22, 23], [23, 24],
                [24, 25], [25, 26], [26, 27], [27, 28]]

        idxs = ['14-15', '15-16', '16-17', '17-18', '18-19',
                '19-20', '20-21', '21-22', '22-23', '23-24',
                '24-25', '25-26', '26-27', '27-28']
        pms = prfs_d['pms']
        stats_d = {'idx': idxs}
        for pm_ in pms:
            stats_d['right-{}'.format(pm_)] = []
            stats_d['false-{}'.format(pm_)] = []
            stats_d['total-{}'.format(pm_)] = []
            for mag_ in mags:
                mag_bin = '{}-{}'.format(mag_[0], mag_[1])
                stats_d['right-{}'.format(pm_)].append(self.data_d[mag_bin][pm_]['right'])
                stats_d['false-{}'.format(pm_)].append(self.data_d[mag_bin][pm_]['false'])
                stats_d['total-{}'.format(pm_)].append(self.data_d[mag_bin][pm_]['total'])

        stats_df = DataFrame(stats_d,
                             columns=['idx', 'right-0', 'false-0', 'total-0',
                                      'right-0.1', 'false-0.1', 'total-0.1',
                                      'right-0.3', 'false-0.3', 'total-0.3',
                                      'right-1.0', 'false-1.0', 'total-1.0',
                                      'right-3.0', 'false-3.0', 'total-3.0',
                                      'right-10.0', 'false-10.0', 'total-10.0'])
        stats_df = stats_df.set_index('idx')

        return stats_df


if __name__ == "__main__":
    FactorsScampPerformance()
