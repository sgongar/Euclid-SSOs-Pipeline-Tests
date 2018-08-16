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

from misc import extract_settings_elvis, check_source, setting_logger
from misc_cats import get_dither


__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def extract_inputs_d():
    """

    :return:
    """
    prfs_dict = extract_settings_elvis()
    inputs_d = {}

    cat_stars_loc = prfs_dict['references']
    cat_stars = fits.open('{}/cat_stars.fits'.format(cat_stars_loc))
    stars_data = Table(cat_stars[1].data)
    stars_df = stars_data.to_pandas()
    stars_idx = range(0, 28474, 1)  # hardcoded - todo!
    stars_df['IDX'] = stars_idx
    inputs_d['stars'] = stars_df

    cat_galaxies_loc = prfs_dict['references']
    cat_galaxies = fits.open('{}/cat_galaxies.fits'.format(cat_galaxies_loc))
    galaxies_data = Table(cat_galaxies[1].data)
    galaxies_df = galaxies_data.to_pandas()
    galaxies_idx = range(0, 143766, 1)  # hardcoded - todo!
    galaxies_df['IDX'] = galaxies_idx
    inputs_d['galaxies'] = galaxies_df

    return inputs_d


def get_norm_speed(o_pm):
    """

    :return:
    """
    speeds_d = {0.01: [0.005, 0.015], 0.03: [0.015, 0.05],
                0.1: [0.05, 0.15], 0.3: [0.15, 0.5],
                1.0: [0.5, 1.5], 3.0: [1.5, 5],
                10.0: [5.0, 15.0], 30.0: [15.0, 50]}

    if o_pm < 0.005:
        pm_norm = 0
    else:
        # pm_norm = 0
        for key_ in speeds_d.keys():
            low = speeds_d[key_][0]
            high = speeds_d[key_][1]
            if low < o_pm <= high:
                pm_norm = key_

    return pm_norm


def get_norm_mag(o_mag):
    """

    :param o_mag:
    :return: mag_bin
    """
    mags = [[14, 15], [15, 16], [16, 17], [17, 18], [18, 19], [19, 20],
            [20, 21], [21, 22], [22, 23], [23, 24], [24, 25], [25, 26],
            [26, 27], [27, 28]]

    mag_bin = ''
    for mag_ in mags:
        if mag_[0] < o_mag < mag_[1]:
            mag_bin = '{}-{}'.format(mag_[0], mag_[1])

    return mag_bin


def get_object(alpha, delta, input_d):
    """

    :param alpha:
    :param delta:
    :return:
    """
    o_df_galaxies = check_source(input_d['galaxies'], alpha, delta,
                                 keys=['ra', 'dec'])
    o_df_stars = check_source(input_d['stars'], alpha, delta,
                              keys=['RA2000(Gaia)', 'DEC2000(Gaia)'])

    if o_df_galaxies.empty is not True:
        return 'galaxies'
    elif o_df_stars.empty is not True:
        return 'star'
    else:
        return 'SSO'


class FalsePositivesScampPerformance:

    def __init__(self):
        """ Me da los valores de salida de todos las estrellas  y galaxias
        presentes en filt 3 obtenidos o no

        """
        self.filter_p_number = 9  # First one with enough data for statistics
        self.prfs_d = extract_settings_elvis()
        self.input_d = extract_inputs_d()

        # False positives dictionary
        self.false_positives = {1: {'RA': [], 'DEC': [], 'MAG': [],
                                    'PM': [], 'CLASS': [], 'OBJECT': []},
                                2: {'RA': [], 'DEC': [], 'MAG': [],
                                    'PM': [], 'CLASS': [], 'OBJECT': []},
                                3: {'RA': [], 'DEC': [], 'MAG': [],
                                    'PM': [], 'CLASS': [], 'OBJECT': []},
                                4: {'RA': [], 'DEC': [], 'MAG': [],
                                    'PM': [], 'CLASS': [], 'OBJECT': []}}

        self.save = True

        logger_name = 'scamp_performance'  # Set as desired
        self.logger = setting_logger(self.prfs_d, logger_name)

        filt_cat = self.gets_filtered_catalog()  # Gets data from filtered
        input_df = self.gets_data()  # Gets data from catalogs
        self.extract_stats(filt_cat, input_df)  # Splits due type

    def gets_filtered_catalog(self):
        """

        :return: filtered_cat, filtered catalog
        """
        filter_n = 'filt__{}.csv'.format(self.filter_p_number)
        filter_o_n = '{}/{}'.format(self.prfs_d['filtered'], filter_n)

        print('Opens filtered catalogue {}'.format(filter_o_n))
        filtered_cat = read_csv('{}'.format(filter_o_n), index_col=0)

        return filtered_cat

    def gets_data(self):
        """ Creates an input dictionary. Each key contains SSOs' information
        for each dither.

        :return: input_dict
        """
        # # For now we only have data for dither 1
        input_df = {1: {}, 2: {}, 3: {}, 4: {}}

        for key_ in input_df.keys():
            # ssos_cat = 'cat_ssos_{}.csv'.format(key_)
            ssos_cat = 'catalogues_input/cat_ssos_{}.csv'.format(key_)
            input_df[key_]['SSOs'] = read_csv(ssos_cat, index_col=0)
            stars_cat = 'tmp_stars/stars.csv'
            input_df[key_]['stars'] = read_csv(stars_cat, index_col=0)
            galaxies_cat = 'tmp_galaxies/galaxies.csv'
            input_df[key_]['galaxies'] = read_csv(galaxies_cat, index_col=0)

        return input_df

    def extract_stats(self, filt_cat, input_df):
        """

        :param filt_cat:
        :param input_df:
        :return:
        """
        # Unique sources (?)
        unique_sources = list(set(filt_cat['SOURCE_NUMBER'].tolist()))

        print('Creating new catalogues from filtered catalogue due type')
        print('Total sources: {}'.format(filt_cat['SOURCE_NUMBER'].size))
        for idx_source_, source_ in enumerate(unique_sources):
            source_df = filt_cat[filt_cat['SOURCE_NUMBER'].isin([source_])]
            # Takes the first value of MAG Series
            o_mag_bin = get_norm_mag(source_df['MEDIAN_MAG_ISO'].iloc[0])
            # Takes the first value of PM Series
            o_pm_norm = get_norm_speed(source_df['PM'].iloc[0])
            # Takes the median value of CLASS_STAR
            o_class_star = source_df['CLASS_STAR'].iloc[0]

            for i, row in enumerate(source_df.itertuples(), 1):
                alpha = row.ALPHA_J2000
                delta = row.DELTA_J2000
                dither_n = get_dither(int(row.CATALOG_NUMBER))
                # Checks object type
                keys = ['RA', 'DEC']  # Catalogue version 2
                test_sso = check_source(input_df[dither_n]['SSOs'],
                                        alpha, delta, keys)
                if test_sso.empty is not True:
                    pass
                else:
                    self.false_positives[dither_n]['RA'].append(alpha)
                    self.false_positives[dither_n]['DEC'].append(delta)
                    self.false_positives[dither_n]['MAG'].append(o_mag_bin)
                    self.false_positives[dither_n]['PM'].append(o_pm_norm)
                    self.false_positives[dither_n]['CLASS'].append(o_class_star)
                    object_type = get_object(alpha, delta, self.input_d)
                    print(o_pm_norm, o_class_star, object_type)
                    self.false_positives[dither_n]['OBJECT'].append(object_type)

        # Regions creation
        for dither_ in self.false_positives.keys():
            alpha_list = self.false_positives[dither_]['RA']
            alpha_serie = Series(alpha_list, name='ALPHA_J2000')
            delta_list = self.false_positives[dither_]['DEC']
            delta_serie = Series(delta_list, name='DELTA_J2000')
            mag_list = self.false_positives[dither_]['MAG']
            mag_serie = Series(mag_list, name='MAG_ISO')
            pm_list = self.false_positives[dither_]['PM']
            pm_serie = Series(pm_list, name='PM')
            object_list = self.false_positives[dither_]['OBJECT']
            object_serie = Series(object_list, name='OBJECT')

            output = concat([alpha_serie, delta_serie, pm_serie], axis=1)
            for pm_ in [0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0]:
                output_pm = output[output['PM'].isin([pm_])]
                cat_name = 'false_positives/false_{}_{}.reg'.format(pm_,
                                                                    dither_)
                output_pm.to_csv(cat_name, index=False, header=False, sep=" ")

        # Catalogue creation
        dither_total_list = []
        alpha_total_list = []
        delta_total_list = []
        mag_total_list = []
        pm_total_list = []
        class_total_list = []
        object_total_list = []
        for dither_ in self.false_positives.keys():
            alpha_list = self.false_positives[dither_]['RA']
            for alpha_ in alpha_list:
                dither_total_list.append(dither_)
                alpha_total_list.append(alpha_)
            delta_list = self.false_positives[dither_]['DEC']
            for delta_ in delta_list:
                delta_total_list.append(delta_)
            mag_list = self.false_positives[dither_]['MAG']
            for mag_ in mag_list:
                mag_total_list.append(mag_)
            pm_list = self.false_positives[dither_]['PM']
            for pm_ in pm_list:
                pm_total_list.append(pm_)
            class_list = self.false_positives[dither_]['CLASS']
            for class_ in class_list:
                class_total_list.append(class_)
            object_list = self.false_positives[dither_]['OBJECT']
            for object_ in object_list:
                object_total_list.append(object_)

        dither_serie = Series(dither_total_list, name='DITHER')
        alpha_serie = Series(alpha_total_list, name='ALPHA_J2000')
        delta_serie = Series(delta_total_list, name='DELTA_J2000')
        mag_serie = Series(mag_total_list, name='MAG_ISO')
        pm_serie = Series(pm_total_list, name='PM')
        class_serie = Series(class_total_list, name='CLASS_STAR')
        object_serie = Series(object_total_list, name='OBJECT')

        output = concat([dither_serie, alpha_serie, delta_serie,
                         mag_serie, pm_serie, class_serie,
                         object_serie], axis=1)
        for pm_ in [0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0]:
                output_pm = output[output['PM'].isin([pm_])]
                cat_name = 'false_positives/false_{}.csv'.format(pm_)
                output_pm.to_csv(cat_name)


if __name__ == "__main__":
    FalsePositivesScampPerformance()
