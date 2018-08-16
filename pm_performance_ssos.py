#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Python script for performance

Versions:
- 0.1

Todo:
    * Improve log messages

"""
from numpy import mean, std
from pandas import concat, read_csv, Series

from misc import all_same, extract_settings
from regions import Create_regions


class PMPerformanceSSOs:
    def __init__(self, logger, mag, sex_cf, scmp_cf):
        """

        """
        self.logger = logger
        self.filter_p_number = 2
        self.prfs_d = extract_settings()

        self.mag = mag
        self.scmp_cf = scmp_cf
        self.sex_cf = sex_cf

        self.save = True

        self.logger.debug('Gets dispersion of output values of proper motion')
        self.data_d = self.creates_output_dict()
        self.check_pm_distribution()

    def check_source(self, catalog_n, o_cat, i_alpha, i_delta):
        """

        :param catalog_n:
        :param o_cat:
        :param i_alpha:
        :param i_delta:
        :return:
        """
        o_df = o_cat[o_cat['CATALOG_NUMBER'].isin([catalog_n])]
        o_df = o_df[o_df['ALPHA_J2000'] + self.prfs_d['tolerance'] > i_alpha]
        o_df = o_df[i_alpha > o_df['ALPHA_J2000'] - self.prfs_d['tolerance']]
        o_df = o_df[o_df['DELTA_J2000'] + self.prfs_d['tolerance'] > i_delta]
        o_df = o_df[i_delta > o_df['DELTA_J2000'] - self.prfs_d['tolerance']]

        return o_df

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
            cat_name = '{}/Cat_20-21_d{}'.format(cat_location, dither)
            input_dict[dither] = '{}.dat'.format(cat_name)
        input_dict = Create_regions(input_dict).check_ssos(self.mag, True)

        return input_dict

    def creates_input_df(self, input_dict):
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

        if self.save:
            input_df.to_csv('inputs.csv')

        return input_df

    def creates_output_dict(self):
        """

        :return: data_d
        """
        data_d = {}
        for pm_ in self.prfs_d['pms']:
            data_d[pm_] = []

        return data_d

    def gets_filtered_catalog(self):
        """

        :return: filtered_cat, filtered catalog
        """
        filter_n = 'filt_{}_{}_{}.csv'.format(self.scmp_cf, self.mag,
                                              self.filter_p_number)
        filter_o_n = '{}/{}/{}/{}/{}'.format(self.prfs_d['filter_dir'],
                                             self.mag, self.sex_cf,
                                             self.scmp_cf, filter_n)
        self.logger.debug('opens filtered catalog {}'.format(filter_n))
        # Cross with filtered data - Opens datafile
        filtered_cat = read_csv('{}'.format(filter_o_n), index_col=0)

        return filtered_cat

    def check_pm_distribution(self):
        """

        :return:
        """
        # Creates an input dictionary with all input sources
        self.logger.debug('Scamp configuration: {}'.format(self.scmp_cf))
        self.logger.debug('Sextractor configuration: {}'.format(self.sex_cf))

        input_dict = self.creates_input_dict()
        input_df = self.creates_input_df(input_dict)

        # Open particular file!
        filter_cat = self.gets_filtered_catalog()

        # Gets unique sources from input data
        unique_sources = list(set(input_df['source'].tolist()))
        sources_n = len(unique_sources)
        self.logger.debug('input sources to be analysed {}'.format(sources_n))
        # Loops over input data (Luca's catalog)
        for idx_source, source_ in enumerate(unique_sources):
            tmp_dict = {'SOURCE_NUMBER': [], 'DF': []}
            i_pm = 0
            # Gets associated data in input catalog
            cat_df = input_df[input_df['source'].isin([source_])]
            # Iterate over each detection of each source
            for i, row in enumerate(cat_df.itertuples(), 1):
                catalog_n = row.catalog
                i_pm = row.pm_values
                i_alpha = row.alpha_j2000
                i_delta = row.delta_j2000

                # Checks if there is a source closed to input one
                o_df = self.check_source(catalog_n, filter_cat,
                                         i_alpha, i_delta)

                # If there is one saves data from input data
                if o_df.empty is not True and o_df['PM'].size == 1:
                    scmp_source = o_df['SOURCE_NUMBER'].iloc[0]
                    tmp_dict['SOURCE_NUMBER'].append(scmp_source)
                    tmp_dict['DF'].append(o_df)
                else:
                    tmp_dict['SOURCE_NUMBER'].append('False')

            flag_detection, sources_number = all_same(tmp_dict['SOURCE_NUMBER'])

            if flag_detection:
                # Checks if there is any false detection
                df = tmp_dict['DF'][0]
                o_pm = df['PM'].iloc[0]
                self.data_d[i_pm].append(o_pm)

        mean_l = []
        std_l = []
        detected_l = []

        for key_ in self.prfs_d['pms']:
            if len(self.data_d[key_]) != 0:
                mean_l.append(mean(self.data_d[key_]))
                std_l.append(std(self.data_d[key_]))
                detected_l.append(len(self.data_d[key_]))
            else:
                mean_l.append('nan')
                std_l.append('nan')
                detected_l.append(0)

        pm_s = Series(self.prfs_d['pms'], name='pm')
        mean_s = Series(mean_l, name='mean')
        std_s = Series(std_l, name='std')
        detected_s = Series(detected_l, name='detected')

        stats_df = concat([pm_s, mean_s, std_s, detected_s], axis=1)
        stats_df.to_csv('pm_distribution.csv')


class SlowPMPerformanceSSOs:
    def __init__(self, logger, mag, sex_cf, scmp_cf):
        """

        """
        self.logger = logger
        self.filter_p_number = 4
        self.prfs_d = extract_settings()

        self.mag = mag
        self.scmp_cf = scmp_cf
        self.sex_cf = sex_cf

        self.save = True

        self.logger.debug('Gets output values of slow objects')
        self.data_d = self.creates_output_dict()
        self.a_image_d = self.creates_output_dict()
        self.b_image_d = self.creates_output_dict()
        self.class_star_d = self.creates_output_dict()
        self.check_pm_distribution()

    def check_source(self, catalog_n, o_cat, i_alpha, i_delta):
        """

        :param catalog_n:
        :param o_cat:
        :param i_alpha:
        :param i_delta:
        :return:
        """
        o_df = o_cat[o_cat['CATALOG_NUMBER'].isin([catalog_n])]
        o_df = o_df[o_df['ALPHA_J2000'] + self.prfs_d['tolerance'] > i_alpha]
        o_df = o_df[i_alpha > o_df['ALPHA_J2000'] - self.prfs_d['tolerance']]
        o_df = o_df[o_df['DELTA_J2000'] + self.prfs_d['tolerance'] > i_delta]
        o_df = o_df[i_delta > o_df['DELTA_J2000'] - self.prfs_d['tolerance']]

        return o_df

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
        input_dict = Create_regions(input_dict).check_ssos(self.mag, True)

        return input_dict

    def creates_input_df(self, input_dict):
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

        if self.save:
            input_df.to_csv('inputs.csv')

        return input_df

    def creates_output_dict(self):
        """

        :return: data_d
        """
        data_d = {}
        for pm_ in self.prfs_d['pms']:
            data_d[pm_] = []

        return data_d

    def gets_filtered_catalog(self):
        """

        :return: filtered_cat, filtered catalog
        """
        filter_n = 'filt_{}_{}_{}.csv'.format(self.scmp_cf, self.mag,
                                              self.filter_p_number)
        filter_o_n = '{}/{}/{}/{}/{}'.format(self.prfs_d['filter_dir'],
                                             self.mag, self.sex_cf,
                                             self.scmp_cf, filter_n)
        self.logger.debug('Opens filtered catalog {}'.format(filter_n))
        # Cross with filtered data - Opens datafile
        filtered_cat = read_csv('{}'.format(filter_o_n), index_col=0)

        return filtered_cat

    def check_pm_distribution(self):
        """

        :return:
        """
        # Creates an input dictionary with all input sources
        self.logger.debug('Scamp configuration: {}'.format(self.scmp_cf))
        self.logger.debug('Sextractor configuration: {}'.format(self.sex_cf))

        input_dict = self.creates_input_dict()
        input_df = self.creates_input_df(input_dict)

        # Open particular file!
        filter_cat = self.gets_filtered_catalog()

        # Gets unique sources from input data
        unique_sources = list(set(input_df['source'].tolist()))
        sources_n = len(unique_sources)
        self.logger.debug('input sources to be analysed {}'.format(sources_n))

        pm_dict = {'SOURCE': [], 'CATALOG_N': [], 'ALPHA_J2000': [],
                   'DELTA_J2000': [], 'CLASS_STAR': []}

        # Loops over input data (Luca's catalog)
        for idx_source, source_ in enumerate(unique_sources):
            tmp_dict = {'SOURCE_NUMBER': [], 'DF': []}
            i_pm = 0
            # Gets associated data in input catalog
            cat_df = input_df[input_df['source'].isin([source_])]
            # Iterate over each detection of each source
            for i, row in enumerate(cat_df.itertuples(), 1):
                catalog_n = row.catalog
                i_pm = row.pm_values
                i_alpha = row.alpha_j2000
                i_delta = row.delta_j2000

                # Checks if there is a source closed to input one
                o_df = self.check_source(catalog_n, filter_cat,
                                         i_alpha, i_delta)

                # If there is one saves data from input data
                if o_df.empty is not True and o_df['PM'].size == 1:
                    scmp_source = o_df['SOURCE_NUMBER'].iloc[0]
                    tmp_dict['SOURCE_NUMBER'].append(scmp_source)
                    tmp_dict['DF'].append(o_df)
                else:
                    tmp_dict['SOURCE_NUMBER'].append('False')

            flag_detection, sources_number = all_same(tmp_dict['SOURCE_NUMBER'])

            if flag_detection:
                # Checks if there is any false detection
                df = tmp_dict['DF'][0]
                a_image = df['A_IMAGE'].iloc[0]
                self.a_image_d[i_pm].append(a_image)
                b_image = df['B_IMAGE'].iloc[0]
                self.b_image_d[i_pm].append(b_image)
                class_star = df['CLASS_STAR'].iloc[0]
                self.class_star_d[i_pm].append(class_star)
                o_pm = df['PM'].iloc[0]
                self.data_d[i_pm].append(o_pm)

                if i_pm == 0.001:
                    full_df = tmp_dict['DF']
                    for df_ in full_df:
                        source = df_['SOURCE_NUMBER'].iloc[0]
                        pm_dict['SOURCE'].append(source)
                        catalog_n = df_['CATALOG_NUMBER'].iloc[0]
                        pm_dict['CATALOG_N'].append(catalog_n)
                        alpha = df_['ALPHA_J2000'].iloc[0]
                        pm_dict['ALPHA_J2000'].append(alpha)
                        delta = df_['DELTA_J2000'].iloc[0]
                        pm_dict['DELTA_J2000'].append(delta)
                        class_star = df_['CLASS_STAR'].iloc[0]
                        pm_dict['CLASS_STAR'].append(class_star)

        from pandas import DataFrame
        pm_df = DataFrame(pm_dict)
        pm_df.to_csv('pm_regions.csv')

        mean_l = []
        std_l = []
        mean_a_image_l = []
        std_a_image_l = []
        mean_b_image_l = []
        std_b_image_l = []
        mean_class_star_l = []
        std_class_star_l = []
        detected_l = []

        for key_ in self.prfs_d['pms']:
            if len(self.data_d[key_]) != 0:
                # print(key_)
                # print(self.a_image_d[key_])
                # print(self.b_image_d[key_])
                # print(self.class_star_d[key_])
                mean_l.append(mean(self.data_d[key_]))
                std_l.append(std(self.data_d[key_]))
                mean_a_image_l.append(mean(self.a_image_d[key_]))
                std_a_image_l.append(std(self.a_image_d[key_]))
                mean_b_image_l.append(mean(self.b_image_d[key_]))
                std_b_image_l.append(std(self.b_image_d[key_]))
                mean_class_star_l.append(mean(self.class_star_d[key_]))
                std_class_star_l.append(std(self.class_star_d[key_]))
                detected_l.append(len(self.data_d[key_]))
            else:
                mean_l.append('nan')
                std_l.append('nan')
                mean_a_image_l.append('nan')
                std_a_image_l.append('nan')
                mean_b_image_l.append('nan')
                std_b_image_l.append('nan')
                mean_class_star_l.append('nan')
                std_class_star_l.append('nan')
                detected_l.append(0)

        pm_s = Series(self.prfs_d['pms'], name='pm')
        mean_s = Series(mean_l, name='mean')
        std_s = Series(std_l, name='std')
        mean_a_image_s = Series(mean_a_image_l, name='a_image-mean')
        std_a_image_s = Series(std_a_image_l, name='a_image-std')
        mean_b_image_s = Series(mean_b_image_l, name='b_image-mean')
        std_b_image_s = Series(std_b_image_l, name='b_image-std')
        mean_class_star_s = Series(mean_class_star_l, name='class_star-mean')
        std_class_star_s = Series(std_class_star_l, name='class_star-std')
        detected_s = Series(detected_l, name='detected')

        stats_df = concat([pm_s, mean_s, std_s, mean_a_image_s, std_a_image_s,
                           mean_b_image_s, std_b_image_s, mean_class_star_s,
                           std_class_star_s, detected_s], axis=1)
        stats_df.to_csv('slow_pm_distribution_{}.csv'.format(self.filter_p_number))
