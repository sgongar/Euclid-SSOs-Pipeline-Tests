#!/usr/bin/python
# -*- coding: utf-8 -*-

"""


Todo:
    *

*GNU Terry Pratchett*
"""

from math import cos, sin, radians

from astropy.io import fits
from astropy.wcs import WCS
from numpy import genfromtxt, float64
from pandas import concat, Series

from misc import extract_settings_luca
from misc_fits import get_fits_d, get_fits_limits


__author__ = "Samuel Gongora-Garcia"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Gongora-Garcia"]
__version__ = "0.1"
__maintainer__ = "Samuel Gongora-Garcia"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


class Create_regions:

    def __init__(self, input_catalog):
        """

        @param input_catalog:
        """
        self.input_catalogue = input_catalog
        self.prfs_d = extract_settings_luca()

    def check_stars(self, mag, complete):
        """

        :param mag:
        :param complete:
        :return:
        """
        input_d = self.input_catalogue

        for dither_ in input_d.keys():
            catalog = genfromtxt(input_d[dither_])

            list_x = catalog[:, 0]
            list_y = catalog[:, 1]
            list_mag = catalog[:, 2]
            list_pm = catalog[:, 3]

            x_values = []
            y_values = []

            stars = range(0, 12471, 1)
            indexes = sorted(stars)

            s1 = Series(list_x, name='X_IMAGE', dtype=float64)
            s2 = Series(list_y, name='Y_IMAGE', dtype=float64)
            s3 = Series(list_mag, name='MAG_VALUES', dtype=float64)
            s4 = Series(list_pm, name='PM_INPUT', dtype=float64)

            sources_df = concat([s1, s2, s3, s4], axis=1)
            sources_df = sources_df.iloc[indexes, :]

            ccd_loc = 'mag_{}_CCD_x0_y0_d1.fits'.format(mag)
            fits_loc = '{}/{}/CCDs/{}'.format(self.prfs_d['fits_dir'], mag,
                                              ccd_loc)
            hdulist = fits.open(fits_loc)
            w = WCS(hdulist[0].header)

            regions_list = []
            for source_num in range(sources_df['X_IMAGE'].as_matrix().size):
                x_value = sources_df['X_IMAGE'].as_matrix()[source_num]
                y_value = sources_df['Y_IMAGE'].as_matrix()[source_num]
                regions_list.append([x_value, y_value])

            input_regions = w.all_pix2world(regions_list, 1)

            alpha_list = []
            delta_list = []
            for idx, regions in enumerate(input_regions):
                alpha_list.append(regions[0])
                delta_list.append(regions[1])
                x_values.append(regions_list[idx][0])
                y_values.append(regions_list[idx][1])

            fits_files_all = get_fits_d(mag_=mag, dither=dither_)

            fits_dict = {}
            for fits_ in fits_files_all:
                fits_loc = '{}/{}/CCDs'.format(self.prfs_d['fits_dir'], mag)
                fits_file = '{}/{}'.format(fits_loc, fits_)
                ccd = fits_[-13:-8]
                fits_dict[ccd] = get_fits_limits(fits_file)

            i = 0
            ccd_list = []
            for alpha_, delta_ in zip(alpha_list, delta_list):
                i += 1
                flag = True
                for key_ in fits_dict.keys():
                    below_ra = fits_dict[key_]['below_ra']
                    above_ra = fits_dict[key_]['above_ra']
                    below_dec = fits_dict[key_]['below_dec']
                    above_dec = fits_dict[key_]['above_dec']
                    alpha_comp = below_ra < alpha_ < above_ra
                    delta_comp = below_dec < delta_ < above_dec
                    if alpha_comp and delta_comp:
                        ccd_list.append(key_)
                        flag = False
                if flag:
                    ccd_list.append('False')

            # Creates a list for all sources
            source_list = range(0, len(alpha_list), 1)
            # Populates a list for all sources with the dither number
            dither_list = []
            for dither_idx in range(len(alpha_list)):
                dither_list.append(dither_)

            cats = [['x0_y0', 1, 1], ['x0_y0', 2, 2], ['x0_y0', 3, 3],
                    ['x0_y0', 4, 4], ['x0_y1', 1, 5], ['x0_y1', 2, 6],
                    ['x0_y1', 3, 7], ['x0_y1', 4, 8], ['x0_y2', 1, 9],
                    ['x0_y2', 2, 10], ['x0_y2', 3, 11], ['x0_y2', 4, 12],
                    ['x1_y0', 1, 13], ['x1_y0', 2, 14], ['x1_y0', 3, 15],
                    ['x1_y0', 4, 16], ['x1_y1', 1, 17], ['x1_y1', 2, 18],
                    ['x1_y1', 3, 19], ['x1_y1', 4, 20], ['x1_y2', 1, 21],
                    ['x1_y2', 2, 22], ['x1_y2', 3, 23], ['x1_y2', 4, 24],
                    ['x2_y0', 1, 25], ['x2_y0', 2, 26], ['x2_y0', 3, 27],
                    ['x2_y0', 4, 28], ['x2_y1', 1, 29], ['x2_y1', 2, 30],
                    ['x2_y1', 3, 31], ['x2_y1', 4, 32], ['x2_y2', 1, 33],
                    ['x2_y2', 2, 34], ['x2_y2', 3, 35], ['x2_y2', 4, 36]]

            cats_list = []
            for dither_, CCD_ in zip(dither_list, ccd_list):
                flag = True
                for cat_ in cats:
                    if cat_[1] == dither_ and cat_[0] == CCD_:
                        flag = False
                        cats_list.append(cat_[2])

                if flag:
                    cats_list.append(False)

            # Creates a serie of Pandas Series
            series_source = Series(source_list, name='source')
            series_alpha_j2000 = Series(alpha_list, name='alpha_j2000')
            series_delta_j2000 = Series(delta_list, name='delta_j2000')
            series_mag = Series(sources_df['MAG_VALUES'].tolist(),
                                name='mag_values')
            series_pm = Series(sources_df['PM_INPUT'].tolist(),
                               name='pm_values')
            series_dither = Series(dither_list, name='dither_values')
            series_ccd = Series(ccd_list, name='CCD')
            series_cat = Series(cats_list, name='catalog')

            if complete:
                sources_df = concat([series_source, series_cat,
                                     series_alpha_j2000, series_delta_j2000,
                                     series_mag, series_pm, series_dither,
                                     series_ccd], axis=1)
                sources_df = sources_df[~sources_df['CCD'].isin(['False'])]
            else:
                sources_df = concat([series_alpha_j2000,
                                     series_delta_j2000], axis=1)

            input_d[dither_] = sources_df

        return input_d

    def check_galaxies(self, mag, complete):
        """

        :param mag:
        :param complete:
        :return:
        """
        input_d = self.input_catalogue

        for dither_ in input_d.keys():
            catalog = genfromtxt(input_d[dither_])

            list_x = catalog[:, 0]
            list_y = catalog[:, 1]
            list_mag = catalog[:, 2]
            list_pm = catalog[:, 3]

            x_values = []
            y_values = []

            galaxies = range(12471, 134892, 1)
            indexes = sorted(galaxies)

            s1 = Series(list_x, name='X_IMAGE', dtype=float64)
            s2 = Series(list_y, name='Y_IMAGE', dtype=float64)
            s3 = Series(list_mag, name='MAG_VALUES', dtype=float64)
            s4 = Series(list_pm, name='PM_INPUT', dtype=float64)

            sources_df = concat([s1, s2, s3, s4], axis=1)
            sources_df = sources_df.iloc[indexes, :]

            ccd_loc = 'mag_{}_CCD_x0_y0_d1.fits'.format(mag)
            fits_loc = '{}/{}/CCDs/{}'.format(self.prfs_d['fits_dir'], mag,
                                              ccd_loc)
            hdulist = fits.open(fits_loc)
            w = WCS(hdulist[0].header)

            regions_list = []
            for source_num in range(sources_df['X_IMAGE'].as_matrix().size):
                x_value = sources_df['X_IMAGE'].as_matrix()[source_num]
                y_value = sources_df['Y_IMAGE'].as_matrix()[source_num]
                regions_list.append([x_value, y_value])

            input_regions = w.all_pix2world(regions_list, 1)

            alpha_list = []
            delta_list = []
            for idx, regions in enumerate(input_regions):
                alpha_list.append(regions[0])
                delta_list.append(regions[1])
                x_values.append(regions_list[idx][0])
                y_values.append(regions_list[idx][1])

            fits_files_all = get_fits_d(mag_=mag, dither=dither_)

            fits_dict = {}
            for fits_ in fits_files_all:
                fits_loc = '{}/{}/CCDs'.format(self.prfs_d['fits_dir'], mag)
                fits_file = '{}/{}'.format(fits_loc, fits_)
                ccd = fits_[-13:-8]
                fits_dict[ccd] = get_fits_limits(fits_file)

            i = 0
            CCD_list = []
            for alpha_, delta_ in zip(alpha_list, delta_list):
                i += 1
                flag = True
                for key_ in fits_dict.keys():
                    below_ra = fits_dict[key_]['below_ra']
                    above_ra = fits_dict[key_]['above_ra']
                    below_dec = fits_dict[key_]['below_dec']
                    above_dec = fits_dict[key_]['above_dec']
                    alpha_comp = below_ra < alpha_ < above_ra
                    delta_comp = below_dec < delta_ < above_dec
                    if alpha_comp and delta_comp:
                        CCD_list.append(key_)
                        flag = False
                if flag:
                    CCD_list.append('False')

            # Creates a list for all sources
            source_list = range(0, len(alpha_list), 1)
            # Populates a list for all sources with the dither number
            dither_list = []
            for dither_idx in range(len(alpha_list)):
                dither_list.append(dither_)

            cats = [['x0_y0', 1, 1], ['x0_y0', 2, 2], ['x0_y0', 3, 3],
                    ['x0_y0', 4, 4], ['x0_y1', 1, 5], ['x0_y1', 2, 6],
                    ['x0_y1', 3, 7], ['x0_y1', 4, 8], ['x0_y2', 1, 9],
                    ['x0_y2', 2, 10], ['x0_y2', 3, 11], ['x0_y2', 4, 12],
                    ['x1_y0', 1, 13], ['x1_y0', 2, 14], ['x1_y0', 3, 15],
                    ['x1_y0', 4, 16], ['x1_y1', 1, 17], ['x1_y1', 2, 18],
                    ['x1_y1', 3, 19], ['x1_y1', 4, 20], ['x1_y2', 1, 21],
                    ['x1_y2', 2, 22], ['x1_y2', 3, 23], ['x1_y2', 4, 24],
                    ['x2_y0', 1, 25], ['x2_y0', 2, 26], ['x2_y0', 3, 27],
                    ['x2_y0', 4, 28], ['x2_y1', 1, 29], ['x2_y1', 2, 30],
                    ['x2_y1', 3, 31], ['x2_y1', 4, 32], ['x2_y2', 1, 33],
                    ['x2_y2', 2, 34], ['x2_y2', 3, 35], ['x2_y2', 4, 36]]

            cats_list = []
            for dither_, CCD_ in zip(dither_list, CCD_list):
                flag = True
                for cat_ in cats:
                    if cat_[1] == dither_ and cat_[0] == CCD_:
                        flag = False
                        cats_list.append(cat_[2])

                if flag:
                    cats_list.append(False)

            # Creates a series of Pandas Series
            series_source = Series(source_list, name='source')
            series_alpha_j2000 = Series(alpha_list, name='alpha_j2000')
            series_delta_j2000 = Series(delta_list, name='delta_j2000')
            series_mag = Series(sources_df['MAG_VALUES'].tolist(),
                                name='mag_values')
            series_pm = Series(sources_df['PM_INPUT'].tolist(),
                               name='pm_values')
            series_dither = Series(dither_list, name='dither_values')
            series_ccd = Series(CCD_list, name='CCD')
            series_cat = Series(cats_list, name='catalog')

            if complete:
                sources_df = concat([series_source, series_cat,
                                     series_alpha_j2000, series_delta_j2000,
                                     series_mag, series_pm, series_dither,
                                     series_ccd], axis=1)
                sources_df = sources_df[~sources_df['CCD'].isin(['False'])]
            else:
                sources_df = concat([series_alpha_j2000,
                                     series_delta_j2000], axis=1)

            input_d[dither_] = sources_df

        return input_d

    def check_ssos(self, mag_, complete):
        """

        :param mag_:
        :param complete:
        :return:
        """
        input_d = self.input_catalogue

        for dither_ in input_d.keys():
            catalog = genfromtxt(input_d[dither_])

            list_x = catalog[:, 0]
            list_y = catalog[:, 1]
            list_mag = catalog[:, 2]
            list_pm = catalog[:, 3]
            list_angle = catalog[:, 4]
            total_length = len(list_x)
            list_alpha_pm = range(0, total_length, 1)
            list_delta_pm = range(0, total_length, 1)

            x_values = []
            y_values = []

            speed_0_001 = range(self.prfs_d['first_sso'], 137378, 73)
            speed_0_003 = range(self.prfs_d['first_sso'] + 10, 137388, 73)
            speed_0_01 = range(self.prfs_d['first_sso'] + 20, 137398, 73)
            speed_0_03 = range(self.prfs_d['first_sso'] + 30, 137408, 73)
            speed_0_1 = range(self.prfs_d['first_sso'] + 40, 137418, 73)
            speed_0_3 = range(self.prfs_d['first_sso'] + 50, 137428, 73)
            speed_1 = range(self.prfs_d['first_sso'] + 60, 137438, 73)
            speed_3 = range(self.prfs_d['first_sso'] + 67, 137445, 73)
            speed_10 = range(self.prfs_d['first_sso'] + 68, 137446, 73)
            speed_30 = range(self.prfs_d['first_sso'] + 69, 137447, 73)

            # Gets proper angles
            list_angle = [(angle_ + 180) for angle_ in list_angle]

            # Gets proper motion
            for index in speed_0_001:
                list_mag[index] = list_mag[index] - 2.5
                list_pm[index] = 0.001
                list_alpha_pm[index] = 0.001 * cos(radians(list_angle[index]))
                list_delta_pm[index] = 0.001 * sin(radians(list_angle[index]))
            for index in speed_0_003:
                list_mag[index] = list_mag[index] - 2.5
                list_pm[index] = 0.003
                list_alpha_pm[index] = 0.003 * cos(radians(list_angle[index]))
                list_delta_pm[index] = 0.003 * sin(radians(list_angle[index]))
            for index in speed_0_01:
                list_mag[index] = list_mag[index] - 2.5
                list_pm[index] = 0.01
                list_alpha_pm[index] = 0.01 * cos(radians(list_angle[index]))
                list_delta_pm[index] = 0.01 * sin(radians(list_angle[index]))
            for index in speed_0_03:
                list_mag[index] = list_mag[index] - 2.5
                list_pm[index] = 0.03
                list_alpha_pm[index] = 0.03 * cos(radians(list_angle[index]))
                list_delta_pm[index] = 0.03 * sin(radians(list_angle[index]))
            for index in speed_0_1:
                list_mag[index] = list_mag[index] - 2.5
                list_pm[index] = 0.1
                list_alpha_pm[index] = 0.1 * cos(radians(list_angle[index]))
                list_delta_pm[index] = 0.1 * sin(radians(list_angle[index]))
            for index in speed_0_3:
                list_mag[index] = list_mag[index] - 2.5
                list_pm[index] = 0.3
                list_alpha_pm[index] = 0.3 * cos(radians(list_angle[index]))
                list_delta_pm[index] = 0.3 * sin(radians(list_angle[index]))
            for index in speed_1:
                list_mag[index] = list_mag[index] - 2.5
                list_pm[index] = 1
                list_alpha_pm[index] = 1 * cos(radians(list_angle[index]))
                list_delta_pm[index] = 1 * sin(radians(list_angle[index]))
            for index in speed_3:
                list_pm[index] = list_pm[index] - 1000
                list_alpha_pm[index] = 3 * cos(radians(list_angle[index]))
                list_delta_pm[index] = 3 * sin(radians(list_angle[index]))
            for index in speed_10:
                list_pm[index] = list_pm[index] - 1000
                list_alpha_pm[index] = 10 * cos(radians(list_angle[index]))
                list_delta_pm[index] = 10 * sin(radians(list_angle[index]))
            for index in speed_30:
                list_pm[index] = list_pm[index] - 1000
                list_alpha_pm[index] = 30 * cos(radians(list_angle[index]))
                list_delta_pm[index] = 30 * sin(radians(list_angle[index]))

            indexes = (speed_0_001 + speed_0_003 + speed_0_01 + speed_0_03 +
                       speed_0_1 + speed_0_3 + speed_1 + speed_3 + speed_10 +
                       speed_30)
            indexes = sorted(indexes)

            s1 = Series(list_x, name='X_IMAGE', dtype=float64)
            s2 = Series(list_y, name='Y_IMAGE', dtype=float64)
            s3 = Series(list_mag, name='MAG_VALUES', dtype=float64)
            s4 = Series(list_pm, name='PM_INPUT', dtype=float64)
            s5 = Series(list_alpha_pm, name='PM_ALPHA', dtype=float64)
            s6 = Series(list_delta_pm, name='PM_DELTA', dtype=float64)
            s7 = Series(list_angle, name='ANGLE_INPUT', dtype=float64)

            sources_df = concat([s1, s2, s3, s4, s5, s6, s7], axis=1)
            sources_df = sources_df.iloc[indexes, :]

            ccd_loc = 'mag_{}_CCD_x0_y0_d1.fits'.format(mag_)
            fits_loc = '{}/{}/CCDs/{}'.format(self.prfs_d['fits_dir'], mag_,
                                              ccd_loc)
            hdulist = fits.open(fits_loc)
            w = WCS(hdulist[0].header)

            regions_list = []
            for source_num in range(sources_df['X_IMAGE'].as_matrix().size):
                x_value = sources_df['X_IMAGE'].as_matrix()[source_num]
                y_value = sources_df['Y_IMAGE'].as_matrix()[source_num]
                regions_list.append([x_value, y_value])

            input_regions = w.all_pix2world(regions_list, 1)

            alpha_list = []
            delta_list = []
            for idx, regions in enumerate(input_regions):
                alpha_list.append(regions[0])
                delta_list.append(regions[1])
                x_values.append(regions_list[idx][0])
                y_values.append(regions_list[idx][1])

            fits_files_all = get_fits_d(mag_, dither=dither_)

            fits_dict = {}
            for fits_ in fits_files_all:
                CCD = fits_[-13:-8]
                fits_file = '{}/{}/CCDs/{}'.format(self.prfs_d['fits_dir'],
                                                   mag_, fits_)
                fits_dict[CCD] = get_fits_limits(fits_file)

            i = 0
            CCD_list = []
            for alpha_, delta_ in zip(alpha_list, delta_list):
                i += 1
                flag = True
                for key_ in fits_dict.keys():
                    below_ra = fits_dict[key_]['below_ra']
                    above_ra = fits_dict[key_]['above_ra']
                    below_dec = fits_dict[key_]['below_dec']
                    above_dec = fits_dict[key_]['above_dec']
                    alpha_comp = below_ra < alpha_ < above_ra
                    delta_comp = below_dec < delta_ < above_dec
                    if alpha_comp and delta_comp:
                        CCD_list.append(key_)
                        flag = False
                if flag:
                    CCD_list.append('False')

            # Creates a list for all sources
            source_list = range(0, len(alpha_list), 1)
            # Populates a list for all sources with the dither number
            dither_list = []
            for dither_idx in range(len(alpha_list)):
                dither_list.append(dither_)

            cats = [['x0_y0', 1, 1], ['x0_y0', 2, 2], ['x0_y0', 3, 3],
                    ['x0_y0', 4, 4], ['x0_y1', 1, 5], ['x0_y1', 2, 6],
                    ['x0_y1', 3, 7], ['x0_y1', 4, 8], ['x0_y2', 1, 9],
                    ['x0_y2', 2, 10], ['x0_y2', 3, 11], ['x0_y2', 4, 12],
                    ['x1_y0', 1, 13], ['x1_y0', 2, 14], ['x1_y0', 3, 15],
                    ['x1_y0', 4, 16], ['x1_y1', 1, 17], ['x1_y1', 2, 18],
                    ['x1_y1', 3, 19], ['x1_y1', 4, 20], ['x1_y2', 1, 21],
                    ['x1_y2', 2, 22], ['x1_y2', 3, 23], ['x1_y2', 4, 24],
                    ['x2_y0', 1, 25], ['x2_y0', 2, 26], ['x2_y0', 3, 27],
                    ['x2_y0', 4, 28], ['x2_y1', 1, 29], ['x2_y1', 2, 30],
                    ['x2_y1', 3, 31], ['x2_y1', 4, 32], ['x2_y2', 1, 33],
                    ['x2_y2', 2, 34], ['x2_y2', 3, 35], ['x2_y2', 4, 36]]

            cats_list = []
            for dither_, CCD_ in zip(dither_list, CCD_list):
                flag = True
                for cat_ in cats:
                    if cat_[1] == dither_ and cat_[0] == CCD_:
                        flag = False
                        cats_list.append(cat_[2])

                if flag:
                    cats_list.append(False)

            # Creates a serie of Pandas Series
            series_source = Series(source_list, name='source')
            series_alpha_j2000 = Series(alpha_list, name='alpha_j2000')
            series_delta_j2000 = Series(delta_list, name='delta_j2000')
            series_mag = Series(sources_df['MAG_VALUES'].tolist(), name='mag_values')
            series_pm = Series(sources_df['PM_INPUT'].tolist(), name='pm_values')
            series_pm_alpha = Series(sources_df['PM_ALPHA'].tolist(), name='pm_alpha')
            series_pm_delta = Series(sources_df['PM_DELTA'].tolist(), name='pm_delta')
            series_angle = Series(sources_df['ANGLE_INPUT'].tolist(),
                           name='angle_values')
            series_dither = Series(dither_list, name='dither_values')
            series_ccd = Series(CCD_list, name='CCD')
            series_cat = Series(cats_list, name='catalog')

            if complete:
                sources_df = concat([series_source, series_cat,
                                     series_alpha_j2000, series_delta_j2000,
                                     series_mag, series_pm, series_pm_alpha,
                                     series_pm_delta, series_angle,
                                     series_dither, series_ccd], axis=1)
                sources_df = sources_df[~sources_df['CCD'].isin(['False'])]
            else:
                sources_df = concat([series_alpha_j2000,
                                     series_delta_j2000], axis=1)

            input_d[dither_] = sources_df

        return input_d
