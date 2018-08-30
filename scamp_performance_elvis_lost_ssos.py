# !/usr/bin/python
# -*- coding: utf-8 -*-

""" Gets
   - factors from magnitude bins

Versions:
- 0.1

Todo:
    * Improve log messages

*GNU Terry Pratchett*
"""

from astropy.io import fits
from astropy.table import Table
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from numpy import array
from pandas import DataFrame, read_csv, Series
from pyds9 import DS9
import statsmodels.api as sm

from misc import extract_settings_elvis, check_source, setting_logger


__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def creates_ssos_dict():
    """

    :return:
    """
    d = {'20-21': {0.01: 0, 0.03: 0, 0.1: 0, 0.3: 0, 1: 0, 3: 0, 10: 0, 30: 0},
         '21-22': {0.01: 0, 0.03: 0, 0.1: 0, 0.3: 0, 1: 0, 3: 0, 10: 0, 30: 0},
         '22-23': {0.01: 0, 0.03: 0, 0.1: 0, 0.3: 0, 1: 0, 3: 0, 10: 0, 30: 0},
         '23-24': {0.01: 0, 0.03: 0, 0.1: 0, 0.3: 0, 1: 0, 3: 0, 10: 0, 30: 0},
         '24-25': {0.01: 0, 0.03: 0, 0.1: 0, 0.3: 0, 1: 0, 3: 0, 10: 0, 30: 0},
         '25-26': {0.01: 0, 0.03: 0, 0.1: 0, 0.3: 0, 1: 0, 3: 0, 10: 0, 30: 0}}

    return d


def get_zoom_level(pm):
    """

    :param pm:
    :return: zoom
    """
    if pm == 30:
        zoom_level = 1
    elif pm == 10:
        zoom_level = 4
    elif pm == 3:
        zoom_level = 8
    elif pm == 1:
        zoom_level = 8
    elif pm == 0.3:
        zoom_level = 8
    elif pm == 0.1:
        zoom_level = 8
    elif pm == 0.03:
        zoom_level = 8
    elif pm == 0.01:
        zoom_level = 8
    else:
        zoom_level = 2

    return zoom_level


def get_chisquared(ra, error_ra, dec, error_dec, epoch):
    """

    :param ra:
    :param error_ra:
    :param dec:
    :param error_dec:
    :param epoch:
    :return:
    """
    fittings = []
    for dimension in [ra, dec]:
        x = array(epoch)
        y = array(dimension)
        sigma = []
        if dimension == ra:
            sigma = error_ra
        if dimension == dec:
            sigma = error_dec
        edim = array([1 / var for var in sigma])
        x = sm.add_constant(x)
        # Model: y~x+c
        model = sm.WLS(y, x, weigths=edim)
        fitted = model.fit()
        fittings.append(fitted.rsquared)

    return fittings


def get_pm(catalog, source_number, catalog_number):
    """

    :param catalog:
    :param source_number:
    :param catalog_number:
    :return:
    """
    cat_df = catalog[catalog['SOURCE_NUMBER'].isin([source_number])]
    cat_df = cat_df[cat_df['CATALOG_NUMBER'].isin([catalog_number])]

    if cat_df.empty is True:
        out_pm = 'lost'
    else:
        out_pm = float(cat_df['PM'].iloc[0])

    return out_pm


def get_mag(catalog, source_number, catalog_number):
    """

    :param catalog:
    :param source_number:
    :param catalog_number:
    :return:
    """
    cat_df = catalog[catalog['SOURCE_NUMBER'].isin([source_number])]
    cat_df = cat_df[cat_df['CATALOG_NUMBER'].isin([catalog_number])]

    if cat_df.empty is True:
        out_mag = 'lost'
    else:
        out_mag = float(cat_df['MAG_ISO'].iloc[0])

    return out_mag


def gets_filtered_catalog():
    """

    :return:
    """
    prfs_d = extract_settings_elvis()

    filter_n = 'filt__9.csv'  # First one with PM
    filter_o_n = '{}/{}'.format(prfs_d['filtered'], filter_n)

    print('Opens filtered catalogue {}'.format(filter_o_n))
    filtered_cat = read_csv('{}'.format(filter_o_n), index_col=0)

    return filtered_cat


def get_cats(dither_n):
    """

    :param dither_n:
    :return: cat_list
    """
    cats_dict = {}
    for i in range(1, 5, 1):
        cats_dict[i] = []

    for cat_ccd in range(0, 144, 4):
        for cat_dither in range(1, 5, 1):
            cats_dict[cat_dither].append(cat_dither + cat_ccd)

    return cats_dict[dither_n]


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


class ScampPerformanceLostOnes:

    def __init__(self):
        """

        """
        self.filter_p_number = 9  # First one with enough data for statistics
        self.prfs_d = extract_settings_elvis()
        self.data_keys = ['IN_SOURCE', 'SOURCE_NUMBER', 'DITHER', 'IN_ALPHA',
                          'OUT_ALPHA', 'IN_DELTA', 'OUT_DELTA', 'IN_PM',
                          'IN_PM_NORM', 'OUT_PM', 'IN_MAG', 'OUT_MAG',
                          'CHI_ALPHA', 'CHI_DELTA']
        self.lost_keys = ['DITHER', 'IN_SOURCE', 'DITHER', 'IN_ALPHA', 'IN_DELTA',
                          'IN_PM', 'IN_MAG']
        self.save = True

        logger_name = 'scamp_performance'  # Set as desired
        self.logger = setting_logger(self.prfs_d, logger_name)

        scamp_cat = self.gets_scamp_catalog()  # Gets data from filtered
        self.filtered_cat = gets_filtered_catalog()
        input_df = self.gets_data() # Gets data from catalogs
        origin = 'filtered'
        lost_d, data_d = self.extract_lost(input_df, scamp_cat)

        data_df = DataFrame(data_d, columns=self.data_keys)
        data_df.to_csv('data.csv')

        """
        for dither_ in range(1, 5, 1):
            lost_df = DataFrame(lost_d, columns=self.lost_keys)
            print(lost_df.columns)
            lost_df = lost_df[lost_df['DITHER'].isin([dither_])]
            lost_df.to_csv('lost_{}.reg'.format(dither_),
                           index=False, header=False, sep=" ")
        """

    def gets_scamp_catalog(self):
        """

        :return: filtered_cat, filtered catalog
        """
        scamp_n = 'full_1.cat'
        scamp_o_n = '{}/{}'.format(self.prfs_d['output_cats'], scamp_n)

        print('Opens Scamp catalogue {}'.format(scamp_o_n))
        full_cat = fits.open(scamp_o_n)  # Opens full catalog
        full_db = Table(full_cat[2].data)  # Converts it to Astropy Table
        print('Converts full Astropy catalog to Pandas format')
        scamp_cat = full_db.to_pandas()  # Converts it to Pandas format

        return scamp_cat

    def gets_data(self):
        """ Creates an input dictionary. Each key contains SSOs' information
        for each dither.

        :return: input_dict
        """
        # # For now we only have data for dither 1
        input_df = {1: {}, 2: {}, 3: {}, 4: {}}

        for key_ in input_df.keys():
            # ssos_cat = 'cat_ssos_{}.csv'.format(key_)
            ssos_cat = 'cats/cat_clean_ssos_{}.csv'.format(key_)
            input_df[key_]['SSOs'] = read_csv(ssos_cat, index_col=0)
            stars_cat = 'tmp_stars/stars.csv'
            input_df[key_]['stars'] = read_csv(stars_cat, index_col=0)
            galaxies_cat = 'tmp_galaxies/galaxies.csv'
            input_df[key_]['galaxies'] = read_csv(galaxies_cat, index_col=0)

        return input_df

    def extract_lost(self, input_df, scamp_cat):
        """

        :param input_df:
        :param scamp_cat:
        :return:
        """
        # Creates a dictionary with the number of input sources
        ssos_d = creates_ssos_dict()
        # Creates a dictionary in order to save all de input data
        total_dict = {}
        # A list of keys contained in the input dataframes
        dict_keys = ['ABMAG', 'DEC', 'DITHER', 'IDX',
                     'RA', 'SOURCE', 'THETA', 'VEL']
        # Populates the input dictionary with lists.
        for key_ in dict_keys:
            total_dict[key_] = []

        # Merges the dataframe of each dither into a single dictionary
        for dither_ in range(1, 5, 1):
            for column_ in input_df[dither_]['SSOs'].columns:
                for value_ in input_df[dither_]['SSOs'][column_].tolist():
                    total_dict[column_].append(value_)
        # Creates a dataframe from the previous dictionary
        total_df = DataFrame(total_dict)
        # Gets the unique sources of the total dataframe
        unique_sources = list(set(total_df['SOURCE'].tolist()))

        # Data dictionary for extracted SSOs
        data_d = {}
        for key_ in self.data_keys:
            data_d[key_] = []

        # Data dictionary for lost SSOs
        lost_d = {}
        for key_ in self.lost_keys:
            lost_d[key_] = []

        # Test variables
        total_unique_sources = 0
        unique_sources_lost = 0
        unique_sources_found = 0
        total_1_0 = 0

        for idx_source_, source_ in enumerate(unique_sources):
            source_df = total_df[total_df['SOURCE'].isin([source_])]
            test = 0
            test_1_0 = 0
            for idx_column_, column_ in enumerate(source_df.itertuples(), 1):
                test += 1

                alpha = column_.RA
                delta = column_.DEC
                pm = column_.VEL
                mag = column_.ABMAG
                dither_n = column_.DITHER
                source = column_.SOURCE

                i_pm_test = get_norm_speed(pm)
                i_mag_test = get_norm_mag(mag)

                keys = ['ALPHA_J2000', 'DELTA_J2000']  # Catalogue version 2
                ccd_df = check_source(scamp_cat, alpha, delta, keys)
                cats = get_cats(dither_n)
                test_df = ccd_df[ccd_df['CATALOG_NUMBER'].isin(cats)]

                total_unique_sources += 1  # Adds one to the total number

                # Populates the dictionary of lost sources
                if test_df.empty is True:
                    unique_sources_lost += 1  # Adds one to the lost sources
                    in_source = source
                    lost_d['IN_SOURCE'].append(in_source)
                    dither = dither_n
                    lost_d['DITHER'].append(dither)
                    in_alpha = alpha
                    lost_d['IN_ALPHA'].append(in_alpha)
                    in_delta = delta
                    lost_d['IN_DELTA'].append(in_delta)
                    in_pm = pm
                    lost_d['IN_PM'].append(in_pm)
                    in_mag = mag
                    lost_d['IN_MAG'].append(in_mag)

                    output = 'NOT'

                # Populates the dictionary of extracted sources
                elif test_df.empty is False:
                    unique_sources_found += 1   # Adds one to the found sources
                    in_source = source
                    data_d['IN_SOURCE'].append(in_source)
                    source_number = int(test_df['SOURCE_NUMBER'].iloc[0])
                    data_d['SOURCE_NUMBER'].append(source_number)
                    dither = dither_n
                    data_d['DITHER'].append(dither)
                    in_alpha = alpha
                    data_d['IN_ALPHA'].append(in_alpha)
                    out_alpha = float(test_df['ALPHA_J2000'].iloc[0])
                    data_d['OUT_ALPHA'].append(out_alpha)
                    in_delta = delta
                    data_d['IN_DELTA'].append(in_delta)
                    out_delta = float(test_df['DELTA_J2000'].iloc[0])
                    data_d['OUT_DELTA'].append(out_delta)
                    in_pm = pm
                    data_d['IN_PM'].append(in_pm)
                    in_pm_norm = get_norm_speed(in_pm)
                    data_d['IN_PM_NORM'].append(in_pm_norm)
                    catalog_number = int(test_df['CATALOG_NUMBER'].iloc[0])
                    out_pm = get_pm(self.filtered_cat, source_number,
                                    catalog_number)
                    data_d['OUT_PM'].append(out_pm)
                    in_mag = mag
                    data_d['IN_MAG'].append(in_mag)
                    out_mag = get_mag(self.filtered_cat, source_number,
                                      catalog_number)
                    data_d['OUT_MAG'].append(out_mag)
                    ra = ccd_df['ALPHA_J2000'].tolist()
                    error_ra = ccd_df['ERRA_WORLD'].tolist()
                    dec = ccd_df['DELTA_J2000'].tolist()
                    error_dec = ccd_df['ERRB_WORLD'].tolist()
                    epoch = ccd_df['EPOCH'].tolist()

                    chi = get_chisquared(ra, error_ra, dec, error_dec, epoch)
                    if chi[0] >= 0.0:
                        chi_ra = chi[0]
                    else:
                        chi_ra = 'null'
                    data_d['CHI_ALPHA'].append(chi_ra)

                    if chi[1] >= 0.0:
                        chi_dec = chi[1]
                    else:
                        chi_dec = 'null'
                    data_d['CHI_DELTA'].append(chi_dec)

                    output = 'OK'

                if i_pm_test == 1.0 and i_mag_test == '23-24':
                    test_1_0 += 1

            if test_1_0 >= 3:
                total_1_0 += 1

            if test >= 3:
                mag_norm = get_norm_mag(mag)
                pm_norm = get_norm_speed(pm)
                ssos_d[mag_norm][pm_norm] += 1


        for mag_key in ['20-21', '21-22', '22-23', '23-24', '24-25', '25-26']:
            for pm_key in [0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30]:
                print(mag_key, pm_key, ssos_d[mag_key][pm_key])

        # Test variables
        print('total_unique_sources {}'.format(total_unique_sources))
        print('unique_sources_lost {}'.format(unique_sources_lost))
        print('unique_sources_found {}'.format(unique_sources_found))

        return lost_d, data_d


class PlotScampPerformanceLostOnes:

    def __init__(self):
        """

        """
        self.plot_size = [16.53, 11.69]
        self.plot_dpi = 100

        data_df = self.read_data()
        self.draw_plots(data_df)

    def read_data(self):
        """

        :return:
        """
        data_df = read_csv('data.csv', index_col=0)
        return data_df

    def draw_plots(self, data_df):
        """

        :param data_df:
        :return:
        """
        # unique_sources = list(set(data_df['IN_SOURCE']))

        selected_pms_df = data_df[data_df['IN_PM_NORM'].isin([0.1, 0.3,
                                                              1.0, 3.0])]
        unique_sources = list(set(selected_pms_df['IN_SOURCE']))

        with PdfPages('data.pdf') as pdf:
            for idx_in_source_, in_source_ in enumerate(unique_sources):
                source_df = data_df[data_df['IN_SOURCE'].isin([in_source_])]
                # MAG_ISO Galaxies
                for idx_d_, d_ in enumerate(source_df['DITHER']):
                    dither_df = source_df[source_df['DITHER'].isin([d_])]
                    if idx_d_ == 0:
                        p_alpha = float(dither_df['IN_ALPHA'].iloc[0])
                        p_delta = float(dither_df['IN_DELTA'].iloc[0])

                    input_regions = 'regions/cat_clean_ssos_ds9_{}.reg'.format(d_)
                    extracted_regions = 'regions/extracted_ds9_{}.reg'.format(d_)
                    input_fits = 'coadd_{}.fits'.format(d_)
                    alpha = float(dither_df['IN_ALPHA'].iloc[0])
                    delta = float(dither_df['IN_DELTA'].iloc[0])
                    source = float(dither_df['IN_SOURCE'].iloc[0])
                    dither = float(dither_df['DITHER'].iloc[0])
                    pm = float(dither_df['IN_PM_NORM'].iloc[0])

                    # Launch ds9 command
                    image = self.get_image(extracted_regions, input_regions,
                                           input_fits, alpha, delta, source,
                                           dither, pm, p_alpha, p_delta)
                    image = image[1:-50, :]  # Removes ds9's legend

                    fig = plt.figure(figsize=self.plot_size, dpi=self.plot_dpi)
                    ax_1 = fig.add_subplot(1, 1, 1)
                    ax_1.set_title('input source: {} - i_pm {}'.format(in_source_, pm))

                    plt.imshow(image)

                    pdf.savefig()  # saves figure
                    plt.clf()  # clear current figure
                    plt.close(fig)  # removes figure

    def get_image(self, extracted_regions, input_regions, input_fits,
                  alpha, delta, source, dither, pm, p_alpha, p_delta):
        """

        :param extracted_regions:
        :param input_regions:
        :param input_fits:
        :param alpha:
        :param delta:
        :param source:
        :param dither:
        :param pm:
        :return:
        """
        print(pm)
        d = DS9()
        # Open file
        open_file = 'file {}'.format(input_fits)
        d.set(open_file)
        # Changes scale
        scale = 'scale zscale'
        d.set(scale)
        # Loads extracted regions
        d.set('regions load {}'.format(extracted_regions))
        # Loads input regions
        d.set('regions load {}'.format(input_regions))
        # Saves image
        pan_object = 'pan to {} {} wcs fk5 degree'.format(alpha, delta)
        d.set(pan_object)
        zoom_value = get_zoom_level(pm)
        zoom_object = 'zoom to {}'.format(zoom_value)
        d.set(zoom_object)
        crosshair = 'crosshair {} {} wcs fk5 degree'.format(p_alpha, p_delta)
        d.set(crosshair)
        save_image = 'saveimage png {}_{}.png 100'.format(source, dither)
        d.set(save_image)

        img = mpimg.imread('{}_{}.png'.format(source, dither))

        return img


if __name__ == "__main__":
    ScampPerformanceLostOnes()
    # PlotScampPerformanceLostOnes()

# # !/usr/bin/python
# # -*- coding: utf-8 -*-
#
# """Python script for performance
#
# Versions:
# - 0.1 Initial release
#
# Todo:
#     *
#
# """
# from datetime import datetime, timedelta
# from decimal import Decimal
# from math import modf
#
# from astropy import units
# from astropy.coordinates import Angle
#
# from matplotlib.backends.backend_pdf import PdfPages
# from matplotlib import cm
# from matplotlib.ticker import FormatStrFormatter
# import matplotlib.image as mpimg
# import matplotlib.pyplot as plt
# import numpy as np
# from pyds9 import DS9
#
# from misc import significant_l, extract_settings
#
# __author__ = "Samuel Góngora García"
# __copyright__ = "Copyright 2017"
# __credits__ = ["Samuel Góngora García"]
# __version__ = "0.1"
# __maintainer__ = "Samuel Góngora García"
# __email__ = "sgongora@cab.inta-csic.es"
# __status__ = "Development"
#
#
# def get_regions(input_fits, sex_cf, prfs_d, mag):
#     """
#
#     :param input_fits:
#     :param sex_cf:
#     :param prfs_d:
#     :param mag:
#     :return:
#     """
#     #
#     regions_dir = '{}/{}/CCDs/{}'.format(prfs_d['fits_dir'], mag, sex_cf)
#     regions_name = '{}.reg'.format(input_fits[-27:-5])
#     regions_file = '{}/{}'.format(regions_dir, regions_name)
#
#     return regions_file
#
#
# def image(input_regions, input_fits, regions, alpha, delta, idx, pm,
#           p_alpha, p_delta, zoom):
#     """
#
#     :param input_regions:
#     :param input_fits:
#     :param regions:
#     :param alpha:
#     :param delta:
#     :param idx:
#     :param pm:
#     :param p_alpha:
#     :param p_delta:
#     :param zoom:
#     :return:
#     """
#     d = DS9()
#     # Open file
#     open_file = 'file {}'.format(input_fits)
#     d.set(open_file)
#     # Changes scale
#     scale = 'scale zscale'
#     d.set(scale)
#     # Loads regions with format
#     d.set('regions load {}'.format(regions))
#     # Loads input regions
#     d.set('regions load {}'.format(input_regions))
#     # Saves image
#     pan_object = 'pan to {} {} wcs fk5 degree'.format(alpha, delta)
#     d.set(pan_object)
#     zoom_value = get_zoom_level(pm, zoom)
#     zoom_object = 'zoom to {}'.format(zoom_value)
#
#     d.set(zoom_object)
#     if idx > 0:
#         crosshair = 'crosshair {} {} wcs fk5 degree'.format(p_alpha, p_delta)
#         d.set(crosshair)
#     else:
#         pass
#     save_image = 'saveimage png {}.png 100'.format(idx)
#     d.set(save_image)
#
#     img = mpimg.imread('{}.png'.format(idx))
#
#     # delete image
#
#     return img
#
#
# def create_labels(datetimes):
#     """
#
#     :param datetimes:
#     :return: labels
#     """
#     labels = []
#     for datetime_ in datetimes:
#         label_ = '{}.{}'.format(datetime_.second,
#                                 datetime_.microsecond)
#         labels.append(label_)
#
#     return labels
#
#
# def round_number(number):
#     """
#     TODO Improve description
#
#
#     :param number:
#     :return:
#     """
#     number_t = significant_l(number)
#
#     # Rounds the float value of difference
#     first_digit = '%.2E' % Decimal(number - number_t)
#
#     if first_digit[0] == '-':
#         first_digit = int(first_digit[1])
#     else:
#         first_digit = int(first_digit[0])
#
#     # Redondea al alza el ultimo digio
#     if first_digit > 5:
#         last_digit = str(number_t)[-1]
#         last_digit = int(last_digit)
#         last_digit += 1
#         number_t = list(str(number_t))
#         number_t[-1] = last_digit
#         number_t = [str(i) for i in number_t]
#         number_t = ''.join(number_t)
#         number_t = float(number_t)
#     else:
#         pass  # do nothing
#
#     return number_t
#
#
# def create_x_ticks(epoch_seconds):
#     """
#
#     :param epoch_seconds:
#     :return: x_ticks
#     """
#     x_difference = float(max(epoch_seconds)) - float(min(epoch_seconds))
#     x_major_stp = (x_difference / 3)
#     x_minor_stp = (x_difference / 6)
#
#     # X-SCALE
#     # Major steps
#     x_major_stps = np.arange((epoch_seconds[0] - x_major_stp * 1),
#                              (epoch_seconds[0] + x_major_stp * 4), x_major_stp)
#
#     # Minor steps
#     x_minor_stps = np.arange((epoch_seconds[0] - x_minor_stp * 2),
#                              (epoch_seconds[0] + x_minor_stp * 8),
#                              x_minor_stp)
#
#     x_ticks = {'major_t': x_major_stps, 'minor_t': x_minor_stps,
#                'major_s': x_major_stp, 'minor_s': x_minor_stp}
#
#     return x_ticks
#
#
# def create_y_ticks(delta_seconds):
#     """
#     step = stp
#
#     :param delta_seconds:
#     :return: y_ticks
#     """
#     divisions = 4  #
#
#     # Gets the major step between ticks thought the difference
#     # between the maximum and the minimum value of alpha
#     difference = float(max(delta_seconds)) - float(min(delta_seconds))
#     major_stp = (difference / divisions)
#     #
#     major_stp = float(round_number(major_stp))
#     minor_stp = (major_stp / 4)
#
#     # Gets maximum decimal position of major step
#     decimals = int(str(major_stp)[::-1].find('.'))
#
#     # Major step list starts two times before and end two times after
#     # known values
#     major_stps = np.arange(round(min(delta_seconds), decimals) - major_stp * 2,
#                            round(max(delta_seconds), decimals) + major_stp * 2,
#                            major_stp)
#     # Minor step list starts eight times before and fend eight times
#     # after know values
#     minor_stps = np.arange(round(min(delta_seconds), decimals) - minor_stp * 8,
#                            round(max(delta_seconds), decimals) + minor_stp * 8,
#                            minor_stp)
#
#     y_ticks = {'major_t': major_stps, 'minor_t': minor_stps,
#                'major_s': major_stp, 'minor_s': minor_stp}
#
#     return y_ticks
#
#
# class PlotConfidence:
#
#     def __init__(self, output_path, source, pm, mode, fitted_d, tmp_d,
#                  fits, mag):
#         """
#
#         :param output_path:
#         :param source:
#         :param pm:
#         :param mode:
#         :param fitted_d:
#         :param tmp_d: a dictionary ['alpha', 'delta',
#                                     'error_a', 'error_b', 'epoch']
#         """
#         self.prfs_d = extract_settings()
#         self.output_path = output_path
#         self.source = source
#         self.pm = pm
#         self.mode = mode  # Can be 'i' or 'o'
#         self.fitted_d = fitted_d
#         self.tmp_d = tmp_d
#         self.fits_files = fits
#         self.mag = mag
#         self.ccds = self.gets_ccd_names()
#
#         # Page size
#         self.plot_size = [16.53, 11.69]
#         self.plot_dpi = 100
#
#         self.plot()
#
#     def gets_ccd_names(self):
#         """
#
#         :return:
#         """
#         ccds = []
#         for ccd in self.fits_files:
#             ccds.append(ccd[-13:-5])
#
#         return ccds
#
#     def plot(self):
#         """
#
#         :return:
#         """
#         pdf_name = '{}/{}_{}_{}.pdf'.format(self.output_path, self.source,
#                                             self.pm, self.mode)
#         with PdfPages(pdf_name) as pdf:
#             # ALPHA PARAMETERS
#             fig = plt.figure(figsize=self.plot_size, dpi=self.plot_dpi)
#             ax = fig.add_subplot(1, 1, 1)
#             chi_squared = float("{0:.6f}".format(self.fitted_d['ra']))
#             chi_squared_title = 'chi_squared: {}'.format(chi_squared)
#
#             pm_alpha_cat = float("{0:.6f}".format(self.tmp_d['i_pm_alpha'][0]))
#             pm_alpha_cat_str = 'catalog {}'.format(pm_alpha_cat)
#
#             pm_alpha_ext = float("{0:.6f}".format(self.tmp_d['o_pm_alpha'][0]))
#             pm_alpha_ext_str = 'extracted {}'.format(pm_alpha_ext)
#
#             o_pm_alpha = float(self.tmp_d['o_pm_alpha'][0])
#             o_pm_alpha_err = float(self.tmp_d['o_pm_alpha_err'][0])
#             o_pm_alpha_sn = o_pm_alpha / o_pm_alpha_err
#             o_pm_alpha_sn = float("{0:.6f}".format(o_pm_alpha_sn))
#             pm_alpha_sn_str = 'SN {}'.format(o_pm_alpha_sn)
#             pm_alpha_ext_err = self.tmp_d['o_pm_alpha_err'][0]
#             pm_alpha_ext_err = float("{0:.6f}".format(pm_alpha_ext_err))
#             pm_alpha_ext_err_str = 'pm_error {}'.format(pm_alpha_ext_err)
#
#             pm_alpha_comp = '{} / {}'.format(pm_alpha_ext_str,
#                                              pm_alpha_cat_str)
#             pm_error_comp = '{} / {}'.format(pm_alpha_sn_str,
#                                              pm_alpha_ext_err_str)
#
#             ax.set_title('{}\n{}\n{}'.format(chi_squared_title, pm_alpha_comp,
#                                              pm_error_comp))
#
#             alpha_tmpstmp = []
#             alpha_seconds = []
#             tmp_hour = []
#             tmp_minute = []
#             for alpha_ in self.tmp_d['alpha']:
#                 a = Angle(alpha_, units.degree)
#                 dms = a.dms
#                 degree = int(dms[0])
#                 tmp_hour.append(degree)
#                 minute = int(dms[1])
#                 tmp_minute.append(minute)
#                 # second = float("{0:.6f}".format(dms[2]))
#                 second = float("{0:.6f}".format(a.arcsecond))
#                 alpha_tmpstmp.append('{}:{}:{}'.format(degree, minute, second))
#                 alpha_seconds.append('{}'.format(second))
#
#             # X-AXIS Shared between input and output catalogs
#             epoch_seconds = []
#             for epoch_ in self.tmp_d['epoch']:
#                 frac, whole = modf(epoch_)
#                 seconds_ = frac * 365.25 * 24 * 60
#                 epoch_seconds.append(seconds_)
#
#             # Creates x ticks (major and minor ones)
#             x_ticks = create_x_ticks(epoch_seconds)
#
#             # Y-AXIS
#             myformat = mdates.DateFormatter('%S.%f')
#             ax.yaxis.set_major_formatter(myformat)
#
#             # format elements in alpha_seconds list to floats
#             alpha_seconds = [float(i) for i in alpha_seconds]  # needed?
#             alpha_seconds = [float("{0:.6f}".format(i)) for i in
#                              alpha_seconds]
#
#             # Check if all hours/minutes are same
#             if len(list(set(tmp_hour))) != 1:
#                 raise Exception
#             if len(list(set(tmp_minute))) != 1:
#                 raise Exception
#
#             y_alpha_ticks = create_y_ticks(alpha_seconds)
#
#             # x-ticks assignation
#             ax.set_xticks(x_ticks['major_t'], minor=False)
#             ax.set_xticks(x_ticks['minor_t'], minor=True)
#             # x-ticks labels
#
#             # y-ticks assignation
#             ax.set_yticks(y_alpha_ticks['major_t'], minor=False)
#             ax.set_yticks(y_alpha_ticks['minor_t'], minor=True)
#             # y-ticks labels
#             empty_string_labels = [''] * len(y_alpha_ticks['major_t'])
#             ax.set_yticklabels(empty_string_labels, minor=False)
#             ax.set_yticklabels(y_alpha_ticks['minor_t'], minor=True)
#
#             # Formats grids
#             ax.grid(b=True, which='major', linestyle='-', linewidth=2)
#             ax.grid(b=True, which='minor', linestyle='--', linewidth=1)
#
#             # Annotations
#             for idx_alpha in range(0, len(alpha_tmpstmp), 1):
#                 x = epoch_seconds[idx_alpha]  # x position
#                 y = alpha_seconds[idx_alpha]  # y position
#                 # variables for position and error representation
#                 alpha_position = alpha_tmpstmp[idx_alpha]
#                 alpha_str = '   alpha {}'.format(alpha_position)
#                 error = self.tmp_d['error_b'][idx_alpha] * 3600  # error-y
#                 error_fmt = float("{0:.6f}".format(error))  # error rounded
#                 error_str = '   error {}'.format(error_fmt)
#                 # annotation
#                 ax.annotate('{}\n{}'.format(alpha_str, error_str),
#                             xy=(x, y), textcoords='data', fontsize=13)
#
#             for idx_alpha in range(0, len(alpha_tmpstmp), 1):
#                 x = epoch_seconds[idx_alpha]  # x position
#                 y = alpha_seconds[idx_alpha]  # y position
#                 error = float(self.tmp_d['error_b'][idx_alpha] * 3600)
#                 ax.errorbar(x, y, yerr=error, fmt='o',
#                             ecolor='g', capthick=2, elinewidth=4)
#
#             # Axis labels creation
#             # x-axis
#             x_label = 'EPOCH (seconds)'
#             ax.set_xlabel(x_label)
#             # y-axis
#             y_label_ra = 'Right ascension (")\n'
#             major_s = y_alpha_ticks['major_s']
#             y_label_major_step = 'major step size {}"\n'.format(major_s)
#             minor_s = y_alpha_ticks['minor_s']
#             y_label_minor_step = 'minor step size {}"'.format(minor_s)
#             ax.set_ylabel('{}{}{}'.format(y_label_ra, y_label_major_step,
#                                           y_label_minor_step))
#
#             # Plots data
#             ax.plot(epoch_seconds, alpha_seconds, 'bs', markersize=6,
#                     label='extracted position')
#             # In order to avoid scientific notation plot should be redrawn
#             ax = plt.gca()
#             ax.get_xaxis().get_major_formatter().set_useOffset(False)
#
#             # Plots a legend
#             plt.legend(loc=0, ncol=2, borderaxespad=0.)
#             plt.draw()
#
#             # Saves the current figure in pdf file
#             pdf.savefig()  # saves current figure
#             plt.close(fig)
#
#             #
#             # DELTA PARAMETERS
#             # Another figure is created in order to plot delta output.
#             # X-axis values are shared between both figures but Y-axis for
#             # declination's output is showed in seconds (floats) instead
#             # datetime objects.
#             #
#             fig = plt.figure(figsize=self.plot_size, dpi=self.plot_dpi)
#             ax = fig.add_subplot(1, 1, 1)
#             chi_squared = float("{0:.6f}".format(self.fitted_d['dec']))
#             chi_squared_title = 'chi_squared: {}'.format(chi_squared)
#
#             pm_delta_cat = float("{0:.6f}".format(self.tmp_d['i_pm_delta'][0]))
#             pm_delta_cat_str = 'catalog {}'.format(pm_delta_cat)
#
#             pm_delta_ext = float("{0:.6f}".format(self.tmp_d['o_pm_delta'][0]))
#             pm_delta_ext_str = 'extracted {}'.format(pm_delta_ext)
#
#             o_pm_delta = float(self.tmp_d['o_pm_delta'][0])
#             o_pm_delta_err = float(self.tmp_d['o_pm_delta_err'][0])
#             o_pm_delta_sn = o_pm_delta / o_pm_delta_err
#             o_pm_delta_sn = float("{0:.6f}".format(o_pm_delta_sn))
#             pm_delta_sn_str = 'SN {}'.format(o_pm_delta_sn)
#             pm_delta_ext_err = self.tmp_d['o_pm_delta_err'][0]
#             pm_delta_ext_err = float("{0:.6f}".format(pm_delta_ext_err))
#             pm_delta_ext_err_str = 'pm_error {}'.format(pm_delta_ext_err)
#
#             pm_delta_comp = '{} / {}'.format(pm_delta_ext_str,
#                                              pm_delta_cat_str)
#             pm_error_comp = '{} / {}'.format(pm_delta_sn_str,
#                                              pm_delta_ext_err_str)
#
#             ax.set_title('{}\n{}\n{}'.format(chi_squared_title, pm_delta_comp,
#                                              pm_error_comp))
#
#             delta_tmpstmp = []
#             delta_seconds = []
#             tmp_degree = []
#             tmp_minute = []
#             for delta_ in self.tmp_d['delta']:
#                 d = Angle(delta_, units.degree)
#                 dms = d.dms
#                 degree = int(dms[0])
#                 tmp_degree.append(degree)
#                 minute = int(dms[1])
#                 tmp_minute.append(minute)
#                 # second = float("{0:.6f}".format(dms[2]))
#                 second = float("{0:.6f}".format(a.arcsecond))
#                 delta_tmpstmp.append('{}:{}.{}'.format(degree, minute, second))
#                 delta_seconds.append('{}'.format(second))
#
#             ax.set_xlabel('EPOCH')
#
#             # Y-SCALE
#             # format elements in alpha_seconds list to floats
#             delta_seconds = [float(i) for i in delta_seconds]
#             delta_seconds = [float("{0:.6f}".format(i)) for i in delta_seconds]
#
#             # Check if all hours/minutes are same
#             if len(list(set(tmp_degree))) != 1:
#                 raise Exception
#             if len(list(set(tmp_minute))) != 1:
#                 raise Exception
#
#             y_delta_ticks = create_y_ticks(delta_seconds)
#
#             # x-ticks assignation
#             ax.set_xticks(x_ticks['major_t'], minor=False)
#             ax.set_xticks(x_ticks['minor_t'], minor=True)
#             # x-ticks labels
#
#             # y-ticks assignation
#             ax.set_yticks(y_delta_ticks['major_t'], minor=False)
#             ax.set_yticks(y_delta_ticks['minor_t'], minor=True)
#             # y-ticks labels
#             empty_string_labels = [''] * len(y_delta_ticks['major_t'])
#             ax.set_yticklabels(empty_string_labels, minor=False)
#             # Converts floats into strings
#             # y_minor_ticks = [str(i) for i in y_minor_ticks]
#             ax.set_yticklabels(y_delta_ticks['minor_t'], minor=True)
#
#             # Formats grids
#             ax.grid(b=True, which='major', linestyle='-', linewidth=2)
#             ax.grid(b=True, which='minor', linestyle='--', linewidth=1)
#
#             for idx_delta in range(0, len(delta_tmpstmp), 1):
#                 x = epoch_seconds[idx_delta]  # x position
#                 y = delta_seconds[idx_delta]  # y position
#                 # variables for position and error representation
#                 delta_position = delta_tmpstmp[idx_delta]
#                 delta_str = '   delta {}'.format(delta_position)
#                 error = self.tmp_d['error_b'][idx_delta] * 3600  # error-y
#                 error_fmt = float("{0:.6f}".format(error))  # error rounded
#                 error_str = '   error {}'.format(error_fmt)
#                 # annotation
#                 ax.annotate('{}\n{}'.format(delta_str, error_str),
#                             xy=(x, y), textcoords='data', fontsize=13)
#
#             for idx_delta in range(0, len(delta_tmpstmp), 1):
#                 x = epoch_seconds[idx_delta]  # x position
#                 y = delta_seconds[idx_delta]  # y position
#                 error = float(self.tmp_d['error_b'][idx_delta] * 3600)
#                 ax.errorbar(x, y, yerr=error, fmt='o',
#                             ecolor='g', capthick=2, elinewidth=4)
#
#             # Label creation
#             y_label_ra = 'Declination (")\n'
#             major_s = y_delta_ticks['major_s']
#             y_label_major_step = 'major step size {}"\n'.format(major_s)
#             minor_s = y_delta_ticks['minor_s']
#             y_label_minor_step = 'minor step size {}"'.format(minor_s)
#
#             ax.set_ylabel('{}{}{}'.format(y_label_ra, y_label_major_step,
#                                           y_label_minor_step))
#
#             ax.plot(epoch_seconds, delta_seconds, 'bs', markersize=6,
#                     label='extracted position')
#
#             ax = plt.gca()
#             ax.get_xaxis().get_major_formatter().set_useOffset(False)
#             plt.legend(loc=0, ncol=2, borderaxespad=0.)
#             plt.draw()
#
#             pdf.savefig()  # saves current figure
#             plt.close(fig)
#
#             #
#             # IMAGES
#             #
#
#             # Gets limits
#             alpha_center = 0
#             delta_center = 0
#             if len(self.tmp_d['alpha']) == 3:
#                 alpha_center = self.tmp_d['alpha'][1]
#                 delta_center = self.tmp_d['delta'][1]
#             elif len(self.tmp_d['alpha']) == 4:
#                 alpha_sum = self.tmp_d['alpha'][1] + self.tmp_d['alpha'][2]
#                 alpha_center = alpha_sum / 2
#                 delta_sum = self.tmp_d['delta'][1] + self.tmp_d['delta'][2]
#                 delta_center = delta_sum / 2
#
#             if alpha_center is 0 or delta_center is 0:
#                 raise Exception
#
#             for idx in range(0, len(self.tmp_d['alpha']), 1):
#                 # Get regions for desired fits file
#                 regions_file = get_regions(self.fits_files[idx],
#                                            self.tmp_d['sex_cf'][0],
#                                            self.prfs_d, self.mag)
#
#                 if idx == 0:
#                     p_alpha = 0
#                     p_delta = 0
#                 else:
#                     p_alpha = self.tmp_d['alpha'][idx - 1]
#                     p_delta = self.tmp_d['delta']
#                 # Creates image
#                 img = image(self.fits_files[idx], regions_file,
#                             alpha_center, delta_center, idx, p_alpha, p_delta)
#
#                 img = img[1:-50, :]
#
#                 fig = plt.figure(figsize=self.plot_size, dpi=self.plot_dpi)
#                 ax = fig.add_subplot(1, 1, 1)
#                 plt.imshow(img)
#
#                 labels = [item.get_text() for item in ax.get_xticklabels()]
#
#                 empty_string_labels = [''] * len(labels)
#                 ax.set_xticklabels(empty_string_labels)
#
#                 labels = [item.get_text() for item in ax.get_yticklabels()]
#
#                 empty_string_labels = [''] * len(labels)
#                 ax.set_yticklabels(empty_string_labels)
#
#                 ax.set_ylabel('right ascension')
#                 ax.set_xlabel('declination')
#
#                 pdf.savefig()
#                 plt.close(fig)
#
#             # Loads images in pdf file
#
#         return True
#
#
# class PlotBothConfidence:
#
#     def __init__(self, output_path, source, pm, mode, fitted_d, tmp_d):
#         """
#
#         :param output_path:
#         :param source:
#         :param pm:
#         :param mode:
#         :param fitted_d:
#         :param tmp_d: a dictionary ['i_alpha', 'i_delta', 'o_alpha', 'o_delta',
#                                     'error_a', 'error_b', 'epoch']
#         """
#         self.output_path = output_path
#         self.source = source
#         self.pm = pm
#         self.mode = mode  # Can be 'i' or 'o'
#         self.fitted_d = fitted_d
#         self.tmp_d = tmp_d
#
#         plot_size = [16.53, 11.69]
#         plot_dpi = 100
#
#         self.plot(plot_size, plot_dpi)
#
#     def plot(self, plot_size, plot_dpi):
#         """
#
#         :param plot_size:
#         :param plot_dpi:
#         :return:
#         """
#         pdf_name = '{}/{}_{}_{}.pdf'.format(self.output_path, self.source,
#                                             self.pm, self.mode)
#         with PdfPages(pdf_name) as pdf:
#             # ALPHA PARAMETERS
#             fig = plt.figure(figsize=plot_size, dpi=plot_dpi)
#             ax = fig.add_subplot(1, 1, 1)
#             chi_squared = float("{0:.6f}".format(self.fitted_d['ra']))
#             chi_squared_title = 'chi_squared: {}'.format(chi_squared)
#
#             pm_alpha_cat = float("{0:.6f}".format(self.tmp_d['i_pm_alpha'][0]))
#             pm_alpha_cat_str = 'catalog {}'.format(pm_alpha_cat)
#
#             pm_alpha_ext = float("{0:.6f}".format(self.tmp_d['o_pm_alpha'][0]))
#             pm_alpha_ext_str = 'extracted {}'.format(pm_alpha_ext)
#
#             o_pm_alpha = float(self.tmp_d['o_pm_alpha'][0])
#             o_pm_alpha_err = float(self.tmp_d['o_pm_alpha_err'][0])
#             o_pm_alpha_sn = o_pm_alpha / o_pm_alpha_err
#             o_pm_alpha_sn = float("{0:.6f}".format(o_pm_alpha_sn))
#             pm_alpha_sn_str = 'SN {}'.format(o_pm_alpha_sn)
#             pm_alpha_ext_err = self.tmp_d['o_pm_alpha_err'][0]
#             pm_alpha_ext_err = float("{0:.6f}".format(pm_alpha_ext_err))
#             pm_alpha_ext_err_str = 'pm_error {}'.format(pm_alpha_ext_err)
#
#             pm_alpha_comp = '{} / {}'.format(pm_alpha_ext_str,
#                                              pm_alpha_cat_str)
#             pm_error_comp = '{} / {}'.format(pm_alpha_sn_str,
#                                              pm_alpha_ext_err_str)
#
#             ax.set_title('{}\n{}\n{}'.format(chi_squared_title, pm_alpha_comp,
#                                              pm_error_comp))
#
#             # input_parameters
#             plot_flag = True
#             i_a_tmpstmp = []
#             i_alpha_seconds = []
#             i_tmp_hour = []
#             i_tmp_minute = []
#             for alpha_ in self.tmp_d['i_alpha']:  # fixme change to dms
#                 a = Angle(alpha_, units.degree)
#                 hms = a.hms
#                 hour = int(hms[0])
#                 i_tmp_hour.append(hour)
#                 minute = int(hms[1])
#                 i_tmp_minute.append(minute)
#                 second = float("{0:.6f}".format(hms[2]))
#                 i_a_tmpstmp.append('{}:{}:{}'.format(hour, minute, second))
#                 i_alpha_seconds.append('{}'.format(second))
#
#             i_datetimes = []
#             for idx, i_datetime_ in enumerate(i_a_tmpstmp):
#                 i_datetime_tmp = datetime.strptime(i_datetime_, "%H:%M:%S.%f")
#                 i_datetimes.append(i_datetime_tmp)
#
#             # output parameters
#             o_a_tmpstmp = []
#             o_alpha_seconds = []
#             o_tmp_hour = []
#             o_tmp_minute = []
#             for alpha_ in self.tmp_d['o_alpha']:
#                 a = Angle(alpha_, units.degree)
#                 hms = a.hms
#                 hour = int(hms[0])
#                 o_tmp_hour.append(hour)
#                 minute = int(hms[1])
#                 o_tmp_minute.append(minute)
#                 second = float("{0:.6f}".format(hms[2]))
#                 o_a_tmpstmp.append('{}:{}:{}'.format(hour, minute, second))
#                 o_alpha_seconds.append('{}'.format(second))
#
#             o_datetimes = []
#             for idx, o_datetime_ in enumerate(o_a_tmpstmp):
#                 o_datetime_tmp = datetime.strptime(o_datetime_, "%H:%M:%S.%f")
#                 o_datetimes.append(o_datetime_tmp)
#
#             # X-AXIS Shared between input and output catalogs
#             epoch_seconds = []
#             for epoch_ in self.tmp_d['epoch']:
#                 frac, whole = modf(epoch_)
#                 seconds_ = frac * 365.25 * 24 * 60
#                 epoch_seconds.append(seconds_)
#
#             # Creates x ticks (major and minor ones)
#             x_ticks = create_x_ticks(epoch_seconds)
#
#             # todo reformat!
#             # Y-AXIS
#             myformat = mdates.DateFormatter('%S.%f')
#             ax.yaxis.set_major_formatter(myformat)
#
#             # format elements in alpha_seconds list to floats
#             i_alpha_seconds = [float(i) for i in i_alpha_seconds]  # needed?
#             i_alpha_seconds = [float("{0:.6f}".format(i)) for i in
#                                i_alpha_seconds]
#             o_alpha_seconds = [float(i) for i in o_alpha_seconds]  # needed?
#             o_alpha_seconds = [float("{0:.6f}".format(i)) for i in
#                                o_alpha_seconds]
#
#             # Check if all hours/minutes are same
#             if len(list(set(i_tmp_hour))) != 1:
#                 plot_flag = False
#             if len(list(set(o_tmp_hour))) != 1:
#                 plot_flag = False
#             if len(list(set(i_tmp_minute))) != 1:
#                 plot_flag = False
#             if len(list(set(o_tmp_minute))) != 1:
#                 plot_flag = False
#
#             if plot_flag:
#                 alpha_seconds = i_alpha_seconds + o_alpha_seconds
#
#                 # Creates y ticks (major and minor ones)
#                 y_ticks = create_y_ticks(alpha_seconds)
#
#                 y_labels = create_labels(y_ticks['minor_t'])
#
#                 # x-ticks assignation
#                 ax.set_xticks(x_ticks['major_t'], minor=False)
#                 ax.set_xticks(x_ticks['minor_t'], minor=True)
#                 # x-ticks labels
#
#                 # y-ticks assignation
#                 ax.set_yticks(y_ticks['major_t'], minor=False)
#                 ax.set_yticks(y_ticks['minor_t'], minor=True)
#                 # y-ticks labels
#                 empty_string_labels = [''] * len(y_ticks['major_t'])
#                 ax.set_yticklabels(empty_string_labels, minor=False)
#                 ax.set_yticklabels(y_labels, minor=True)
#
#                 # Formats grids
#                 ax.grid(b=True, which='major', linestyle='-', linewidth=2)
#                 ax.grid(b=True, which='minor', linestyle='--', linewidth=1)
#
#                 # Annotations
#                 for idx_datetime_, datetime_ in enumerate(i_datetimes):
#                     # Format alpha
#                     hour = datetime_.hour
#                     minute = datetime_.minute
#                     second = datetime_.second
#                     msecond = datetime_.microsecond
#                     alpha_str = '   alpha {}:{}:{}:{}'.format(hour, minute,
#                                                               second, msecond)
#
#                     # Annotate position and error associated
#                     ax.annotate('{}'.format(alpha_str),
#                                 xy=(epoch_seconds[idx_datetime_], datetime_),
#                                 textcoords='data', fontsize=13)
#
#                 # Annotations
#                 for idx_datetime_, datetime_ in enumerate(o_datetimes):
#                     # Format alpha
#                     hour = datetime_.hour
#                     minute = datetime_.minute
#                     second = datetime_.second
#                     msecond = datetime_.microsecond
#                     alpha_str = '   alpha {}:{}:{}:{}'.format(hour, minute,
#                                                               second, msecond)
#                     # Format error
#                     error_seconds = self.tmp_d['error_a'][idx_datetime_] * 3600
#                     error = float("{0:.6f}".format(error_seconds))
#                     error_str = '   error  {}"'.format(error)
#                     # Annotate position and error associated
#                     ax.annotate('{}\n{}'.format(alpha_str, error_str),
#                                 xy=(epoch_seconds[idx_datetime_], datetime_),
#                                 textcoords='data', fontsize=13)
#                 for idx_datetime_, datetime_ in enumerate(o_datetimes):
#                     error_seconds = self.tmp_d['error_a'][idx_datetime_] * 3600
#                     errorbar_size = timedelta(0, error_seconds)
#                     ax.errorbar(epoch_seconds[idx_datetime_], datetime_,
#                                 yerr=errorbar_size,
#                                 ecolor='g', capthick=2, elinewidth=4)
#
#                 # Axis labels creation
#                 # x-axis
#                 x_label = 'EPOCH (seconds)'
#                 ax.set_xlabel(x_label)
#                 # y-axis
#                 y_label_ra = 'Right ascension (")\n'
#                 major_s = y_ticks['major_s']
#                 y_label_major_step = 'major step size {}"\n'.format(major_s)
#                 minor_s = y_ticks['minor_s']
#                 y_label_minor_step = 'minor step size {}"'.format(minor_s)
#                 ax.set_ylabel('{}{}{}'.format(y_label_ra, y_label_major_step,
#                                               y_label_minor_step))
#
#                 # Plots data
#                 ax.plot(epoch_seconds, i_datetimes, 'bs', markersize=6,
#                         label='catalog position')
#                 ax.plot(epoch_seconds, o_datetimes, 'rs', markersize=6,
#                         label='extracted position')
#                 ax.plot(epoch_seconds, i_datetimes, linewidth=1)
#                 # In order to avoid scientific notation plot should be redrawn
#                 ax = plt.gca()
#                 ax.get_xaxis().get_major_formatter().set_useOffset(False)
#
#                 # Plots a legend
#                 plt.legend(loc=0, ncol=2, borderaxespad=0.)
#                 plt.draw()
#
#                 # Saves the current figure in pdf file
#                 pdf.savefig()  # saves current figure
#             plt.clf()  # clear current figure
#             plt.close(fig)
#
#             #
#             # DELTA PARAMETERS
#             # Another figure is created in order to plot delta output.
#             # X-axis values are shared between both figures but Y-axis for
#             # declination's output is showed in seconds (floats) instead
#             # datetime objects.
#             #
#             fig = plt.figure(figsize=plot_size, dpi=plot_dpi)
#             ax = fig.add_subplot(1, 1, 1)
#             chi_squared = float("{0:.6f}".format(self.fitted_d['dec']))
#             chi_squared_title = 'chi_squared: {}'.format(chi_squared)
#
#             pm_delta_cat = float("{0:.6f}".format(self.tmp_d['i_pm_delta'][0]))
#             pm_delta_cat_str = 'catalog {}'.format(pm_delta_cat)
#
#             pm_delta_ext = float("{0:.6f}".format(self.tmp_d['o_pm_delta'][0]))
#             pm_delta_ext_str = 'extracted {}'.format(pm_delta_ext)
#
#             o_pm_delta = float(self.tmp_d['o_pm_delta'][0])
#             o_pm_delta_err = float(self.tmp_d['o_pm_delta_err'][0])
#             o_pm_delta_sn = o_pm_delta / o_pm_delta_err
#             o_pm_delta_sn = float("{0:.6f}".format(o_pm_delta_sn))
#             pm_delta_sn_str = 'SN {}'.format(o_pm_delta_sn)
#             pm_delta_ext_err = self.tmp_d['o_pm_delta_err'][0]
#             pm_delta_ext_err = float("{0:.6f}".format(pm_delta_ext_err))
#             pm_delta_ext_err_str = 'pm_error {}'.format(pm_delta_ext_err)
#
#             pm_delta_comp = '{} / {}'.format(pm_delta_ext_str,
#                                              pm_delta_cat_str)
#             pm_error_comp = '{} / {}'.format(pm_delta_sn_str,
#                                              pm_delta_ext_err_str)
#
#             ax.set_title('{}\n{}\n{}'.format(chi_squared_title, pm_delta_comp,
#                                              pm_error_comp))
#
#             plot_flag = True
#             i_d_tmpstmp = []
#             i_delta_seconds = []
#             i_tmp_degree = []
#             i_tmp_minute = []
#             for delta_ in self.tmp_d['i_delta']:
#                 d = Angle(delta_, units.degree)
#                 dms = d.dms
#                 degree = int(dms[0])
#                 i_tmp_degree.append(degree)
#                 minute = int(dms[1])
#                 i_tmp_minute.append(minute)
#                 second = float("{0:.6f}".format(dms[2]))
#                 i_d_tmpstmp.append('{}:{}.{}'.format(degree, minute, second))
#                 i_delta_seconds.append('{}'.format(second))
#
#             o_d_tmpstmp = []
#             o_delta_seconds = []
#             o_tmp_degree = []
#             o_tmp_minute = []
#             for delta_ in self.tmp_d['o_delta']:
#                 d = Angle(delta_, units.degree)
#                 dms = d.dms
#                 degree = int(dms[0])
#                 o_tmp_degree.append(degree)
#                 minute = int(dms[1])
#                 o_tmp_minute.append(minute)
#                 second = float("{0:.6f}".format(dms[2]))
#                 o_d_tmpstmp.append('{}:{}.{}'.format(degree, minute, second))
#                 o_delta_seconds.append('{}'.format(second))
#
#             ax.set_xlabel('EPOCH')
#
#             # Y-SCALE
#             # format elements in alpha_seconds list to floats
#             i_delta_seconds_ = []
#             for idx_, i_delta_second_ in enumerate(i_delta_seconds):
#                 i_delta_tmp = float("{0:.6f}".format(float(i_delta_second_)))
#                 i_delta_seconds_.append(i_delta_tmp)
#             i_delta_seconds = i_delta_seconds_  # it's really needed?
#
#             o_delta_seconds_ = []
#             for idx_, o_delta_second_ in enumerate(o_delta_seconds):
#                 o_delta_tmp = float("{0:.6f}".format(float(o_delta_second_)))
#                 o_delta_seconds_.append(o_delta_tmp)
#             o_delta_seconds = o_delta_seconds_  # it's really needed?
#
#             delta_seconds = i_delta_seconds + o_delta_seconds
#
#             # Check if all hours/minutes are same
#             if len(list(set(i_tmp_degree))) != 1:
#                 plot_flag = False
#             if len(list(set(i_tmp_minute))) != 1:
#                 plot_flag = False
#             if len(list(set(o_tmp_degree))) != 1:
#                 plot_flag = False
#             if len(list(set(o_tmp_minute))) != 1:
#                 plot_flag = False
#
#             if plot_flag:
#                 y_ticks = create_y_ticks(delta_seconds)
#
#                 # x-ticks assignation
#                 ax.set_xticks(x_ticks['major_t'], minor=False)
#                 ax.set_xticks(x_ticks['minor_t'], minor=True)
#                 # x-ticks labels
#
#                 # y-ticks assignation
#                 ax.set_yticks(y_ticks['major_t'], minor=False)
#                 ax.set_yticks(y_ticks['minor_t'], minor=True)
#                 # y-ticks labels
#                 empty_string_labels = [''] * len(y_ticks['major_t'])
#                 ax.set_yticklabels(empty_string_labels, minor=False)
#                 # Converts floats into strings
#                 # y_minor_ticks = [str(i) for i in y_minor_ticks]
#                 ax.set_yticklabels(y_ticks['minor_t'], minor=True)
#
#                 # Formats grids
#                 ax.grid(b=True, which='major', linestyle='-', linewidth=2)
#                 ax.grid(b=True, which='minor', linestyle='--', linewidth=1)
#
#                 for idx_datetime_, datetime_ in enumerate(i_datetimes):
#                     x = epoch_seconds[idx_datetime_]  # x position
#                     y = i_delta_seconds[idx_datetime_]  # y position
#                     # variables for position and error representation
#                     delta_position = i_d_tmpstmp[idx_datetime_]
#                     delta_str = '   delta {}'.format(delta_position)
#                     # annotation
#                     ax.annotate('{}'.format(delta_str),
#                                 xy=(x, y), textcoords='data', fontsize=13)
#
#                 for idx_datetime_, datetime_ in enumerate(o_datetimes):
#                     x = epoch_seconds[idx_datetime_]  # x position
#                     y = o_delta_seconds[idx_datetime_]  # y position
#                     # variables for position and error representation
#                     delta_position = o_d_tmpstmp[idx_datetime_]
#                     delta_str = '   delta {}'.format(delta_position)
#                     error = self.tmp_d['error_b'][idx_datetime_] * 3600
#                     error_fmt = float("{0:.6f}".format(error))  # error rounded
#                     error_str = '   error {}'.format(error_fmt)
#                     # annotation
#                     ax.annotate('{}\n{}'.format(delta_str, error_str),
#                                 xy=(x, y), textcoords='data', fontsize=13)
#
#                 for idx_datetime_, datetime_ in enumerate(o_datetimes):
#                     x = epoch_seconds[idx_datetime_]  # x position
#                     y = o_delta_seconds[idx_datetime_]  # y position
#                     error = float(self.tmp_d['error_b'][idx_datetime_] * 3600)
#                     ax.errorbar(x, y, yerr=error,
#                                 ecolor='g', capthick=2, elinewidth=4)
#
#                 # Label creation
#                 y_label_ra = 'Declination (")\n'
#                 major_s = y_ticks['major_s']
#                 y_label_major_step = 'major step size {}"\n'.format(major_s)
#                 minor_s = y_ticks['minor_s']
#                 y_label_minor_step = 'minor step size {}"'.format(minor_s)
#                 ax.set_ylabel('{}{}{}'.format(y_label_ra, y_label_major_step,
#                                               y_label_minor_step))
#
#                 # Plots data
#                 ax.plot(epoch_seconds, i_delta_seconds, 'bs', markersize=6,
#                         label='extracted position')
#                 ax.plot(epoch_seconds, o_delta_seconds, 'rs', markersize=6,
#                         label='extracted position')
#                 ax.plot(epoch_seconds, i_delta_seconds, linewidth=1)
#                 # In order to avoid scientific notation plot should be redrawn
#                 ax = plt.gca()
#                 ax.get_xaxis().get_major_formatter().set_useOffset(False)
#                 plt.legend(loc=0, ncol=2, borderaxespad=0.)
#                 plt.draw()
#
#                 pdf.savefig()  # saves current figure
#
#             plt.close(fig)
#
#
# class PlotError:
#
#     def __init__(self, output_path, err_d, fits_, mag, ok, fitted_d):
#         """
#
#         :param output_path:
#         :param err_d:
#         """
#         self.prfs_d = extract_settings()
#         self.output_path = output_path
#         self.err_d = err_d
#         self.fits_files = fits_
#         self.ccds = self.gets_ccd_names(fits_)
#         self.mag = mag
#         self.ok = ok
#         self.fitted_d = fitted_d
#
#         # Page size
#         self.plot_size = [16.53, 11.69]
#         self.plot_dpi = 100
#
#         self.plot()
#
#     def gets_ccd_names(self, fits_):
#         """
#
#         :return:
#         """
#         ccds = []
#         for ccd in fits_:
#             ccds.append(ccd[-13:-5])
#
#         return ccds
#
#     def plot(self):
#         """
#
#         :return:
#         """
#         pdf_name = '{}/{}.pdf'.format(self.output_path, self.err_d['source'][0])
#
#         with PdfPages(pdf_name) as pdf:
#             # ALPHA PARAMETERS
#             fig = plt.figure(figsize=self.plot_size, dpi=self.plot_dpi)
#             ax_1 = plt.subplot2grid((1, 5), (0, 0), colspan=4)
#
#             #
#             ax_1.xaxis.set_major_formatter(FormatStrFormatter('%.6f'))
#             ax_1.yaxis.set_major_formatter(FormatStrFormatter('%.6f'))
#             #
#
#             source_str = 'source {}'.format(self.err_d['source'][0])
#             pm_str = 'input_pm {}'.format(self.err_d['i_pm'][0])
#             # ok_str = 'ok {}'.format(self.ok)
#             ax_1.set_title('{}\n{}'.format(source_str, pm_str))
#
#             # alpha coordinates
#             i_alpha_tmpstmp = []
#             i_alpha_seconds = []
#             i_alpha_arcseconds = []
#             tmp_hour = []
#             tmp_minute = []
#             for alpha_ in self.err_d['i_alpha']:
#                 a = Angle(alpha_, units.degree)
#                 dms = a.dms
#                 degree = int(dms[0])
#                 tmp_hour.append(degree)
#                 minute = int(dms[1])
#                 tmp_minute.append(minute)
#                 second = float("{0:.6f}".format(dms[2]))
#                 arcsecond = float("{0:.6f}".format(a.arcsecond))
#                 i_alpha_tmpstmp.append('{}:{}:{}'.format(degree, minute,
#                                                          second))
#                 i_alpha_seconds.append('{}'.format(second))
#                 i_alpha_arcseconds.append(arcsecond)
#
#             o_alpha_tmpstmp = []
#             o_alpha_seconds = []
#             tmp_hour = []
#             tmp_minute = []
#             for alpha_ in self.err_d['o_alpha']:
#                 if alpha_ is not False:
#                     a = Angle(alpha_, units.degree)
#                     dms = a.dms
#                     degree = int(dms[0])
#                     tmp_hour.append(degree)
#                     minute = int(dms[1])
#                     tmp_minute.append(minute)
#                     second = float("{0:.6f}".format(dms[2]))
#                     o_alpha_tmpstmp.append('{}:{}:{}'.format(degree, minute,
#                                                              second))
#                     o_alpha_seconds.append('{}'.format(second))
#
#             alpha_seconds = i_alpha_seconds + o_alpha_seconds
#
#             # delta coordinates
#             i_delta_tmpstmp = []
#             i_delta_seconds = []
#             tmp_degree = []
#             tmp_minute = []
#             for delta_ in self.err_d['i_delta']:
#                 d = Angle(delta_, units.degree)
#                 dms = d.dms
#                 degree = int(dms[0])
#                 tmp_degree.append(degree)
#                 minute = int(dms[1])
#                 tmp_minute.append(minute)
#                 second = float("{0:.6f}".format(dms[2]))
#                 i_delta_tmpstmp.append(
#                     '{}:{}.{}'.format(degree, minute, second))
#                 i_delta_seconds.append('{}'.format(second))
#
#             o_delta_tmpstmp = []
#             o_delta_seconds = []
#             tmp_degree = []
#             tmp_minute = []
#             for delta_ in self.err_d['o_delta']:
#                 if delta_ is not False:
#                     d = Angle(delta_, units.degree)
#                     dms = d.dms
#                     degree = int(dms[0])
#                     tmp_degree.append(degree)
#                     minute = int(dms[1])
#                     tmp_minute.append(minute)
#                     second = float("{0:.6f}".format(dms[2]))
#                     o_delta_tmpstmp.append('{}:{}.{}'.format(degree, minute,
#                                                              second))
#                     o_delta_seconds.append('{}'.format(second))
#
#             delta_seconds = i_delta_seconds + o_delta_seconds
#
#             i_alpha_seconds = [float(i) for i in i_alpha_seconds]  # needed?
#             i_alpha_seconds = [float("{0:.6f}".format(i)) for i in
#                                i_alpha_seconds]
#             i_delta_seconds = [float(i) for i in i_delta_seconds]  # needed?
#             i_delta_seconds = [float("{0:.6f}".format(i)) for i in
#                                i_delta_seconds]
#
#             o_alpha_seconds = [float(i) for i in o_alpha_seconds]  # needed?
#             o_alpha_seconds = [float("{0:.6f}".format(i)) for i in
#                                o_alpha_seconds]
#             o_delta_seconds = [float(i) for i in o_delta_seconds]  # needed?
#             o_delta_seconds = [float("{0:.6f}".format(i)) for i in
#                                o_delta_seconds]
#
#             alpha_seconds = [float(i) for i in alpha_seconds]  # needed?
#             alpha_seconds = [float("{0:.6f}".format(i)) for i in alpha_seconds]
#             delta_seconds = [float(i) for i in delta_seconds]  # needed?
#             delta_seconds = [float("{0:.6f}".format(i)) for i in delta_seconds]
#
#             # Plots data
#             blue_colors = cm.Blues(np.linspace(0, 1, len(i_alpha_seconds) + 1))
#             for idx in range(0, len(i_alpha_seconds), 1):
#                 ax_1.scatter(i_alpha_seconds[idx],
#                              i_delta_seconds[idx],
#                              c=blue_colors[idx + 1], s=32)
#             """
#             ax_1.plot(i_alpha_seconds, i_delta_seconds, 'bs', markersize=6,
#                       label='catalog position')
#             """
#             # print('o_alpha {}'.format(o_alpha_seconds))
#             # print('o_delta {}'.format(o_delta_seconds))
#             # print('i_alpha {}'.format(i_alpha_seconds))
#             # print('i_delta {}'.format(i_delta_seconds))
#             """
#             ax_1.plot(o_alpha_seconds, o_delta_seconds, 'rs', markersize=6,
#                       label='extracted position')
#             """
#             red_colors = cm.Reds(np.linspace(0, 1, len(o_alpha_seconds) + 1))
#             for idx in range(0, len(o_alpha_seconds), 1):
#                 ax_1.scatter(o_alpha_seconds[idx],
#                              o_delta_seconds[idx],
#                              c=red_colors[idx + 1], s=32)
#
#             try:
#                 x_ticks = create_y_ticks(alpha_seconds)
#                 y_ticks = create_y_ticks(delta_seconds)
#             except ZeroDivisionError:
#                 print('alpha: {}'.format(self.err_d['i_alpha']))
#                 print('delta: {}'.format(self.err_d['i_delta']))
#                 print('source {}'.format(self.err_d['source'][0]))
#                 raise Exception
#
#             # x-ticks assignation
#             ax_1.set_xticks(x_ticks['major_t'], minor=False)
#             ax_1.set_xticklabels(x_ticks['major_t'])
#             ax_1.set_xticks(x_ticks['minor_t'], minor=True)
#
#             # y-ticks assignation
#             ax_1.set_yticks(y_ticks['major_t'], minor=False)
#             ax_1.set_yticklabels(y_ticks['major_t'])
#             ax_1.set_yticks(y_ticks['minor_t'], minor=True)
#
#             # Formats grids
#             ax_1.grid(b=True, which='major', linestyle='-', linewidth=2)
#             ax_1.grid(b=True, which='minor', linestyle='--', linewidth=1)
#
#             # x-axis
#             x_label_ra = 'Right ascension \n'
#             major_s = x_ticks['major_s']
#             x_label_major_step = 'major step size {}"\n'.format(major_s)
#             minor_s = x_ticks['minor_s']
#             x_label_minor_step = 'minor step size {}"'.format(minor_s)
#             ax_1.set_xlabel('{}{}{}'.format(x_label_ra, x_label_major_step,
#                                             x_label_minor_step))
#             # x-axis
#             y_label_ra = 'Declination \n'
#             major_s = y_ticks['major_s']
#             y_label_major_step = 'major step size {}"\n'.format(major_s)
#             minor_s = y_ticks['minor_s']
#             y_label_minor_step = 'minor step size {}"'.format(minor_s)
#             ax_1.set_ylabel('{}{}{}'.format(y_label_ra, y_label_major_step,
#                                             y_label_minor_step))
#
#             # Plots a legend
#             plt.legend(loc=0, ncol=2, borderaxespad=0.)
#             plt.setp(ax_1.get_xticklabels(), visible=True)
#             plt.setp(ax_1.get_yticklabels(), visible=True)
#             plt.draw()
#
#             ax_2 = plt.subplot2grid((1, 5), (0, 4))
#
#             i_pm = float("{0:.6f}".format(self.err_d['i_pm'][1]))
#             i_pm_alpha = float("{0:.6f}".format(self.err_d['i_pm_alpha'][1]))
#             i_pm_delta = float("{0:.6f}".format(self.err_d['i_pm_delta'][1]))
#             o_pm = float("{0:.6f}".format(self.err_d['o_pm'][1]))
#             o_pm_alpha = float("{0:.6f}".format(self.err_d['o_pm_alpha'][1]))
#             o_pm_delta = float("{0:.6f}".format(self.err_d['o_pm_delta'][1]))
#             if type(self.fitted_d['ra']) is str:
#                 chi_ra = 'empty'
#             else:
#                 chi_ra = float("{0:.6f}".format(self.fitted_d['ra']))
#             if type(self.fitted_d['dec']) is str:
#                 chi_dec = 'empty'
#             else:
#                 chi_dec = float("{0:.6f}".format(self.fitted_d['dec']))
#
#             table_ = [['', 'catalog values'],
#                       ['cat_pm', i_pm],
#                       ['cat_pm_alpha', i_pm_alpha],
#                       ['cat_pm_delta', i_pm_delta],
#                       ['', 'extracted values'],
#                       ['ext_pm', o_pm],
#                       ['ext_pm_alpha', o_pm_alpha],
#                       ['ext_pm_delta', o_pm_delta],
#                       ['', 'chi squared'],
#                       ['chi_squared_ra', chi_ra],
#                       ['chi_squared_dec', chi_dec]]
#
#             # ax_2.axis('tight')
#             ax_2.axis('off')
#             # ax_2.axis('on')
#             data_table = ax_2.table(cellText=table_, colLabels=None,
#                                     loc='center')
#             data_table.set_fontsize(14)
#
#             # fig.tight_layout()
#
#             # Saves the current figure in pdf file
#             pdf.savefig()  # saves current figure
#             plt.close(fig)
#
#             #
#             # IMAGES
#             #
#
#             # Gets limits
#             alpha_center = 0
#             delta_center = 0
#             if len(self.err_d['i_alpha']) == 2:
#                 alpha_sum = self.err_d['i_alpha'][0] + self.err_d['i_alpha'][1]
#                 alpha_center = alpha_sum / 2
#                 delta_sum = self.err_d['i_delta'][0] + self.err_d['i_delta'][1]
#                 delta_center = delta_sum / 2
#             elif len(self.err_d['i_alpha']) == 3:
#                 alpha_center = self.err_d['i_alpha'][1]
#                 delta_center = self.err_d['i_delta'][1]
#             elif len(self.err_d['i_alpha']) == 4:
#                 alpha_sum = self.err_d['i_alpha'][1] + self.err_d['i_alpha'][2]
#                 alpha_center = alpha_sum / 2
#                 delta_sum = self.err_d['i_delta'][1] + self.err_d['i_delta'][2]
#                 delta_center = delta_sum / 2
#
#             if alpha_center is 0 or delta_center is 0:
#                 print('ERROR')
#                 print("self.err_d['i_alpha'] {}".format(self.err_d['i_alpha']))
#                 print("self.err_d['i_delta'] {}".format(self.err_d['i_delta']))
#
#                 raise Exception
#
#             for idx in range(0, len(self.err_d['i_alpha']), 1):
#                 # Get regions for desired fits file
#                 regions_file = get_regions(self.fits_files[idx],
#                                            self.err_d['sex_cf'][0],
#                                            self.prfs_d, self.mag)
#
#                 # Creates image
#                 dither = self.fits_files[idx][-6:-5]
#                 i_regs = '{}/dither_{}.reg'.format(self.prfs_d['dithers_out'],
#                                                    dither)
#                 if idx == 0:
#                     p_alpha = 0
#                     p_delta = 0
#                 else:
#                     p_alpha = self.err_d['i_alpha'][idx - 1]
#                     p_delta = self.err_d['i_delta'][idx - 1]
#
#                 # Zoomed image
#                 zoom = True
#                 img = image(i_regs, self.fits_files[idx], regions_file,
#                             alpha_center, delta_center, idx,
#                             float(self.err_d['i_pm'][0]), p_alpha, p_delta,
#                             zoom)
#
#                 img = img[1:-50, :]  # Removes ds9's legend
#
#                 fig = plt.figure(figsize=self.plot_size, dpi=self.plot_dpi)
#                 ax_1 = fig.add_subplot(1, 1, 1)
#                 ax_1.set_title('Dither {} - Enlarged view '.format(dither))
#
#                 plt.imshow(img)
#
#                 labels = [item.get_text() for item in ax_1.get_xticklabels()]
#
#                 empty_string_labels = [''] * len(labels)
#                 ax_1.set_xticklabels(empty_string_labels)
#
#                 labels = [item.get_text() for item in ax_1.get_yticklabels()]
#
#                 empty_string_labels = [''] * len(labels)
#                 ax_1.set_yticklabels(empty_string_labels)
#
#                 ax_1.set_xlabel('right ascension')
#                 ax_1.set_ylabel('declination')
#
#                 pdf.savefig()
#                 plt.close(fig)
#
#                 # Wide image
#                 zoom = False
#                 img = image(i_regs, self.fits_files[idx], regions_file,
#                             alpha_center, delta_center, idx,
#                             float(self.err_d['i_pm'][0]), p_alpha, p_delta,
#                             zoom)
#
#                 img = img[1:-50, :]  # Removes ds9's legend
#
#                 fig = plt.figure(figsize=self.plot_size, dpi=self.plot_dpi)
#                 ax_1 = fig.add_subplot(1, 1, 1)
#                 ax_1.set_title('Dither {} - Wide view '.format(dither))
#
#                 plt.imshow(img)
#
#                 labels = [item.get_text() for item in ax_1.get_xticklabels()]
#
#                 empty_string_labels = [''] * len(labels)
#                 ax_1.set_xticklabels(empty_string_labels)
#
#                 labels = [item.get_text() for item in ax_1.get_yticklabels()]
#
#                 empty_string_labels = [''] * len(labels)
#                 ax_1.set_yticklabels(empty_string_labels)
#
#                 ax_1.set_xlabel('right ascension')
#                 ax_1.set_ylabel('declination')
#
#                 pdf.savefig()
#                 plt.close(fig)
#
#
