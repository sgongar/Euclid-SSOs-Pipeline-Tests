#!/usr/bin/python
# -*- coding: utf-8 -*-

"""  Utilities for catalogues management.

Todo:
    * Improve log messages
    * Improve usability
"""
from os import listdir

from astropy.io import fits
from astropy.table import Table
from pandas import concat, read_csv

from misc import extract_settings_elvis


__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def get_dither(catalog_n):
    """

    :param catalog_n:
    :return: cat_list
    """
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


def get_cat(ccd):
    """ Returns the catalogue of a chosen CCD.

    :param ccd: The chosen CCD.
    :return: The sextractor catalogue.
    """
    cats = []
    cat_n = 0
    idx = 0

    for x_ in range(1, 7, 1):
        for y_ in range(1, 7, 1):
            for d_ in range(1, 5, 1):
                cat_name = 'x{}_y{}'.format(x_, y_)
                cats.append([cat_name, d_, idx])

                idx += 1

    ccd_position = ccd[4:9]
    dither_n = ccd[11:12]

    for cat_ in cats:
        if str(ccd_position) == cat_[0] and int(dither_n) == cat_[1]:
            cat_n = cat_[2]

    return cat_n


def get_cats(dither):
    """ Returns a list with all catalogues in the catalogues directory of a
    chosen dither.

    :param dither: The chosen dither.
    :return: A list with the catalogues files.
    """
    prfs_d = extract_settings_elvis()  # Gets the preferences dictionary.
    cats_list = []  # Creates an empty list of catalogues.

    # Gets all the catalogues files.
    cats = listdir('{}'.format(prfs_d['fits_dir']))
    for cat_ in cats:
        if cat_[-4:] == '.cat':
            cats_list.append(cat_)

    # Gets the catalogues of the chosen dither.
    list_out = []
    for cat_ in cats_list:
        if cat_[-5:-4] == str(dither):
            list_out.append(cat_)

    return list_out


def gets_data():
    """ Creates an input dictionary. Each key contains SSOs' information
    for each dither.

    :return: input_dict
    """
    # # For now we only have data for dither 1
    input_df = {1: {}, 2: {}, 3: {}, 4: {}}

    for key_ in input_df.keys():
        # Uses clean ones instead total ones
        ssos_cat = 'catalogues_input/cat_clean_ssos_{}.csv'.format(key_)
        input_df[key_]['SSOs'] = read_csv(ssos_cat, index_col=0)
        stars_cat = 'catalogues_detected/stars.csv'
        input_df[key_]['stars'] = read_csv(stars_cat, index_col=0)
        galaxies_cat = 'catalogues_detected/galaxies.csv'
        input_df[key_]['galaxies'] = read_csv(galaxies_cat, index_col=0)

    return input_df


def extract_cats_d():
    """

    :return:
    """
    cats_d = {}
    prfs_dict = extract_settings_elvis()
    for dither in range(1, 5, 1):
        cats_d[dither] = {}
        cats = get_cats(dither)
        for cat_name in cats:
            hdu_list = fits.open('{}/{}'.format(prfs_dict['fits_dir'],
                                                cat_name))
            cat_data = Table(hdu_list[2].data)
            cat_df = cat_data.to_pandas()  # Converts to Pandas format
            cat_number = get_cat(cat_name)  # Gets cat's number from cat's name
            cats_d[dither][cat_name] = cat_df

            cat_list = [cat_number] * cat_df['NUMBER'].size
            cats_d[dither][cat_name]['CATALOG_NUMBER'] = cat_list

    return cats_d


def create_full_cats(cats_d):
    """

    :param cats_d:
    :return:
    """
    full_d = {}

    for dither in range(1, 5, 1):
        dither_l = []
        for key_ in cats_d[dither].keys():
            dither_l.append(cats_d[dither][key_])
        full_d[dither] = concat(dither_l, ignore_index=True)
        full_idx = range(0, full_d[dither]['NUMBER'].size, 1)
        full_d[dither]['IDX'] = full_idx

    return full_d


def extract_stars_df():
    """

    :return:
    """
    prfs_dict = extract_settings_elvis()
    cat_stars_loc = prfs_dict['references']
    cat_stars = fits.open('{}/cat_stars.fits'.format(cat_stars_loc))
    stars_data = Table(cat_stars[1].data)
    stars_df = stars_data.to_pandas()
    stars_idx = range(0, 28474, 1)  # hardcoded - todo!
    stars_df['IDX'] = stars_idx

    return stars_df


def extract_galaxies_df():
    """

    :return:
    """
    prfs_dict = extract_settings_elvis()
    cat_galaxies_loc = prfs_dict['references']
    cat_galaxies = fits.open('{}/cat_galaxies.fits'.format(cat_galaxies_loc))
    galaxies_data = Table(cat_galaxies[1].data)
    galaxies_df = galaxies_data.to_pandas()
    galaxies_idx = range(0, 143766, 1)
    galaxies_df['IDX'] = galaxies_idx

    return galaxies_df


def extract_ssos_df():
    """

    :return:
    """
    prfs_dict = extract_settings_elvis()
    ssos_df = read_csv('{}/ssos_cat.txt'.format(prfs_dict['references']),
                       delim_whitespace=True)
    ssos_source = range(0, ssos_df['RA'].size, 1)
    ssos_df['SOURCE'] = ssos_source

    print(ssos_df['IDX'].size)

    ssos_idx = range(0, 999, 1)
    ssos_df['IDX'] = ssos_idx

    return ssos_df
