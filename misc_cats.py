#!/usr/bin/python
# -*- coding: utf-8 -*-

"""  Utilities for catalogues management.

Todo:
    * Improve log messages
    * Improve usability
"""
from os import listdir

from pandas import read_csv

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
