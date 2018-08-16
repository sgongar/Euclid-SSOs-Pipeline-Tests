#!/usr/bin/python
# -*- coding: utf-8 -*-

""" Creates a bunch of region files from sextractor catalogues.
Output can be save as CCD regions file or dither regions file.

Versions:
- 0.1: Initial release. Split from ssos_catalog_creation.py
- 0.2: Rename. Description improved

Todo:
    * Unit testing.

*GNU Terry Pratchett*
"""
from astropy.io import fits
from astropy.table import Table
from pandas import concat, Series

from misc import extract_settings_elvis
from misc_cats import get_cats


__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.2"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def create_regions_by_ccd():
    """ Creates a regions file for each CCD.

    :return: True if everything is alright.
    """
    for dither_ in range(1, 5, 1):
        cat_list = get_cats(dither_)
        for cat_ in cat_list:
            catalog = fits.open('{}/{}'.format(prfs_d['fits_dir'], cat_))
            catalog_data = Table(catalog[2].data).to_pandas()

            alpha_list = Series(catalog_data['ALPHA_J2000'].tolist(),
                                name='ALPHA_J2000')
            delta_list = Series(catalog_data['DELTA_J2000'].tolist(),
                                name='DELTA_J2000')

            positions_table = concat([alpha_list, delta_list], axis=1)
            positions_table.to_csv('{}.reg'.format(cat_[:-4]),
                                   index=False, header=False, sep=" ")

    return True


def create_regions_by_dither():
    """ Creates a regions file for each dither.

    :return: True if everything is alright.
    """
    for dither_ in range(1, 5, 1):
        alpha_list = []
        delta_list = []
        cat_list = get_cats(dither_)
        for cat_ in cat_list:
            catalog = fits.open('{}/{}'.format(prfs_d['fits_dir'], cat_))
            catalog_data = Table(catalog[2].data).to_pandas()

            alpha_tmp = catalog_data['ALPHA_J2000'].tolist()
            delta_tmp = catalog_data['DELTA_J2000'].tolist()

            for alpha_ in alpha_tmp:
                alpha_list.append(alpha_)
            for delta_ in delta_tmp:
                delta_list.append(delta_)

        alpha_series = Series(alpha_list, name='ALPHA_J2000')
        delta_series = Series(delta_list, name='DELTA_J2000')
        positions_table = concat([alpha_series, delta_series], axis=1)
        positions_table.to_csv('regions/extracted_{}.reg'.format(dither_),
                               index=False, header=False, sep=" ")

    return True


if __name__ == "__main__":
    prfs_d = extract_settings_elvis()

    if create_regions_by_dither():
        pass
    else:
        raise Exception
