#!/usr/bin/python
# -*- coding: utf-8 -*-

""" Creates a catalog populated of galaxies from sextracted catalogs
of single CCDs images.

Versions:
- 0.1: Initial release. Split from stars_catalog_creation.py

Information:
- cat: -> hdu_list catalogue
- data: -> Table formatted data
- df: -> dataframe formatted data

Todo:
    * Unit testing.

*GNU Terry Pratchett*
"""
from astropy.io import fits
from astropy.table import Table
from pandas import concat

from misc import extract_settings_elvis

__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def extract_galaxies_df():
    """

    :return:
    """
    cat_galaxies_loc = prfs_dict['references']
    cat_galaxies = fits.open('{}/cat_galaxies.fits'.format(cat_galaxies_loc))
    galaxies_data = Table(cat_galaxies[1].data)
    galaxies_df = galaxies_data.to_pandas()
    galaxies_idx = range(0, 143766, 1)
    galaxies_df['IDX'] = galaxies_idx

    return galaxies_df


def create_catalog():
    """

    :return:
    """
    save = True
    galaxies_df = extract_galaxies_df()

    if save:
        galaxies_df.to_csv('catalogues_input/galaxies.csv')


if __name__ == "__main__":
    prfs_dict = extract_settings_elvis()

    catalogue = create_catalog()
