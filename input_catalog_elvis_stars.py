#!/usr/bin/python
# -*- coding: utf-8 -*-

"""

Versions:
- 0.1: Initial release.

Information:
- cat: -> hdu_list catalogue
- data: -> Table formatted data
- df: -> dataframe formatted data

Todo:
    *

*GNU Terry Pratchett*
"""
from astropy.io import fits
from astropy.table import Table

from misc import extract_settings_elvis

__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def extract_stars_df():
    """

    :return:
    """
    cat_stars_loc = prfs_dict['references']
    cat_stars = fits.open('{}/cat_stars.fits'.format(cat_stars_loc))
    stars_data = Table(cat_stars[1].data)
    stars_df = stars_data.to_pandas()
    stars_idx = range(0, 28474, 1)  # hardcoded - todo!
    stars_df['IDX'] = stars_idx

    return stars_df


def create_catalog():
    """

    :return:
    """
    save = True
    stars_df = extract_stars_df()

    if save:
        stars_df.to_csv('catalogues_input/stars.csv')


if __name__ == "__main__":
    prfs_dict = extract_settings_elvis()

    catalogue = create_catalog()
