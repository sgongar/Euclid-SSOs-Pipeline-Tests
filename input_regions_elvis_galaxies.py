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

    positions_table = concat([galaxies_df['ra'], galaxies_df['dec']], axis=1)
    if save:
        positions_table.to_csv('regions_input/galaxies.reg', index=False,
                               header=False, sep=" ")


if __name__ == "__main__":
    prfs_dict = extract_settings_elvis()

    catalogue = create_catalog()
