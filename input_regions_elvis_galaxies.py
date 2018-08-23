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
from pandas import concat

from misc import extract_settings_elvis
from misc_cats import extract_galaxies_df


__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


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
