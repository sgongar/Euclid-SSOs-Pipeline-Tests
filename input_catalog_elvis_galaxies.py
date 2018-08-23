#!/usr/bin/python
# -*- coding: utf-8 -*-

""" Creates a catalog populated of galaxies from sextracted catalogs
of single CCDs images.

Versions:
- 0.1: Initial release. Split from stars_catalog_creation.py

Information:
-

Todo:
    * Unit testing.

*GNU Terry Pratchett*
"""
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

    if save:
        galaxies_df.to_csv('catalogues_input/galaxies.csv')


if __name__ == "__main__":
    prfs_dict = extract_settings_elvis()

    catalogue = create_catalog()
