#!/usr/bin/python
# -*- coding: utf-8 -*-

"""

Versions:
- 0.1: Initial release.

Information:
-

Todo:
    *

*GNU Terry Pratchett*
"""
from misc import extract_settings_elvis
from misc_cats import extract_stars_df

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
    stars_df = extract_stars_df()

    if save:
        stars_df.to_csv('catalogues_input/stars.csv')


if __name__ == "__main__":
    prfs_dict = extract_settings_elvis()

    catalogue = create_catalog()
