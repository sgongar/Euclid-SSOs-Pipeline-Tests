# !/usr/bin/python
# -*- coding: utf-8 -*-

"""

Versions:
- 0.1

Todo:
    * Improve log messages

*GNU Terry Pratchett*
"""

from astropy.io import fits
from astropy.table import Table
from pandas import concat, read_csv, Series

from misc import extract_settings_elvis, check_source, setting_logger
from misc import get_norm_mag, get_norm_speed
from misc_cats import get_dither


__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


class PlotFalseMovement:
    def __init__(self):
        """
        """
        self.prfs_d = extract_settings_elvis()
        self.read_catalog()

    def read_catalogue(self):
        """

        :return: catalogue
        """
        print(self.prfs_d)


if __name__ == "__main__":
    PlotFalseMovement()
