# !/usr/bin/python
# -*- coding: utf-8 -*-

""" Gets
   - factors from magnitude bins

Versions:
- 0.1

Todo:
    * Improve log messages

*GNU Terry Pratchett*
"""

from pandas import concat, DataFrame, read_csv, Series

from misc import extract_settings_elvis, setting_logger


__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def main_class_pm():
    mags = ['20-21', '21-22', '22-23', '23-24', '24-25', '25-26']
    catalog_loc = '/pcdisk/holly/sgongora/Dev/Euclid-tests/performance/stats/total.csv'
    catalog = read_csv(catalog_loc, index_col=0)

    for pm_ in prfs_dict['pms']:
        pm_cat = catalog[catalog['pm'].isin([pm_])]
        for mag_ in mags:
            mag_cat = pm_cat[pm_cat['mag'].isin([mag_])]
            if pm_cat.empty is False and mag_cat.empty is False:
                mag_cat.to_csv('stats_{}_{}.csv'.format(mag_, pm_))


if __name__ == "__main__":
    prfs_dict = extract_settings_elvis()

    main_class_pm()
