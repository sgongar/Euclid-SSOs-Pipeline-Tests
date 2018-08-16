#!/usr/bin/python
# -*- coding: utf-8 -*-

"""

Versions:
- 0.1: Initial release. Split from ssos_catalog_creation.py

Information:
- cat: -> hdu_list catalogue
- data: -> Table formatted data
- df: -> dataframe formatted data

Todo:
    *

*GNU Terry Pratchett*

"""
from pandas import read_csv

from misc import extract_settings_elvis

__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


# Reads first catalogue dither
# Propagate over dither

def create_empty_catalog_dict():
    """

    :return: cat_d
    """
    cat_d = {'IDX': [], 'DITHER': [], 'CATALOG_NUMBER': [], 'ALPHA_J2000': [],
             'DELTA_J2000': [], 'PM': [], 'MAG': [], 'PA': []}

    return cat_d


def extract_ssos_df():
    """

    :return:
    """
    columns = ['ALPHA_J2000', 'DELTA_J2000', 'PM', 'PA',
               'MAG', 'MAG', 'MAG', 'MAG']
    # TODO actual catalogues is only valid for dither -> 1
    # TODO create catalogues for dithers 2, 3 and 4
    cat_ssos = read_csv('{}/ssos_cat.txt'.format(prfs_dict['references']),
                        delim_whitespace=True, header=None, names=columns)

    ssos_idx = range(0, cat_ssos['ALPHA_J2000'].size, 1)
    cat_ssos['IDX'] = ssos_idx
    ssos_df = cat_ssos[['IDX', 'ALPHA_J2000', 'DELTA_J2000',
                        'PM', 'PA', 'MAG']]

    print(ssos_df)

    return ssos_df


def propagate_dithers():
    """

    :param inputs_d:
    :return:
    """
    ssos_d = create_empty_catalog_dict()
    ssos_df = extract_ssos_df()
    unique_sources = ssos_df['IDX']
    # Move over sources and dithers
    for idx_sso, sso_ in enumerate(unique_sources):
        source_df = ssos_df[ssos_df['IDX'].isin([sso_])]
        # Dither 1
        idx_source = source_df['IDX'].iloc[0]
        ssos_d['IDX'].append(idx_source)
        dither_source = 1
        ssos_d['DITHER'].append(dither_source)
        alpha_source = source_df['ALPHA_J2000'].iloc[0]
        ssos_d['ALPHA_J2000'].append(alpha_source)
        delta_source = source_df['DELTA_J2000'].iloc[0]
        ssos_d['DELTA_J2000'].append(delta_source)
        pm_source = source_df['PM'].iloc[0]
        ssos_d['PM'].append(pm_source)
        pa_source = source_df['PA'].iloc[0]
        ssos_d['PA'].append(pa_source)
        mag_source = source_df['MAG'].iloc[0]
        ssos_d['MAG'].append(mag_source)


if __name__ == "__main__":
    prfs_dict = extract_settings_elvis()

    propagate_dithers()
