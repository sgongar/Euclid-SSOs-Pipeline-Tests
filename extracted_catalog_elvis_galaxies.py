#!/usr/bin/python
# -*- coding: utf-8 -*-

""" Creates a catalogue populated of galaxies from sextracted catalogues
of single CCDs images.

Versions:
- 0.1: Initial release. Split from stars_catalog_creation.py

Information:
-

Behaviour:


Todo:
    * POO implementation.
    * Unit testing.

*GNU Terry Pratchett*
"""
from multiprocessing import Process
from sys import stdout

from pandas import concat, DataFrame, read_csv

from misc import check_distance, check_source, extract_settings_elvis
from misc_cats import extract_cats_d, create_full_cats, extract_galaxies_df

__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def create_empty_catalog_dict():
    """

    :return: cat_d
    """
    cat_d = {'DITHER': [], 'CATALOG_NUMBER': [], 'X_WORLD': [], 'Y_WORLD': [],
             'MAG_AUTO': [], 'A_IMAGE': [], 'B_IMAGE': [], 'THETA_IMAGE': [],
             'ERRA_IMAGE': [], 'ERRB_IMAGE': [], 'MAGERR_AUTO': [],
             'ERRA_WORLD': [], 'ERRB_WORLD': [], 'ERRTHETA_WORLD': [],
             'CLASS_STAR': []}

    return cat_d


def create_catalog():
    """

    :return:
    """
    cats_d = extract_cats_d()  # extracts dataframes from catalogues
    full_d = create_full_cats(cats_d)  # creates dataframe from CCDs catalogues
    galaxies_df = extract_galaxies_df()
    save = True

    unique_sources = galaxies_df['IDX']
    total_galaxies = galaxies_df['IDX'].size

    sub_list_size = total_galaxies / 18

    sub_list_l = []
    for idx_sub_list in range(0, 18, 1):
        if idx_sub_list != (18 - 1):
            idx_down = sub_list_size * idx_sub_list
            idx_up = sub_list_size * (idx_sub_list + 1)
            sub_list_l.append(unique_sources[idx_down:idx_up])
        else:
            idx_down = sub_list_size * idx_sub_list
            sub_list_l.append(unique_sources[idx_down:])

    areas_j = []
    for idx_l in range(0, 18, 1):
        areas_p = Process(target=create_galaxies_catalog_thread,
                          args=(idx_l, sub_list_l[idx_l], galaxies_df, full_d))
        areas_j.append(areas_p)
        areas_p.start()

    active_areas = list([job.is_alive() for job in areas_j])
    while True in active_areas:
        active_areas = list([job.is_alive() for job in areas_j])
        pass

    # Merges areas
    # Merges catalogs
    galaxies_list = []
    for idx_csv in range(0, 18, 1):
        galaxies_ = read_csv('tmp_galaxies/galaxies_{}.csv'.format(idx_csv),
                             index_col=0)
        galaxies_list.append(galaxies_)

    galaxies_df = concat(galaxies_list)

    if save:
        galaxies_df.to_csv('catalogues_detected/galaxies.csv')

    return galaxies_df


def create_galaxies_catalog_thread(idx_l, sub_list, galaxies_df, full_d):
    """

    :param idx_l:
    :param sub_list:
    :param galaxies_df:
    :param full_d:
    :return:
    """
    save = True
    keys = ['ALPHA_J2000', 'DELTA_J2000']

    cat_d = create_empty_catalog_dict()
    total_thread = len(sub_list)
    stdout.write('total galaxies {} of thread {}\n'.format(total_thread,
                                                           idx_l))
    for idx, galaxy in enumerate(sub_list):
        source_df = galaxies_df[galaxies_df['IDX'].isin([galaxy])]
        alpha = source_df['ra'].iloc[0]
        delta = source_df['dec'].iloc[0]

        source_d = create_empty_catalog_dict()
        for dither in range(1, 5, 1):
            o_df = check_source(full_d[dither], alpha, delta, keys)
            if o_df.empty is not True:
                # Returns the index of the closest found source
                index = check_distance(o_df, alpha, delta)
                o_df = o_df.iloc[[index]]
                source_d['DITHER'].append(dither)

                catalog_number = int(o_df['CATALOG_NUMBER'].iloc[0])
                source_d['CATALOG_NUMBER'].append(catalog_number)

                x_world = float(o_df['X_WORLD'].iloc[0])
                source_d['X_WORLD'].append(x_world)

                y_world = float(o_df['Y_WORLD'].iloc[0])
                source_d['Y_WORLD'].append(y_world)

                mag_auto = float(o_df['MAG_AUTO'].iloc[0])
                source_d['MAG_AUTO'].append(mag_auto)

                magerr_auto = float(o_df['MAGERR_AUTO'].iloc[0])
                source_d['MAGERR_AUTO'].append(magerr_auto)

                a_image = float(o_df['A_IMAGE'].iloc[0])
                source_d['A_IMAGE'].append(a_image)

                b_image = float(o_df['B_IMAGE'].iloc[0])
                source_d['B_IMAGE'].append(b_image)

                theta_image = float(o_df['THETA_IMAGE'].iloc[0])
                source_d['THETA_IMAGE'].append(theta_image)

                erra_image = float(o_df['ERRA_IMAGE'].iloc[0])
                source_d['ERRA_IMAGE'].append(erra_image)

                errb_image = float(o_df['ERRB_IMAGE'].iloc[0])
                source_d['ERRB_IMAGE'].append(errb_image)

                erra_world = float(o_df['ERRA_WORLD'].iloc[0])
                source_d['ERRA_WORLD'].append(erra_world)

                errb_world = float(o_df['ERRB_WORLD'].iloc[0])
                source_d['ERRB_WORLD'].append(errb_world)

                errtheta_world = float(o_df['ERRTHETA_WORLD'].iloc[0])
                source_d['ERRTHETA_WORLD'].append(errtheta_world)

                class_star = float(o_df['CLASS_STAR'].iloc[0])
                source_d['CLASS_STAR'].append(class_star)

        if len(source_d['DITHER']) != 0:
            for key_ in source_d.keys():
                for value_ in source_d[key_]:
                    cat_d[key_].append(value_)

    # for key_ in cat_d.keys():
    #     print(key_, len(cat_d[key_]))

    cat_df = DataFrame(cat_d, columns=['DITHER', 'CATALOG_NUMBER',
                                       'X_WORLD', 'Y_WORLD', 'MAG_AUTO',
                                       'MAGERR_AUTO', 'A_IMAGE', 'B_IMAGE',
                                       'THETA_IMAGE', 'ERRA_IMAGE',
                                       'ERRB_IMAGE', 'ERRA_WORLD',
                                       'ERRB_WORLD', 'ERRTHETA_WORLD',
                                       'CLASS_STAR'])
    if save:
        cat_df.to_csv('tmp_galaxies/galaxies_{}.csv'.format(idx_l))


if __name__ == "__main__":
    prfs_dict = extract_settings_elvis()

    catalogue = create_catalog()
