#!/usr/bin/python
# -*- coding: utf-8 -*-

""" Creates a catalog populated of ssos from sextracted catalogs
of single CCDs images.

Versions:
- 0.1: Initial release.

Information:
-

Todo:
    * Get out columns definition. Too long for a single function.
    * Improve variable nomenclature.
    * Explanations are still not so easy to understand.
    * POO implementation?
    * Unit testing.

*GNU Terry Pratchett*
"""
from math import cos, sin
from multiprocessing import Process
from sys import stdout

from astropy.units import degree
from astropy.coordinates import SkyCoord
from numpy import median, pi
from pandas import concat, DataFrame, read_csv, Series

from images_management_elvis import get_borders
from misc_cats import extract_cats_d, create_full_cats, extract_ssos_df
from misc import extract_settings_elvis, check_distance, check_source

__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def create_output_catalog_dict():
    """

    :return: cat_d
    """
    cat_d = {'IDX': [], 'SOURCE': [], 'DITHER': [], 'RA': [], 'DEC': [],
             'VEL': [], 'ABMAG': [], 'THETA': [], 'MAG_AUTO': [],
             'A_IMAGE': [], 'B_IMAGE': [], 'THETA_IMAGE': [],
             'ERRA_IMAGE': [], 'ERRB_IMAGE': [], 'MAGERR_AUTO': [],
             'ERRA_WORLD': [], 'ERRB_WORLD': [], 'ERRTHETA_WORLD': [],
             'CLASS_STAR': []}

    return cat_d


def create_empty_catalog_dict():
    """

    :return: cat_d
    """
    cat_d = {'IDX': [], 'SOURCE': [], 'DITHER': [], 'RA': [], 'DEC': [],
             'VEL': [], 'ABMAG': [], 'THETA': []}

    return cat_d


def propagate_dithers():
    """

    :return:
    """
    ssos_d = create_empty_catalog_dict()
    ssos_df = extract_ssos_df()
    unique_sources = ssos_df['SOURCE']
    idx = 0
    # Move over sources and dithers
    for idx_source, source_ in enumerate(unique_sources):
        source_df = ssos_df[ssos_df['SOURCE'].isin([source_])]
        # Dither 1
        ssos_d['IDX'].append(idx)
        idx_source = source_df['SOURCE'].iloc[0]
        ssos_d['SOURCE'].append(idx_source)
        dither_source = 1
        ssos_d['DITHER'].append(dither_source)
        pm_source = source_df['VEL'].iloc[0]
        ssos_d['VEL'].append(pm_source)
        pa_source = source_df['THETA'].iloc[0]
        ssos_d['THETA'].append(pa_source)
        mag_source = source_df['ABMAG'].iloc[0]
        ssos_d['ABMAG'].append(mag_source)

        dither_time = (565.0 / 2) / 3600.0
        alpha_source = source_df['RA'].iloc[0]
        alpha_increment_per_hour = cos(
            float(pa_source * pi / 180.0)) * float(pm_source)
        alpha_increment_per_dither = alpha_increment_per_hour * dither_time
        alpha_source = alpha_source + (alpha_increment_per_dither / 3600.0)
        ssos_d['RA'].append(alpha_source)
        delta_source = source_df['DEC'].iloc[0]
        delta_increment_per_hour = sin(
            float(pa_source * pi / 180.0)) * float(pm_source)
        delta_increment_per_dither = delta_increment_per_hour * dither_time
        delta_source = delta_source + (delta_increment_per_dither / 3600.0)
        ssos_d['DEC'].append(delta_source)

        for dither_source in range(2, 5, 1):
            # dither_time is equal to fraction of hour
            dither_time = 1003.0/3600.0

            idx += 1
            alpha_increment_per_hour = cos(float(pa_source*pi/180.0)) * float(pm_source)
            alpha_increment_per_dither = alpha_increment_per_hour * dither_time
            alpha_source = alpha_source + (alpha_increment_per_dither / 3600.0)
            delta_increment_per_hour = sin(float(pa_source*pi/180.0)) * float(pm_source)
            delta_increment_per_dither = delta_increment_per_hour * dither_time
            delta_source = delta_source + (delta_increment_per_dither / 3600.0)

            ssos_d['IDX'].append(idx)
            ssos_d['SOURCE'].append(idx_source)
            ssos_d['DITHER'].append(dither_source)
            ssos_d['RA'].append(alpha_source)
            ssos_d['DEC'].append(delta_source)
            ssos_d['VEL'].append(pm_source)
            ssos_d['THETA'].append(pa_source)
            ssos_d['ABMAG'].append(mag_source)

        idx += 1

    sso_cat = DataFrame(ssos_d)

    return sso_cat


def filter_by_position(sso_df):
    """

    sso_clean_df columns:
    - ABMAG
    - DEC
    - DITHER
    - IDX
    - RA
    - SOURCE
    - THETA
    - VEL

    :param sso_df:
    :return: sso_clean_df:
    """
    borders_d = get_borders()

    right_sources = []

    unique_sources = list(set(sso_df['SOURCE']))

    for idx_source_, source_ in enumerate(unique_sources):
        source_df = sso_df[sso_df['SOURCE'].isin([source_])]
        for idx_dither_, row in enumerate(source_df.itertuples(), 1):
            alpha = row.RA
            delta = row.DEC
            source_coords = SkyCoord(ra=alpha * degree, dec=delta * degree,
                                     equinox='J2021.5')
            source_coords_ecliptic = source_coords.barycentrictrueecliptic
            lon_e = float(source_coords_ecliptic.lon.degree)
            lat_e = float(source_coords_ecliptic.lat.degree)
            for ccd_ in borders_d[idx_dither_].keys():
                borders = borders_d[idx_dither_][ccd_]
                lon_comp = borders['below_lon'] < lon_e < borders['above_lon']
                lat_comp = borders['below_lat'] < lat_e < borders['above_lat']
                comp = lon_comp and lat_comp

                if comp:
                    right_sources.append(row.IDX)

    # Removes non visible sources
    sso_clean_df = sso_df[sso_df['IDX'].isin(right_sources)]

    return sso_clean_df


def create_catalog():
    """

    :return:
    """
    cats_d = extract_cats_d()  # extracts dataframes from catalogues
    full_d = create_full_cats(cats_d)  # creates dataframe from CCDs catalogues
    ssos_df = propagate_dithers()
    ssos_clean_df = filter_by_position(ssos_df)

    unique_sources = list(set(ssos_clean_df['SOURCE'].tolist()))
    total_ssos = len(list(set(ssos_clean_df['SOURCE'].tolist())))

    sub_list_size = total_ssos / 10
    sub_list_l = []
    for idx_sub_list in range(0, 10, 1):
        if idx_sub_list != (10 - 1):
            idx_down = sub_list_size * idx_sub_list
            idx_up = sub_list_size * (idx_sub_list + 1)
            sub_list_l.append(unique_sources[idx_down:idx_up])
        else:
            idx_down = sub_list_size * idx_sub_list
            sub_list_l.append(unique_sources[idx_down:])

    areas_j = []
    for idx_l in range(0, 10, 1):
        areas_p = Process(target=create_ssos_catalog_thread,
                          args=(idx_l, sub_list_l[idx_l], ssos_clean_df,
                                full_d))
        areas_j.append(areas_p)
        areas_p.start()

    active_areas = list([job.is_alive() for job in areas_j])
    while True in active_areas:
        active_areas = list([job.is_alive() for job in areas_j])
        pass

    # Merges areas
    # Merges catalogs
    ssos_list = []
    for idx_csv in range(0, 10, 1):
        ssos_ = read_csv('tmp_ssos/ssos_{}.csv'.format(idx_csv),
                         index_col=0)
        ssos_list.append(ssos_)

    ssos_df = concat(ssos_list)

    median_b_image_list = []
    median_mag_auto_list = []
    for source_ in list(set(ssos_df['SOURCE'].tolist())):
        source_df = ssos_df[ssos_df['SOURCE'].isin([source_])]
        median_b_image = median(source_df['B_IMAGE'].tolist())
        median_mag_auto = median(source_df['MAG_AUTO'].tolist())
        for idx in range(source_df['IDX'].size):
            median_b_image_list.append(median_b_image)
            median_mag_auto_list.append(median_mag_auto)

    ssos_df['MEDIAN_B_IMAGE'] = Series(median_b_image_list,
                                       name='MEDIAN_B_IMAGE')
    ssos_df['MEDIAN_MAG_AUTO'] = Series(median_mag_auto_list,
                                        name='MEDIAN_MAG_AUTO')

    ssos_df.to_csv('catalogues_detected/ssos.csv')

    return ssos_df


def create_ssos_catalog_thread(idx_l, sub_list, ssos_df, full_d):
    """

    :param idx_l:
    :param sub_list:
    :param ssos_df:
    :param full_d:
    :return:
    """
    keys = ['ALPHA_J2000', 'DELTA_J2000']

    # Creates an empty catalog for all detected sources
    cat_d = create_output_catalog_dict()
    total_thread = len(sub_list)
    stdout.write('total SSOs {} of thread {}\n'.format(total_thread, idx_l))

    for idx, sso in enumerate(sub_list):
        source_df = ssos_df[ssos_df['SOURCE'].isin([sso])]

        for dither_ in source_df['DITHER'].tolist():
            dither_df = source_df[source_df['DITHER'].isin([dither_])]
            # Gets alpha/delta of actual source
            alpha = float(dither_df['RA'].iloc[0])
            delta = float(dither_df['DEC'].iloc[0])

            o_df = check_source(full_d[dither_], alpha, delta, keys)
            if o_df.empty is not True:
                # Returns the index of the closest found source
                index = check_distance(o_df, alpha, delta)

                o_df = o_df.iloc[[index]]

                idx_ = int(dither_df['IDX'].iloc[0])
                cat_d['IDX'].append(idx_)

                source = int(dither_df['SOURCE'].iloc[0])
                cat_d['SOURCE'].append(source)

                cat_d['DITHER'].append(dither_)

                alpha_j2000 = float(o_df['ALPHA_J2000'].iloc[0])
                cat_d['RA'].append(alpha_j2000)

                delta_j2000 = float(o_df['DELTA_J2000'].iloc[0])
                cat_d['DEC'].append(delta_j2000)

                vel = float(dither_df['VEL'].iloc[0])
                cat_d['VEL'].append(vel)

                abmag = float(dither_df['ABMAG'].iloc[0])
                cat_d['ABMAG'].append(abmag)

                theta = float(dither_df['THETA'].iloc[0])
                cat_d['THETA'].append(theta)

                mag_auto = float(o_df['MAG_AUTO'].iloc[0])
                cat_d['MAG_AUTO'].append(mag_auto)

                a_image = float(o_df['A_IMAGE'].iloc[0])
                cat_d['A_IMAGE'].append(a_image)

                b_image = float(o_df['B_IMAGE'].iloc[0])
                cat_d['B_IMAGE'].append(b_image)

                theta_image = float(o_df['THETA_IMAGE'].iloc[0])
                cat_d['THETA_IMAGE'].append(theta_image)

                erra_image = float(o_df['ERRA_IMAGE'].iloc[0])
                cat_d['ERRA_IMAGE'].append(erra_image)

                errb_image = float(o_df['ERRB_IMAGE'].iloc[0])
                cat_d['ERRB_IMAGE'].append(errb_image)

                magerr_auto = float(o_df['MAGERR_AUTO'].iloc[0])
                cat_d['MAGERR_AUTO'].append(magerr_auto)

                erra_world = float(o_df['ERRA_WORLD'].iloc[0])
                cat_d['ERRA_WORLD'].append(erra_world)

                errb_world = float(o_df['ERRB_WORLD'].iloc[0])
                cat_d['ERRB_WORLD'].append(errb_world)

                errtheta_world = float(o_df['ERRTHETA_WORLD'].iloc[0])
                cat_d['ERRTHETA_WORLD'].append(errtheta_world)

                class_star = float(o_df['CLASS_STAR'].iloc[0])
                cat_d['CLASS_STAR'].append(class_star)

    cat_df = DataFrame(cat_d, columns=['IDX', 'SOURCE', 'DITHER', 'RA', 'DEC',
                                       'VEL', 'ABMAG', 'THETA', 'MAG_AUTO',
                                       'A_IMAGE', 'B_IMAGE', 'THETA_IMAGE',
                                       'ERRA_IMAGE', 'ERRB_IMAGE',
                                       'MAGERR_AUTO', 'ERRA_WORLD',
                                       'ERRB_WORLD', 'ERRTHETA_WORLD',
                                       'CLASS_STAR'])

    cat_df.to_csv('tmp_ssos/ssos_{}.csv'.format(idx_l))


if __name__ == "__main__":
    prfs_dict = extract_settings_elvis()

    catalogue = create_catalog()
