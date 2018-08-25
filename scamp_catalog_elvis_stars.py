#!/usr/bin/python
# -*- coding: utf-8 -*-

""" Creates a catalog populated of stars from sextracted catalogs
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
from multiprocessing import Process
from sys import stdout

from pandas import concat, DataFrame, read_csv

from misc_cats import extract_cats_d, create_full_cats, extract_stars_df
from misc_cats import create_scamp_df
from misc import extract_settings_elvis, check_distance, check_source

__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def get_cat(self, cat_n):
    """

    :param cat_n:
    :return: cat_file
    """
    cats = ['empty_cat']

    for x_ in range(1, 7, 1):
        for y_ in range(1, 7, 1):
            for d_ in range(1, 5, 1):
                cat_name = 'CCD_x{}_y{}_d{}.cat'.format(x_, y_, d_)
                cats.append(cat_name)

    return cats[cat_n]


def test_function():
    """
    """
    from astropy.io import fits
    from astropy.table import Table

    prfs_dict = extract_settings_elvis()

    cats_number = 144
    cat_d = {}
    for cat_n in range(1, cats_number + 1, 1):
        cat_file = get_cat(cat_n)
        cat_data = fits.open('{}/{}'.format(prfs_dict['fits_dir'],
                                            cat_file))

        ccd_df = Table(cat_data[2].data)
        # self.logger.debug('CCD catalog {} to Pandas'.format(cat_n))
        cat_d[cat_n] = ccd_df.to_pandas()

    print(cat_n)


def create_empty_catalog_dict():
    """

    :return: cat_d
    """
    cat_d = {'DITHER': [], 'CATALOG_NUMBER': [], 'X_WORLD': [], 'Y_WORLD': [],
             'MAG_AUTO': [], 'A_IMAGE': [], 'B_IMAGE': [], 'THETA_IMAGE': [],
             'ERRA_IMAGE': [], 'ERRB_IMAGE': [], 'MAGERR_AUTO': [],
             'ERRA_WORLD': [], 'ERRB_WORLD': [], 'ERRTHETA_WORLD': [],
             'CLASS_STAR': [], 'PM': []}

    return cat_d


def create_catalog():
    """

    :return:
    """
    test_function()

    raise Exception
    stars_df = extract_stars_df()
    cats_d = extract_cats_d()  # extracts dataframes from catalogues
    full_d = create_full_cats(cats_d)  # creates dataframe from CCDs catalogues
    scamp_df = create_scamp_df()

    print(scamp_df)

    unique_sources = stars_df['IDX']
    total_stars = stars_df['IDX'].size

    sub_list_size = total_stars / 18

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
        areas_p = Process(target=create_stars_catalog_thread,
                          args=(idx_l, sub_list_l[idx_l],
                                stars_df, full_d, scamp_df))
        areas_j.append(areas_p)
        areas_p.start()

    active_areas = list([job.is_alive() for job in areas_j])
    while True in active_areas:
        active_areas = list([job.is_alive() for job in areas_j])
        pass

    # Merges areas
    # Merges catalogs
    stars_list = []
    for idx_csv in range(0, 18, 1):
        stars_ = read_csv('tmp_stars/stars_{}.csv'.format(idx_csv),
                          index_col=0)
        stars_list.append(stars_)

    stars_df = concat(stars_list)
    stars_df.to_csv('catalogues_detected/stars.csv')

    return stars_df


def create_stars_catalog_thread(idx_l, sub_list, stars_df, full_d, scamp_df):
    """

    :param idx_l:
    :param sub_list:
    :param stars_df:
    :param full_d:
    :return:
    """
    keys = ['ALPHA_J2000', 'DELTA_J2000']

    cat_d = create_empty_catalog_dict()
    total_thread = len(sub_list)
    stdout.write('total stars {} of thread {}\n'.format(total_thread, idx_l))
    for idx, star in enumerate(sub_list):
        source_df = stars_df[stars_df['IDX'].isin([star])]
        alpha = source_df['RA2000(Gaia)'].iloc[0]
        delta = source_df['DEC2000(Gaia)'].iloc[0]

        source_d = create_empty_catalog_dict()
        for dither in range(1, 5, 1):
            sex_df = check_source(full_d[dither], alpha, delta, keys)
            scamp_df = check_source(df, alpha, delta, keys)
            # Should check source in Scamp too!
            if sex_df.empty is not True and scamp_df.empty is not True:
                # Returns the index of the closest found source
                index = check_distance(sex_df, alpha, delta)
                sex_df = sex_df.iloc[[index]]

                source_d['DITHER'].append(dither)

                catalog_number = int(sex_df['CATALOG_NUMBER'].iloc[0])
                source_d['CATALOG_NUMBER'].append(catalog_number)

                x_world = float(sex_df['X_WORLD'].iloc[0])
                source_d['X_WORLD'].append(x_world)

                y_world = float(sex_df['Y_WORLD'].iloc[0])
                source_d['Y_WORLD'].append(y_world)

                mag_auto = float(sex_df['MAG_AUTO'].iloc[0])
                source_d['MAG_AUTO'].append(mag_auto)

                magerr_auto = float(sex_df['MAGERR_AUTO'].iloc[0])
                source_d['MAGERR_AUTO'].append(magerr_auto)

                a_image = float(sex_df['A_IMAGE'].iloc[0])
                source_d['A_IMAGE'].append(a_image)

                b_image = float(sex_df['B_IMAGE'].iloc[0])
                source_d['B_IMAGE'].append(b_image)

                theta_image = float(sex_df['THETA_IMAGE'].iloc[0])
                source_d['THETA_IMAGE'].append(theta_image)

                erra_image = float(sex_df['ERRA_IMAGE'].iloc[0])
                source_d['ERRA_IMAGE'].append(erra_image)

                errb_image = float(sex_df['ERRB_IMAGE'].iloc[0])
                source_d['ERRB_IMAGE'].append(errb_image)

                erra_world = float(sex_df['ERRA_WORLD'].iloc[0])
                source_d['ERRA_WORLD'].append(erra_world)

                errb_world = float(sex_df['ERRB_WORLD'].iloc[0])
                source_d['ERRB_WORLD'].append(errb_world)

                errtheta_world = float(sex_df['ERRTHETA_WORLD'].iloc[0])
                source_d['ERRTHETA_WORLD'].append(errtheta_world)

                class_star = float(sex_df['CLASS_STAR'].iloc[0])
                source_d['CLASS_STAR'].append(class_star)

        if len(source_d['DITHER']) != 0:
            for key_ in source_d.keys():
                for value_ in source_d[key_]:
                    cat_d[key_].append(value_)

    cat_df = DataFrame(cat_d, columns=['DITHER', 'CATALOG_NUMBER',
                                       'X_WORLD', 'Y_WORLD', 'MAG_AUTO',
                                       'MAGERR_AUTO', 'A_IMAGE', 'B_IMAGE',
                                       'THETA_IMAGE', 'ERRA_IMAGE',
                                       'ERRB_IMAGE', 'ERRA_WORLD',
                                       'ERRB_WORLD', 'ERRTHETA_WORLD',
                                       'CLASS_STAR'])
    cat_df.to_csv('tmp_stars/stars_{}.csv'.format(idx_l))


if __name__ == "__main__":
    prfs_dict = extract_settings_elvis()

    catalogue = create_catalog()
