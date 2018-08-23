#!/usr/bin/python
# -*- coding: utf-8 -*-

""" Creates regions populated of stars from sextracted catalogs
of single CCDs images.

Versions:
- 0.1: Initial release. Split from check.py
       Recreated for ELViS analysis pipeline.
- 0.2: Single input catalogue and multiple input catalogues added.
- 0.3: Check method for catalogue creation added.
- 0.4: Easiest way implemented. Rewritten.
- 0.5: Now can saves catalogues by dither.

Information:
- cat: -> hdu_list catalogue
- data: -> Table formatted data
- df: -> dataframe formatted data

Todo:
    * Get out columns definition. Too long for a single function.
    * Creates units tests.
    * Improve variable nomenclature.
    * Explanations are still not so easy to understand.
    * POO implementation?

*GNU Terry Pratchett*
"""
from multiprocessing import Process
from subprocess import Popen
from sys import stdout

from astropy.io import fits
from astropy.table import Table
from pandas import concat, DataFrame, read_csv

from misc_cats import extract_cats_d, create_full_cats
from misc import create_sextractor_dict
from misc import extract_settings_elvis, check_distance, check_source

__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.5"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def extract_inputs_d():
    """

    :return:
    """
    inputs_d = {}

    cat_stars_loc = prfs_dict['references']
    cat_stars = fits.open('{}/cat_stars.fits'.format(cat_stars_loc))
    stars_data = Table(cat_stars[1].data)
    stars_df = stars_data.to_pandas()
    stars_idx = range(0, 28474, 1)  # hardcoded - todo!
    stars_df['IDX'] = stars_idx
    inputs_d['stars'] = stars_df

    cat_galaxies_loc = prfs_dict['references']
    cat_galaxies = fits.open('{}/cat_galaxies.fits'.format(cat_galaxies_loc))
    galaxies_data = Table(cat_galaxies[1].data)
    galaxies_df = galaxies_data.to_pandas()
    galaxies_idx = range(0, 143766, 1)  # hardcoded - todo!
    galaxies_df['IDX'] = galaxies_idx
    inputs_d['galaxies'] = galaxies_df

    return inputs_d


def create_empty_catalog_dict():
    """

    :return: cat_d
    """
    cat_d = {'DITHER': [], 'CATALOG_NUMBER': [], 'X_WORLD': [], 'Y_WORLD': [],
             'MAG_AUTO': [], 'MAGERR_AUTO': [], 'ERRA_WORLD': [],
             'ERRB_WORLD': [], 'ERRTHETA_WORLD': []}

    return cat_d


def create_catalog():
    """

    :return:
    """
    cats_d = extract_cats_d()  # extracts dataframes from catalogues
    full_d = create_full_cats(cats_d)  # creates dataframe from CCDs catalogues
    inputs_d = extract_inputs_d()
    save = True

    unique_sources = inputs_d['stars']['IDX']
    total_stars = inputs_d['stars']['IDX'].size

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
                          args=(idx_l, sub_list_l[idx_l], inputs_d, full_d))
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
    positions_table = concat([stars_df['X_WORLD'],
                              stars_df['Y_WORLD']], axis=1)
    if save:
        positions_table.to_csv('regions_detected/stars.reg', index=False,
                               header=False, sep=" ")

    return stars_df


def create_stars_catalog_thread(idx_l, sub_list, inputs_d, full_d):
    """

    :param idx_l:
    :param sub_list:
    :param inputs_d:
    :param full_d:
    :return:
    """
    save = True
    keys = ['ALPHA_J2000', 'DELTA_J2000']

    cat_d = create_empty_catalog_dict()
    total_thread = len(sub_list)
    stdout.write('total stars {} of thread {}\n'.format(total_thread, idx_l))
    for idx, star in enumerate(sub_list):
        stars_df = inputs_d['stars']
        source_df = stars_df[stars_df['IDX'].isin([star])]
        alpha = source_df['RA2000(Gaia)'].iloc[0]
        delta = source_df['DEC2000(Gaia)'].iloc[0]

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

                erra_world = float(o_df['ERRA_WORLD'].iloc[0])
                source_d['ERRA_WORLD'].append(erra_world)

                errb_world = float(o_df['ERRB_WORLD'].iloc[0])
                source_d['ERRB_WORLD'].append(errb_world)

                errtheta_world = float(o_df['ERRTHETA_WORLD'].iloc[0])
                source_d['ERRTHETA_WORLD'].append(errtheta_world)

        if len(source_d['DITHER']) != 0:
            for key_ in source_d.keys():
                for value_ in source_d[key_]:
                    cat_d[key_].append(value_)

    cat_df = DataFrame(cat_d, columns=['DITHER', 'CATALOG_NUMBER', 'X_WORLD',
                                       'Y_WORLD', 'MAG_AUTO', 'MAGERR_AUTO',
                                       'ERRA_WORLD', 'ERRB_WORLD',
                                       'ERRTHETA_WORLD'])
    if save:
        cat_df.to_csv('tmp_stars/stars_{}.csv'.format(idx_l))


def write_stars_catalog(stars_df):
    """

    :param stars_df:
    :return:
    """
    # Stars catalogue creation
    # todo - create full_coadd.cat
    # todo - launch swarp
    # todo - launch sextractor

    analysis_d, len_dicts = create_sextractor_dict(0, False)

    try:
        test_cat_name = '{}/coadd.cat'.format(prfs_dict['references'])
        test_coadd_cat = fits.open(test_cat_name)
    except IOError:
        # Create full coadded catalogue
        create_coadded_image()
        extract_catalogue(analysis_d)
        test_cat_name = '{}/coadd.cat'.format(prfs_dict['references'])
        test_coadd_cat = fits.open(test_cat_name)

    # c1 = fits.Column(name='NUMBER', format='1J', disp='I10',
    #                  array=stars_df['NUMBER'])
    # Kron-like elliptical aperture magnitude
    c1 = fits.Column(name='MAG_AUTO', format='1E', unit='mag',
                     disp='F8.4', array=stars_df['MAG_AUTO'])
    # RMS error for AUTO magnitude
    c2 = fits.Column(name='MAGERR_AUTO', format='1E', unit='mag',
                     disp='F8.4', array=stars_df['MAGERR_AUTO'])
    # Barycenter position along world x axis
    c3 = fits.Column(name='X_WORLD', format='1D', unit='deg', disp='E18.10',
                     array=stars_df['X_WORLD'])
    # Barycenter position along world y axis
    c4 = fits.Column(name='Y_WORLD', format='1D', unit='deg', disp='E18.10',
                     array=stars_df['Y_WORLD'])
    # World RMS position error along major axis
    c5 = fits.Column(name='ERRA_WORLD', format='1E', unit='deg',
                     disp='G12.7', array=stars_df['ERRA_WORLD'])
    # World RMS position error along minor axis
    c6 = fits.Column(name='ERRB_WORLD', format='1E', unit='deg',
                     disp='G12.7', array=stars_df['ERRB_WORLD'])
    # Error ellipse pos.angle(CCW / world - x)
    c7 = fits.Column(name='ERRTHETA_WORLD', format='1E', unit='deg',
                     disp='F6.2', array=stars_df['ERRTHETA_WORLD'])

    col_defs = fits.ColDefs([c1, c2, c3, c4, c5, c6, c7])

    tb_hdu = fits.BinTableHDU.from_columns(col_defs)

    test_coadd_cat[2] = tb_hdu
    test_coadd_cat[2].header['EXTNAME'] = 'LDAC_OBJECTS'

    newcat_name = '{}/stars_catalogue.cat'.format(prfs_dict['references'])
    test_coadd_cat.writeto(newcat_name, overwrite=True)


def create_coadded_image():
    """

    :return:
    """
    sw_1 = 'swarp CCD_*{}.fits'.format(prfs_dict['fits_dir'])

    cmd = sw_1

    swarp_p = Popen(cmd, shell=True)
    swarp_p.wait()


def extract_catalogue(analysis_dict):
    """

    :param analysis_dict:
    :return:
    """
    coadded_fits = '{}/coadd.fits'.format(prfs_dict['fits_dir'])
    coadded_cat = '{}/coadd.cat'.format(prfs_dict['fits_dir'])
    sx_1 = 'sex -c {} {}'.format(prfs_dict['conf_sex'], coadded_fits)
    sx_2 = ' -CATALOG_NAME {}'.format(coadded_cat)
    sx_3 = ' -PARAMETERS_NAME {}'.format(prfs_dict['params_sex'])
    sx_4 = ' -STARNNW_NAME {}'.format(prfs_dict['neural_sex'])
    sx_5 = ' -DETECT_MINAREA {}'.format(analysis_dict['detect_minarea'])
    sx_6 = ' -DETECT_THRESH {}'.format(analysis_dict['detect_thresh'])
    sx_7 = ' -ANALYSIS_THRESH {}'.format(analysis_dict['analysis_thresh'])
    sx_8 = ' -DEBLEND_NTHRESH {}'.format(analysis_dict['deblend_nthresh'])
    sx_9 = ' -DEBLEND_MINCONT {}'.format(analysis_dict['deblend_mincount'])
    sx_10 = ' -FILTER_NAME {}'.format(analysis_dict['filter'])

    cmd = sx_1 + sx_2 + sx_3 + sx_4 + sx_5 + sx_6 + sx_7 + sx_8 + sx_9 + sx_10

    sextractor_p = Popen(cmd, shell=True)
    sextractor_p.wait()


if __name__ == "__main__":
    prfs_dict = extract_settings_elvis()

    catalogue = create_catalog()
    write_stars_catalog(catalogue)
