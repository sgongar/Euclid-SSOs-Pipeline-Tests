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
from astropy.units import degree
from astropy.coordinates import SkyCoord
from pandas import concat, read_csv, Series

from images_management_elvis import get_borders
from misc import extract_settings_elvis, check_source
from misc_cats import get_cats


__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


class FactorsScampPerformance:

    def __init__(self):
        """ Me da los valores de salida de todos las estrellas  y galaxias
        presentes en filt 3 obtenidos o no

        """
        self.filter_p_number = 3  # First one with enough data for statistics
        self.prfs_d = extract_settings_elvis()

        ccds = True
        filtered = False
        scamp = False

        self.input_df = read_csv('tmp_galaxies/galaxies.csv', index_col=0)
        self.filt_cat = self.gets_filtered_catalog()  # Gets data from filtered

        if ccds:
            cats_d = self.extract_cats()
            self.extract_stats_ccds(cats_d)
        elif filtered:
            pass
            # self.extract_stats_filt()  # not implemented yet
        elif scamp:
            pass
            # self.extract_stats_scamp(input_df)  # not implemented yet
        else:
            pass

    def extract_cats(self):
        """

        :return:
        """
        cats_d = {}
        for dither_ in range(1, 5, 1):
            cat_list = get_cats(dither_)
            cats_d[dither_] = {}
            for cat_ in cat_list:
                cat = fits.open('{}/{}'.format(self.prfs_d['fits_dir'], cat_))
                cat_df = Table(cat[2].data).to_pandas()
                cats_d[dither_][cat_] = cat_df

        return cats_d

    def gets_filtered_catalog(self):
        """

        :return: filtered_cat, filtered catalog
        """
        filter_n = 'filt__{}.csv'.format(self.filter_p_number)
        # filter_n = 'filt__full_1.csv'
        filter_o_n = '{}/{}'.format(self.prfs_d['filtered'], filter_n)

        print('Opens filtered catalogue {}'.format(filter_o_n))
        filtered_cat = read_csv('{}'.format(filter_o_n), index_col=0)

        return filtered_cat

    def extract_stats_ccds(self, cats_d):
        """

        :param cats_d:
        :param filt_cat:
        :return:
        """
        # Unique sources (?)
        test_dict = {'PM': [], 'A_IMAGE': [], 'B_IMAGE': [],
                     'CLASS_STAR': [], 'MAG_ISO': [], 'ELLIPTICITY': []}

        idx = 0
        borders_d = get_borders()

        print('total {}'.format(self.input_df['DITHER'].size))
        for i, row in enumerate(self.input_df.itertuples(), 1):
            print('idx: {}'.format(idx))
            idx += 1
            dither = int(row.DITHER)
            i_alpha = float(row.X_WORLD)
            i_delta = float(row.Y_WORLD)

            pm = ''
            a_image = ''
            b_image = ''
            class_star = ''
            mag_iso = ''
            ellipticity = ''
            test = False

            # Looks for where the source should be
            source_coords = SkyCoord(ra=i_alpha * degree, dec=i_delta * degree,
                                     equinox='J2021.5')
            source_coords_ecliptic = source_coords.barycentrictrueecliptic
            lon_e = float(source_coords_ecliptic.lon.degree)
            lat_e = float(source_coords_ecliptic.lat.degree)
            for ccd_ in borders_d[dither].keys():
                borders = borders_d[dither][ccd_]
                lon_comp = borders['below_lon'] < lon_e < borders['above_lon']
                lat_comp = borders['below_lat'] < lat_e < borders['above_lat']
                comp = lon_comp and lat_comp

                if comp:
                    catalogue = '{}.cat'.format(ccd_[:-5])
                    out_df = check_source(cats_d[dither][catalogue],
                                          i_alpha, i_delta,
                                          keys=['ALPHA_J2000', 'DELTA_J2000'])

                    if out_df.empty is not True:
                        alpha = float(out_df['ALPHA_J2000'].iloc[0])
                        delta = float(out_df['DELTA_J2000'].iloc[0])
                        pm_df = check_source(self.filt_cat, alpha, delta,
                                             keys=['ALPHA_J2000',
                                                   'DELTA_J2000'])
                        if pm_df.empty is not True:
                            test = True
                            pm = float(pm_df['PM'].iloc[0])
                            a_image = float(pm_df['MEDIAN_A_IMAGE'].iloc[0])
                            b_image = float(pm_df['MEDIAN_B_IMAGE'].iloc[0])
                            class_star = float(pm_df['MEDIAN_CLASS_STAR'].iloc[0])
                            mag_iso = float(pm_df['MEDIAN_MAG_ISO'].iloc[0])
                            ellipticity = float(pm_df['MEDIAN_ELLIPTICITY'].iloc[0])

            if test:
                test_dict['PM'].append(pm)
                test_dict['A_IMAGE'].append(a_image)
                test_dict['B_IMAGE'].append(b_image)
                test_dict['CLASS_STAR'].append(class_star)
                test_dict['MAG_ISO'].append(mag_iso)
                test_dict['ELLIPTICITY'].append(ellipticity)

        pm_list = Series(test_dict['PM'], name='PM')
        a_image_list = Series(test_dict['A_IMAGE'], name='A_IMAGE')
        b_image_list = Series(test_dict['B_IMAGE'], name='B_IMAGE')
        class_star_list = Series(test_dict['CLASS_STAR'], name='CLASS_STAR')
        mag_iso_list = Series(test_dict['MAG_ISO'], name='MAG_ISO')
        ellipticity_list = Series(test_dict['ELLIPTICITY'], name='ELLIPTICITY')

        positions_table = concat([pm_list, a_image_list, b_image_list,
                                  class_star_list, mag_iso_list,
                                  ellipticity_list], axis=1)
        positions_table.to_csv('false_catalogues/galaxies_catalogue.csv')
        positions_table.to_csv('true_catalogues/galaxies_catalogue.csv')


    # def extract_stats_filt(self, filt_cat, input_df):
    #     """
    #
    #     :param filt_cat:
    #     :param input_df:
    #     :return:
    #     """
    #     # Unique sources (?)
    #     unique_sources = list(set(input_df['SOURCE'].tolist()))
    #
    #     ok = 0
    #     no = 0
    #     filt_cat_size = filt_cat['SOURCE_NUMBER'].size
    #     total_input_size = input_df['SOURCE'].size
    #     print('sources in input catalogue: {}'.format(total_input_size))
    #     print('sources in filtered catalogue {}'.format(filt_cat_size))
    #
    #     for source_ in unique_sources:
    #         source_df = input_df[input_df['SOURCE'].isin([source_])]
    #         print(source_df.columns)
    #         for i, row in enumerate(source_df.itertuples(), 1):
    #             dither_df = source_df[source_df['DITHER'].isin([row.DITHER])]
    #             i_alpha = float(dither_df['RA'].iloc[0])
    #             i_delta = float(dither_df['DEC'].iloc[0])
    #             out_df = check_source(filt_cat, i_alpha, i_delta,
    #                                   keys=['ALPHA_J2000', 'DELTA_J2000'])
    #
    #             if out_df.empty:
    #                 no += 1
    #             else:
    #                 ok += 1
    #
    #     print('detected: {}'.format(ok))
    #     print('not detected: {}'.format(no))


if __name__ == "__main__":
    FactorsScampPerformance()
