# !/usr/bin/python
# -*- coding: utf-8 -*-

""" Gets
    - 'median_a_image'
    - 'median_erra_image'
    - 'median_b_image'
    - 'median_errb_image'
    - 'median_class_star'
    - 'ellipticity'
    - 'median_mag_iso'
    - 'median_magerr_iso'
    - 'median_flux_iso'
    from scamp's output. Saves them to different csv files.

Versions:
- 0.1

Todo:
    * Improve log messages

*GNU Terry Pratchett*
"""
from multiprocessing import Process
from os import getcwd, remove

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from pandas import concat, DataFrame, read_csv

from misc import extract_settings_elvis, check_source, setting_logger

__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def get_dither(catalog_n):
    """

    :param catalog_n:
    :return: dither_n
    """
    dither_n = 0

    if catalog_n <= 36:
        dither_n = 1
    elif 36 < catalog_n <= 72:
        dither_n = 2
    elif 72 < catalog_n <= 108:
        dither_n = 3
    elif 108 < catalog_n <= 144:
        dither_n = 4

    return dither_n


def create_output_dicts():
    # prfs_d = extract_settings_elvis()

    # Creates a dictionary
    ssos_d = {'median_a_image': [], 'median_erra_image': [],
              'median_b_image': [], 'median_errb_image': [],
              'median_class_star': [],
              'ellipticity': [], 'median_mag_iso': [],
              'median_magerr_iso': [], 'median_flux_iso': [],
              'pm': [], 'pmerr': [], 'fwhm_image': []}
    # for idx in range(0, len(prfs_d['pms']), 1):
    #     ssos_d['median_mag_iso'].append([])
    #     ssos_d['median_magerr_iso'].append([])
    #     ssos_d['median_a_image'].append([])
    #     ssos_d['median_erra_image'].append([])
    #     ssos_d['median_b_image'].append([])
    #     ssos_d['median_errb_image'].append([])
    #     ssos_d['median_class_star'].append([])
    #     ssos_d['ellipticity'].append([])
    #     ssos_d['output_pm'].append([])
    #     ssos_d['median_flux_iso'].append([])

    stars_d = {'median_a_image': [], 'median_erra_image': [],
               'median_b_image': [], 'median_errb_image': [],
               'median_class_star': [], 'ellipticity': [],
               'median_mag_iso': [], 'median_magerr_iso': [],
               'median_flux_iso': [], 'pm': [], 'pmerr': [], 'fwhm_image': []}
    # for idx in range(0, len(prfs_d['pms']), 1):
    #     stars_d['median_mag_iso'].append([])
    #     stars_d['median_magerr_iso'].append([])
    #     stars_d['median_a_image'].append([])
    #     stars_d['median_erra_image'].append([])
    #     stars_d['median_b_image'].append([])
    #     stars_d['median_errb_image'].append([])
    #     stars_d['median_class_star'].append([])
    #     stars_d['ellipticity'].append([])
    #     stars_d['output_pm'].append([])
    #     stars_d['median_flux_iso'].append([])

    galaxies_d = {'median_a_image': [], 'median_erra_image': [],
                  'median_b_image': [], 'median_errb_image': [],
                  'median_class_star': [], 'ellipticity': [],
                  'median_mag_iso': [], 'median_magerr_iso': [],
                  'median_flux_iso': [], 'pm': [], 'pmerr': [],
                  'fwhm_image': []}
    # for idx in range(0, len(prfs_d['pms']), 1):
    #     galaxies_d['median_mag_iso'].append([])
    #     galaxies_d['median_magerr_iso'].append([])
    #     galaxies_d['median_a_image'].append([])
    #     galaxies_d['median_erra_image'].append([])
    #     galaxies_d['median_b_image'].append([])
    #     galaxies_d['median_errb_image'].append([])
    #     galaxies_d['median_class_star'].append([])
    #     galaxies_d['ellipticity'].append([])
    #     galaxies_d['output_pm'].append([])
    #     galaxies_d['median_flux_iso'].append([])

    lost_d = {'median_a_image': [], 'median_erra_image': [],
              'median_b_image': [], 'median_errb_image': [],
              'median_class_star': [],
              'ellipticity': [], 'median_mag_iso': [],
              'median_magerr_iso': [], 'median_flux_iso': [],
              'pm': [], 'pmerr': [], 'fwhm_image': []}

    output_d = {'ssos': ssos_d, 'stars': stars_d, 'galaxies': galaxies_d,
                'lost': lost_d}

    return output_d


class TotalScampPerformance:

    def __init__(self):
        """ Me da los valores de salida de todos las estrellas  y galaxias
        presentes en filt 3 obtenidos o no

        """
        self.filter_p_number = 3  # First one with enough data for statistics
        self.prfs_d = extract_settings_elvis()
        self.data_d = create_output_dicts()

        self.save = True

        logger_name = 'scamp_performance'  # Set as desired
        self.logger = setting_logger(self.prfs_d, logger_name)

        filt_cat = self.gets_filtered_catalog()  # Gets data from filtered
        input_df = self.gets_data()  # Gets data from catalogs

        unique_sources = list(set(filt_cat['SOURCE_NUMBER'].tolist()))
        sub_list_size = len(unique_sources) / self.prfs_d['cores_number']

        sub_list_l = []
        for idx_sub_list in range(0, self.prfs_d['cores_number'], 1):
            if idx_sub_list != (self.prfs_d['cores_number'] - 1):
                idx_down = sub_list_size * idx_sub_list
                idx_up = sub_list_size * (idx_sub_list + 1)
                sub_list_l.append(unique_sources[idx_down:idx_up])
            else:
                idx_down = sub_list_size * idx_sub_list
                sub_list_l.append(unique_sources[idx_down:])

        areas_j = []
        for idx_l in range(0, self.prfs_d['cores_number'], 1):
            areas_p = Process(target=self.splits_data,
                              args=(idx_l, sub_list_l[idx_l],
                                    filt_cat, input_df,))
            areas_j.append(areas_p)
            areas_p.start()

        active_areas = list([job.is_alive() for job in areas_j])
        while True in active_areas:
            active_areas = list([job.is_alive() for job in areas_j])
            pass

        # Merges catalogues
        for key_ in ['stars', 'galaxies', 'ssos', 'lost']:
            csv_list = []
            for idx_csv in range(0, self.prfs_d['cores_number'], 1):
                csv_ = read_csv('cat_{}_{}.csv'.format(idx_csv, key_),
                                index_col=0)
                csv_list.append(csv_)

            full_df = concat(csv_list)
            full_df.to_csv('full_cats/cat_{}.csv'.format(key_))

        # Removes old catalogues
        for key_ in ['stars', 'galaxies', 'ssos', 'lost']:
            for idx_csv in range(0, self.prfs_d['cores_number'], 1):
                remove('{}/cat_{}_{}.csv'.format(getcwd(), idx_csv, key_))

    def gets_filtered_catalog(self):
        """

        :return: filtered_cat, filtered catalog
        """
        filter_n = 'filt__{}.csv'.format(self.filter_p_number)
        filter_o_n = '{}/{}'.format(self.prfs_d['filtered'], filter_n)

        print('Opens filtered catalogue {}'.format(filter_o_n))
        filtered_cat = read_csv('{}'.format(filter_o_n), index_col=0)

        return filtered_cat

    def gets_data(self):
        """ Creates an input dictionary. Each key contains SSOs' information
        for each dither.

        :return: input_dict
        """
        # # For now we only have data for dither 1
        input_df = {1: {}, 2: {}, 3: {}, 4: {}}

        for key_ in input_df.keys():
            ssos_cat = 'cats/cat_clean_ssos_{}.csv'.format(key_)
            input_df[key_]['SSOs'] = read_csv(ssos_cat, index_col=0)
            stars_cat = 'tmp_stars/stars.csv'
            input_df[key_]['stars'] = read_csv(stars_cat, index_col=0)
            galaxies_cat = 'tmp_galaxies/galaxies.csv'
            input_df[key_]['galaxies'] = read_csv(galaxies_cat, index_col=0)

        return input_df

    def splits_data(self, idx_l, unique_sources, filt_cat, input_df):
        """ Splits filtered catalogue due object type.
        This method creates a dictionary with three different keys.
        Each one (SSOs, stars, galaxies) it is a Dataframe with all valuable
        data.

        :param unique_sources:
        :param filt_cat:
        :param input_df:
        :return:
        """
        print('Creating new catalogues from filtered catalogue due type')
        print('Total sources: {} of thread: {}'.format(len(unique_sources),
                                                       idx_l))
        for source_ in unique_sources:
            source_df = filt_cat[filt_cat['SOURCE_NUMBER'].isin([source_])]

            for i, row in enumerate(source_df.itertuples(), 1):
                dither_n = get_dither(int(row.CATALOG_NUMBER))
                # Checks object type
                alpha = source_df['ALPHA_J2000'].iloc[0]
                delta = source_df['DELTA_J2000'].iloc[0]
                keys = ['RA', 'DEC']
                test_sso = check_source(input_df[dither_n]['SSOs'],
                                        alpha, delta, keys)
                keys = ['X_WORLD', 'Y_WORLD']
                test_star = check_source(input_df[dither_n]['stars'],
                                         alpha, delta, keys)

                keys = ['X_WORLD', 'Y_WORLD']
                test_galaxy = check_source(input_df[dither_n]['galaxies'],
                                           alpha, delta, keys)
                if test_sso.empty is not True:
                    key_ = 'ssos'
                elif test_star.empty is not True:
                    key_ = 'stars'
                elif test_galaxy.empty is not True:
                    key_ = 'galaxies'
                else:
                    key_ = 'lost'

                a_image = row.MEDIAN_A_IMAGE
                self.data_d[key_]['median_a_image'].append(a_image)
                erra_image = row.MEDIAN_ERRA_IMAGE
                self.data_d[key_]['median_erra_image'].append(erra_image)
                b_image = row.MEDIAN_B_IMAGE
                self.data_d[key_]['median_b_image'].append(b_image)
                errb_image = row.MEDIAN_ERRB_IMAGE
                self.data_d[key_]['median_errb_image'].append(errb_image)
                class_star = row.MEDIAN_CLASS_STAR
                self.data_d[key_]['median_class_star'].append(class_star)
                ellipticity = row.ELLIPTICITY
                self.data_d[key_]['ellipticity'].append(ellipticity)
                mag_iso = row.MEDIAN_MAG_ISO
                self.data_d[key_]['median_mag_iso'].append(mag_iso)
                magerr_iso = row.MEDIAN_MAGERR_ISO
                self.data_d[key_]['median_magerr_iso'].append(magerr_iso)
                flux_iso = row.MEDIAN_FLUX_ISO
                self.data_d[key_]['median_flux_iso'].append(flux_iso)
                pm = row.PM
                self.data_d[key_]['pm'].append(pm)
                pmerr = row.PMERR
                self.data_d[key_]['pmerr'].append(pmerr)
                fwhm_image = row.FWHM_IMAGE
                self.data_d[key_]['fwhm_image'].append(fwhm_image)

        for type_key in self.data_d.keys():
            data_df = DataFrame(self.data_d[type_key])
            data_df.to_csv('cat_{}_{}.csv'.format(idx_l, type_key))


class PlotTotalScampPerformance:

    def __init__(self):
        """
        Read data from files
        Creates a page per data
        """
        self.dpi = 100
        self.keys_to_plot = ['median_a_image', 'median_b_image',
                             'median_class_image', 'ellipticity',
                             'median_mag_iso', 'median_flux_iso']
        self.data_d = {}

        objects = ['cat_ssos', 'cat_stars', 'cat_galaxies']
        for object_ in objects:
            self.read_data(object_)

        self.plot()

    def read_data(self, object_):
        """

        :param object_:
        :return:
        """
        self.data_d[object_] = read_csv('{}.csv'.format(object_), index_col=0)

    def plot(self):
        """

        :return:
        """
        for key_ in self.data_d.keys():
            pdf_name = '{}.pdf'.format(key_)
            data_d = self.data_d[key_]

            with PdfPages(pdf_name) as pdf:
                # MAG_ISO Galaxies
                fig_1 = plt.figure(figsize=(16.53, 11.69), dpi=self.dpi)
                ax_1 = fig_1.add_subplot(1, 1, 1)
                ax_1.set_title('MAG_ISO - Galaxies')

                ax_1.scatter(data_d['median_b_image'],
                             data_d['median_a_image'],
                             label='a_image', c='b')

                ax_1.legend(loc=4)
                ax_1.grid(True)

                pdf.savefig()  # saves figure
                plt.clf()  # clear current figure
                plt.close(fig_1)  # removes figure


if __name__ == "__main__":
    TotalScampPerformance()
    # PlotTotalScampPerformance()
