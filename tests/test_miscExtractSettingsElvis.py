#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Python script

Versions:

Problems:
    * setup has been ignored

Todo:
    * Improve log messages
    *

"""
import os
import sys
from mock import MagicMock

from unittest import TestCase, main

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from errors import BadSettings

import misc


__author__ = "Samuel Gongora-Garcia"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Gongora-Garcia"]
__version__ = "0.1"
__maintainer__ = "Samuel Gongora-Garcia"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


class TestExtractSettingsElvis(TestCase):
    """

    """

    def create_right_settings_file(self):
        """
        @return:
        """
        testFile = open(".settings_ELViS.ini", "w")
        testFile.write("[Version]\n"
                       "centos_version: /media/sf_CarpetaCompartida/luca_data/\n"
                       "cab_version: /dev/shm/\n"
                       "\n"
                       "[HomeDirs]\n"
                       "centos_home: /home/user/Work/Projects/pipeline_elvis\n"
                       "cab_home: /pcdisk/holly/sgongora/Dev/Euclid-tests/pipeline_elvis\n"
                       "\n"
                       "[ImagesDirs]\n"
                       "fits_dir: ELViS_v3/CCDs\n"
                       "fpas_dir: ELViS_v3/FPAs\n"
                       "\n"
                       "[ImagesTime]\n"
                       "time_1: 2021-06-26T09:00:00.00000\n"
                       "time_2: 2021-06-26T09:16:43.00000\n"
                       "time_3: 2021-06-26T09:33:26.00000\n"
                       "time_4: 2021-06-26T09:50:09.00000\n"
                       "\n"
                       "[CatsDirs]\n"
                       "output_cats: ELViS_v3/outputs\n"
                       "references: ELViS_v3/cats\n"
                       "filtered: ELViS_v3/filtered\n"
                       "\n"
                       "[CatsOrganization]\n"
                       "first_star: 1\n"
                       "first_galaxy: 12471\n"
                       "first_sso: 134895\n"
                       "\n"
                       "[OutputDirs]\n"
                       "results_dir: /results\n"
                       "catalogs_dir: /results/catalogs\n"
                       "filter_dir: /results/filtered\n"
                       "\n"
                       "[ConfigDirs]\n"
                       "conf_scamp: /configuration_files/scamp.conf\n"
                       "conf_sex: /configuration_files/sextractor.conf\n"
                       "params_sex: /configuration_files/sextractor.params\n"
                       "neural_sex: /configuration_files/sextractor.nnw\n"
                       "params_cat: /configuration_files/sextractor_cat.params\n"
                       "logger_config: /configuration_files/logging.conf\n"
                       "\n"
                       "[Verbosity]\n"
                       "verbose_swarp: QUIET\n"
                       "verbose_sextractor: QUIET\n"
                       "verbose_scamp: QUIET\n"
                       "\n"
                       "[Misc]\n"
                       "detections: 3\n"
                       "pm_low: 0.0005\n"
                       "pm_up: 30\n"
                       "pm_sn: 2\n"
                       "confidences: 5\n"
                       "r_fit: 0.90\n"
                       "cores_number: 0\n"
                       "tolerance: 0.000138889\n")
        testFile.close()

    def create_right_settings_file_cores_define(self):
        """
        @return:
        """
        testFile = open(".settings_ELViS.ini", "w")
        testFile.write("[Version]\n"
                       "centos_version: /media/sf_CarpetaCompartida/luca_data/\n"
                       "cab_version: /dev/shm/\n"
                       "\n"
                       "[HomeDirs]\n"
                       "centos_home: /home/user/Work/Projects/pipeline_elvis\n"
                       "cab_home: /pcdisk/holly/sgongora/Dev/Euclid-tests/pipeline_elvis\n"
                       "\n"
                       "[ImagesDirs]\n"
                       "fits_dir: ELViS_v3/CCDs\n"
                       "fpas_dir: ELViS_v3/FPAs\n"
                       "\n"
                       "[ImagesTime]\n"
                       "time_1: 2021-06-26T09:00:00.00000\n"
                       "time_2: 2021-06-26T09:16:43.00000\n"
                       "time_3: 2021-06-26T09:33:26.00000\n"
                       "time_4: 2021-06-26T09:50:09.00000\n"
                       "\n"
                       "[CatsDirs]\n"
                       "output_cats: ELViS_v3/outputs\n"
                       "references: ELViS_v3/cats\n"
                       "filtered: ELViS_v3/filtered\n"
                       "\n"
                       "[CatsOrganization]\n"
                       "first_star: 1\n"
                       "first_galaxy: 12471\n"
                       "first_sso: 134895\n"
                       "\n"
                       "[OutputDirs]\n"
                       "results_dir: /results\n"
                       "catalogs_dir: /results/catalogs\n"
                       "filter_dir: /results/filtered\n"
                       "\n"
                       "[ConfigDirs]\n"
                       "conf_scamp: /configuration_files/scamp.conf\n"
                       "conf_sex: /configuration_files/sextractor.conf\n"
                       "params_sex: /configuration_files/sextractor.params\n"
                       "neural_sex: /configuration_files/sextractor.nnw\n"
                       "params_cat: /configuration_files/sextractor_cat.params\n"
                       "logger_config: /configuration_files/logging.conf\n"
                       "\n"
                       "[Verbosity]\n"
                       "verbose_swarp: QUIET\n"
                       "verbose_sextractor: QUIET\n"
                       "verbose_scamp: QUIET\n"
                       "\n"
                       "[Misc]\n"
                       "detections: 3\n"
                       "pm_low: 0.0005\n"
                       "pm_up: 30\n"
                       "pm_sn: 2\n"
                       "confidences: 5\n"
                       "r_fit: 0.90\n"
                       "cores_number: 10\n"
                       "tolerance: 0.000138889\n")
        testFile.close()

    def setUp(self):
        """

        :return:
        """
        pass

    def test_right_settings_file_for_centos(self):
        """

        :return:
        """
        misc.get_os = MagicMock(return_value='centos')
        self.create_right_settings_file()

        return self.assertIs(type(misc.extract_settings_elvis()), dict)

    def test_right_settings_file_for_centos_set_cores_number(self):
        """

        :return:
        """
        misc.get_os = MagicMock(return_value='centos')
        self.create_right_settings_file_cores_define()

        prfs_d = misc.extract_settings_elvis()

        return self.assertIs(prfs_d['cores_number'], 10)

    def test_right_settings_file_for_cab(self):
        """

        :return:
        """
        misc.get_os = MagicMock(return_value='cab')
        self.create_right_settings_file()

        return self.assertIs(type(misc.extract_settings_elvis()), dict)

    def test_right_settings_file_for_cab_set_cores_number(self):
        """

        :return:
        """
        misc.get_os = MagicMock(return_value='cab')
        self.create_right_settings_file_cores_define()

        prfs_d = misc.extract_settings_elvis()

        return self.assertIs(prfs_d['cores_number'], 10)

    def test_right_settings_file_wrong_os(self):
        """

        :return:
        """
        misc.get_os = MagicMock(return_value='wrongOS')
        self.create_right_settings_file()

        return self.assertRaises(BadSettings, misc.extract_settings_elvis)

    def tearDown(self):
        """

        :return:
        """
        os.remove(".settings_ELViS.ini")


if __name__ == '__main__':
    main()
