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

from unittest import TestCase, main
from mock import MagicMock
from pandas import DataFrame

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import misc


__author__ = "Samuel Gongora-Garcia"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Gongora-Garcia"]
__version__ = "0.1"
__maintainer__ = "Samuel Gongora-Garcia"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


class TestCheckSource(TestCase):
    """

    """

    def setUp(self):
        """

        :return:
        """
        dict_ = {'ALPHA_J2000': [0.999861111, 1.000138889],
                 'DELTA_J2000': [0.999861111, 1.000138889]}
        self.df_ = DataFrame(dict_)

    def test_check_in_source(self):
        """
        create_configurations should return:
        [[30, 0.01, 1.5, 4, 'models/gauss_2.0_5x5.conv']], 1

        prfs_d = extract_settings_elvis()
        # How many times the cross search radius will be increased
        tolerance_factor = 2

        df = df[df[keys[0]] + prfs_d['tolerance'] * tolerance_factor > i_alpha]
        df = df[i_alpha > df[keys[0]] - prfs_d['tolerance'] * tolerance_factor]
        df = df[df[keys[1]] + prfs_d['tolerance'] * tolerance_factor > i_delta]
        df = df[i_delta > df[keys[1]] - prfs_d['tolerance'] * tolerance_factor]

        :return:
        """
        misc.extract_settings_elvis = MagicMock()
        misc.extract_settings_elvis.return_value = {'tolerance': 0.000138889}
        i_alpha = 1.0
        i_delta = 1.0
        keys = ['ALPHA_J2000', 'DELTA_J2000']
        test = misc.check_source(self.df_, i_alpha, i_delta, keys)

        return self.assertFalse(test.empty)

    def test_check_out_source(self):
        """
        create_configurations should return:
        [[30, 0.01, 1.5, 4, 'models/gauss_2.0_5x5.conv']], 1

        prfs_d = extract_settings_elvis()
        # How many times the cross search radius will be increased
        tolerance_factor = 2

        df = df[df[keys[0]] + prfs_d['tolerance'] * tolerance_factor > i_alpha]
        df = df[i_alpha > df[keys[0]] - prfs_d['tolerance'] * tolerance_factor]
        df = df[df[keys[1]] + prfs_d['tolerance'] * tolerance_factor > i_delta]
        df = df[i_delta > df[keys[1]] - prfs_d['tolerance'] * tolerance_factor]

        :return:
        """
        misc.extract_settings_elvis = MagicMock()
        misc.extract_settings_elvis.return_value = {'tolerance': 0.000138889}
        i_alpha = 2.0
        i_delta = 2.0
        keys = ['ALPHA_J2000', 'DELTA_J2000']
        test = misc.check_source(self.df_, i_alpha, i_delta, keys)

        return self.assertTrue(test.empty)

    def tearDown(self):
        """

        :return:
        """
        pass


if __name__ == '__main__':
    main()









