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

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import misc


__author__ = "Samuel Gongora-Garcia"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Gongora-Garcia"]
__version__ = "0.1"
__maintainer__ = "Samuel Gongora-Garcia"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


class TestCreateConfigurations(TestCase):
    """

    """

    def setUp(self):
        """

        :return:
        """
        pass

    def test_sextractor(self):
        """
        create_configurations should return:
        [[30, 0.01, 1.5, 4, 'models/gauss_2.0_5x5.conv']], 1

        :return:
        """
        mode = {'type': 'sextractor'}
        configurations = misc.create_configurations(mode)

        statement_1 = self.assertIs(type(configurations[0]), list)
        statement_2 = self.assertIs(type(configurations[0][0]), list)
        statement_3 = self.assertIs(configurations[1], 1)

        return statement_1 and statement_2 and statement_3

    def test_scamp(self):
        """

        :return:
        """
        mode = {'type': 'scamp'}
        configurations = misc.create_configurations(mode)

        statement_1 = self.assertIs(type(configurations[0]), list)
        statement_2 = self.assertIs(type(configurations[0][0]), list)
        statement_3 = self.assertIs(configurations[1], 1)

        return statement_1 and statement_2 and statement_3

    def tearDown(self):
        """

        :return:
        """
        pass


if __name__ == '__main__':
    main()
