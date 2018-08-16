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
from errors import InvalidScampConfiguration
import misc

__author__ = "Samuel Gongora-Garcia"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Gongora-Garcia"]
__version__ = "0.1"
__maintainer__ = "Samuel Gongora-Garcia"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


class TestCreateScampDict(TestCase):
    """

    """

    def setUp(self):
        """

        :return:
        """
        pass

    def test_right_index_scamp_configuration(self):
        """

        :return:
        """
        scamp_dict = misc.create_scamp_dict(0)

        statement_1 = self.assertIs(type(scamp_dict[0]), dict)
        statement_2 = self.assertIs(scamp_dict[1], 1)

        return statement_1 and statement_2

    def test_configuration_for_catalogue(self):
        """

        :return:
        """
        self.assertRaises(InvalidScampConfiguration,
                          misc.create_scamp_dict, 1)

    def tearDown(self):
        """

        :return:
        """
        pass


if __name__ == '__main__':
    main()
