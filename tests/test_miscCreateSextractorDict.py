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


class TestCreateSextractorDict(TestCase):
    """

    """

    def setUp(self):
        """

        :return:
        """
        pass

    def test_configuration_for_catalogue(self):
        """

        :return:
        """
        sex_dict = misc.create_sextractor_dict(1, True)

        statement_1 = self.assertIs(type(sex_dict[0]), dict)
        statement_2 = self.assertIs(sex_dict[1], 1)

        return statement_1 and statement_2

    def test_configuration_for_images(self):
        """

        :return:
        """
        sex_dict = misc.create_sextractor_dict(0, False)

        statement_1 = self.assertIs(type(sex_dict[0]), dict)
        statement_2 = self.assertIs(sex_dict[1], 1)

        return statement_1 and statement_2

    def tearDown(self):
        """

        :return:
        """
        pass


if __name__ == '__main__':
    main()
