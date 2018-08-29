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

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from errors import WrongOS
import catalog_creation_elvis

__author__ = "Samuel Gongora-Garcia"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Gongora-Garcia"]
__version__ = "0.1"
__maintainer__ = "Samuel Gongora-Garcia"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


class TestCatalogCreationElvisFunctions(TestCase):
    """

    """

    def setUp(self):
        """

        :return:
        """
        pass

    def test_create_empty_catalog_keys(self):
        """

        :return:
        """
        test = catalog_creation_elvis.create_empty_catalog_dict()

        statement_1 = self.assertIs(len(test.keys()), 9)
        statement_2 = self.assertIs(type(test), dict)

        return statement_1 and statement_2

    def tearDown(self):
        """

        :return:
        """
        pass


if __name__ == '__main__':
    main()
