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
import misc
import platform

__author__ = "Samuel Gongora-Garcia"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Gongora-Garcia"]
__version__ = "0.1"
__maintainer__ = "Samuel Gongora-Garcia"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


class TestGetOs(TestCase):
    """

    """

    def setUp(self):
        """

        :return:
        """
        pass

    def test_fedora_23(self):
        """

        :return:
        """
        platform.platform = MagicMock(return_value='fedora-23')

        self.assertIs(misc.get_os(), 'cab')

    def test_fedora_19(self):
        """

        :return:
        """
        platform.platform = MagicMock(return_value='fedora-19')

        self.assertIs(misc.get_os(), 'cab')

    def test_centos(self):
        """

        :return:
        """
        platform.platform = MagicMock(return_value='centos')

        self.assertIs(misc.get_os(), 'centos')

    def test_exception(self):
        """

        :return:
        """
        platform.platform = MagicMock(return_value='wrongOS')

        self.assertRaises(WrongOS, misc.get_os)

    def tearDown(self):
        """

        :return:
        """
        pass


if __name__ == '__main__':
    main()
