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
from project_warnings import TooFastSource

import misc
from misc import get_norm_speed


__author__ = "Samuel Gongora-Garcia"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Gongora-Garcia"]
__version__ = "0.1"
__maintainer__ = "Samuel Gongora-Garcia"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


class TestGetNormSpeed(TestCase):
    """

    """

    def setUp(self):
        """

        :return:
        """
        pass

    def test_right_speed(self):
        """

        :return:
        """
        right_speed = 0.01

        return self.assertEqual(misc.get_norm_speed(right_speed), 0.01)

    def test_zero_speed(self):
        """

        :return:
        """
        zero_speed = 0.001

        return self.assertEqual(misc.get_norm_speed(zero_speed), 0)

    def test_too_fast_source(self):
        """

        :return:
        """
        zero_speed = 100

        return self.assertRaises(TooFastSource, get_norm_speed, zero_speed)

    def tearDown(self):
        """

        :return:
        """
        pass


if __name__ == '__main__':
    main()
