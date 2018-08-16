#!/usr/bin/python
# -*- coding: utf-8 -*-

""" Utilities for fits images management.

Todo:
    * Improve log messages
    * Improve usability
"""
from os import listdir

from astropy.io import fits
from astropy.wcs import WCS

from misc import extract_settings_elvis, extract_settings_luca


__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def get_fits_d(mag_, dither):
    """ Returns a list with all fits files in the fits directory of a chosen
    dither and magnitude bin.

    :param mag_: The chosen magnitude bin.
    :param dither: The chosen dither.
    :return: A list with the fits files.
    """
    prfs_d = extract_settings_luca()
    fits_list = []

    files = listdir('{}/{}/CCDs/'.format(prfs_d['fits_dir'], mag_))
    for file_ in files:
        if file_[:1] == 'm' and file_[-5:] == '.fits':
            fits_list.append(file_)

    list_out = []
    for file_ in fits_list:
        if file_[-6:-5] == str(dither):
            list_out.append(file_)

    return list_out


def get_fits(dither):
    """  Returns a list with all fits files in the fits directory of a chosen
    dither.

    :param dither: The chosen dither.
    :return: The list with the fits files.
    """
    prfs_d = extract_settings_elvis()  # Gets the preferences dictionary.
    fits_list = []  # Creates an empty list of fits files.

    # Gets all fits files.
    files = listdir('{}/'.format(prfs_d['fits_dir']))
    for file_ in files:
        if file_[-5:] == '.fits':
            fits_list.append(file_)

    # Gets the images of the chosen dither.
    fits_out = []
    for file_ in fits_list:
        if file_[-6:-5] == str(dither):  # The
            fits_out.append(file_)

    return fits_out


def get_fits_limits(fits_image):
    """

    :param fits_image: The chosen image.
    :return: a dict with ra/dec limits above_ra, below_ra,
    above_dec_, below_dec
    """

    data, header = fits.getdata(fits_image, header=True)
    w = WCS(fits_image)

    above_x, above_y = header['NAXIS1'], header['NAXIS2']
    above_ra, above_dec = w.all_pix2world(above_x, above_y, 0)

    below_ra, below_dec = w.all_pix2world(0, 0, 0)

    # Useless, it is really important?
    # c = SkyCoord(ra=[above_ra, below_ra] * degree,
    #              dec=[above_dec, below_dec] * degree)

    ra = [above_ra, below_ra]
    dec = [above_dec, below_dec]

    limits = {'below_ra': float(min(ra)), 'above_ra': float(max(ra)),
              'below_dec': float(min(dec)), 'above_dec': float(max(dec))}

    return limits



