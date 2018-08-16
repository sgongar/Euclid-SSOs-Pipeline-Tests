#!/usr/bin/python
# -*- coding: utf-8 -*-

"""

Versions:
- 0.1: Initial release. Split from ssos_catalog_creation.py

Information:
- cat: -> hdu_list catalogue
- data: -> Table formatted data
- df: -> dataframe formatted data

Todo:
    *

*GNU Terry Pratchett*

"""
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.units import degree
from astropy.wcs import WCS

from misc import extract_settings_elvis
from misc_fits import get_fits

__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def get_borders():
    """

    :return:
    """
    prfs_d = extract_settings_elvis()
    borders_d = {}
    for dither_ in range(1, 5, 1):
        borders_d[dither_] = {}
        fits_list = get_fits(dither_)
        for fits_ in fits_list:
            limits = get_fits_limits('{}/{}'.format(prfs_d['fits_dir'], fits_))

            borders_d[dither_][fits_] = limits

    return borders_d


def get_fits_limits(fits_image):
    """ todo - to another new file?

    @param fits_image: fits image

    @return limits: a dict with ra/dec limits above_ra, below_ra,
                    above_dec_, below_dec
    """
    # logger.info('getting limits of {} image'.format(fits_image))

    data, header = fits.getdata(fits_image, header=True)
    w = WCS(fits_image)

    above_x, above_y = header['NAXIS1'], header['NAXIS2']
    above_ra, above_dec = w.all_pix2world(above_x, above_y, 0)

    below_ra, below_dec = w.all_pix2world(0, 0, 0)

    above_coords = SkyCoord(ra=above_ra * degree,
                            dec=above_dec * degree,
                            equinox='J2021.5')
    below_coords = SkyCoord(ra=below_ra * degree,
                            dec=below_dec * degree,
                            equinox='J2021.5')

    above_coords_converted = above_coords.geocentrictrueecliptic
    below_coords_converted = below_coords.geocentrictrueecliptic

    lon = [above_coords_converted.lon.degree,
           below_coords_converted.lon.degree]
    lat = [above_coords_converted.lat.degree,
           below_coords_converted.lat.degree]

    limits = {'below_lon': float(min(lon)), 'above_lon': float(max(lon)),
              'below_lat': float(min(lat)), 'above_lat': float(max(lat))}

    # check position
    # sometimes some values could be higher when are tagged as "lowest"
    return limits
