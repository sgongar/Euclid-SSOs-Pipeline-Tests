#!/usr/bin/python
# -*- coding: utf-8 -*-

"""

Content:
    * get_os
    * conf_map
    * extract_settings_luca

Todo:
    * Improve log messages
    * Improve usability
"""
from math import hypot

from multiprocessing import cpu_count
from ConfigParser import ConfigParser
import platform
from logging import getLogger, config

from errors import BadSettings, InvalidScampConfiguration, WrongOS
from project_warnings import TooFastSource

__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def get_norm_mag(o_mag):
    """

    :param o_mag:
    :return: mag_bin
    """
    mags = [[11, 12], [12, 13], [13, 14], [14, 15], [15, 16], [16, 17],
            [17, 18], [18, 19], [19, 20], [20, 21], [21, 22], [22, 23],
            [23, 24], [24, 25], [25, 26], [26, 27], [27, 28], [28, 29]]

    mag_bin = ''
    for mag_ in mags:
        if mag_[0] < o_mag < mag_[1]:
            mag_bin = '{}-{}'.format(mag_[0], mag_[1])

    return mag_bin


def get_norm_speed(o_pm):
    """

    :return:
    """
    pm_norm = 0
    speeds_d = {0.01: [0.005, 0.015], 0.03: [0.015, 0.05],
                0.1: [0.05, 0.15], 0.3: [0.15, 0.5],
                1.0: [0.5, 1.5], 3.0: [1.5, 5],
                10.0: [5.0, 15.0], 30.0: [15.0, 50]}

    if o_pm < 0.005:
        pm_norm = 0
    elif o_pm > 50.0:
        raise TooFastSource
    else:
        # pm_norm = 0
        for key_ in speeds_d.keys():
            low = speeds_d[key_][0]
            high = speeds_d[key_][1]
            if low < o_pm <= high:
                pm_norm = key_

    return pm_norm


def create_sextractor_dict(conf_num, cat_conf):
    """

    @param conf_num:
    @param cat_conf:

    @return analysis_d:
    """

    if cat_conf:
        configurations = [2, 0.1, 5, 4, 'models/gauss_2.0_5x5.conv']
        len_conf = 1
    else:
        mode = {'type': 'sextractor'}  # hard-coded
        configurations, len_conf = create_configurations(mode)

    analysis_l = []

    if type(configurations[0]) is list:
        for configuration in range(len(configurations)):
            temp_list = [configurations[configuration][0],
                         configurations[configuration][1],
                         configurations[configuration][2],
                         configurations[configuration][2],
                         configurations[configuration][3],
                         configurations[configuration][4]]
            analysis_l.append(temp_list)
        analysis_d = {'deblend_nthresh': analysis_l[conf_num][0],
                      'deblend_mincount': analysis_l[conf_num][1],
                      'detect_thresh': analysis_l[conf_num][2],
                      'analysis_thresh': analysis_l[conf_num][3],
                      'detect_minarea': analysis_l[conf_num][4],
                      'filter': analysis_l[conf_num][5]}
    else:
        analysis_d = {'deblend_nthresh': configurations[0],
                      'deblend_mincount': configurations[1],
                      'detect_thresh': configurations[2],
                      'analysis_thresh': configurations[2],
                      'detect_minarea': configurations[3],
                      'filter': configurations[4]}

    return analysis_d, len_conf


def create_configurations(mode):
    """ creates a list of configuration lists merging
        different input parameters
    :param mode: can be 'sextractor' or 'scamp'
    :return:
    """
    if mode['type'] == 'sextractor':
        l_deblending = [30]
        l_mincount = [0.01]
        l_threshold = [1.5]

        l_area = [4]
        l_filter_name = ['models/gauss_2.0_5x5.conv']

        configurations = []
        for deblending in l_deblending:
            for mincount in l_mincount:
                for threshold in l_threshold:
                    for area in l_area:
                        for filt in l_filter_name:
                            configurations.append([deblending, mincount,
                                                   threshold, area, filt])
        configurations_len = len(configurations)
        return configurations, configurations_len

    elif mode['type'] == 'scamp':
        l_crossid_radius = [10]  # [10] seconds
        l_pixscale_maxerr = [1.1]  # [1.2] scale-factor
        l_posangle_maxerr = [0.5]  # [0.5, 2.5] degrees
        l_position_maxerr = [0.04]

        configurations = []

        for crossid in l_crossid_radius:
            for pixscale in l_pixscale_maxerr:
                for posangle in l_posangle_maxerr:
                    for position in l_position_maxerr:
                        configurations.append([crossid, pixscale,
                                               posangle, position])

        configurations_len = len(configurations)

        return configurations, configurations_len


def create_scamp_dict(conf_num):
    """

    :param conf_num:
    :return:
    """

    scamp_list = []
    mode = {'type': 'scamp'}
    configurations, len_conf = create_configurations(mode)

    for conf in configurations:
        temp_list = [conf[0], conf[1], conf[2], conf[3]]
        scamp_list.append(temp_list)

    try:
        scamp_dict = {'crossid_radius': scamp_list[conf_num][0],
                      'pixscale_maxerr': scamp_list[conf_num][1],
                      'posangle_maxerr': scamp_list[conf_num][2],
                      'position_maxerr': scamp_list[conf_num][3]}
    except IndexError:
        raise InvalidScampConfiguration

    return scamp_dict, len_conf


def setting_logger(prfs_d, logger_name):
    """ sets-up a logger object ready to be used

    TODO improve logger definition

    @return logger:
    """
    config.fileConfig(prfs_d['logger_config'])

    # TODO implement logger level setting
    """
    if argv[4] == '-INFO':
        logger = getLogger("main_process").setLevel(INFO)
    elif argv[4] == '-DEBUG':
        logger = getLogger("main_process").setLevel(DEBUG)
    else:
        raise Exception
    """
    logger = getLogger(logger_name)
    logger.info("Pipeline started")
    # logger.FileHandler('spam.log')

    return logger


def get_os():
    """ a function that gets the current operative system
    for now works in Debian, Fedora and Ubuntu shell (Microsoft version)

    @return os_system: a string which contains the operative system name
    """

    if 'fedora-23' in platform.platform():
        os_system = 'test'
    elif 'debian' in platform.platform():
        os_system = 'debian'
    elif 'Ubuntu' in platform.platform():
        os_system = 'ubuntu'
    elif 'fedora-26' in platform.platform():
        os_system = 'fedora'
    elif 'fedora-19' in platform.platform():
        os_system = 'cab'
    elif 'centos' in platform.platform():
        os_system = 'centos'
    else:
        raise WrongOS

    return os_system


def conf_map(config_, section):
    """

    @param config_:
    @param section:

    @return dict1:
    """
    dict1 = {}
    options = config_.options(section)
    for option in options:
        try:
            dict1[option] = config_.get(section, option)
            if dict1[option] == -1:
                print('skip: {}'.format(option))
        except KeyError:
            print('exception on {}'.format(option))
            dict1[option] = None
    return dict1


def extract_settings_elvis():
    """ Creates a dictionary with all the configuration parameters
        at this moment configuration file location is fixed at main directory

    @return prfs_d: a dictionary which contains all valuable data.
    """
    cf = ConfigParser()
    cf.read(".settings_ELViS.ini")  # Location of the configuration parameters

    prfs_d = {}
    os_version = get_os()  # Gets the current operating system.

    if os_version == 'centos':
        prfs_d['version'] = conf_map(cf, "Version")['centos_version']
        prfs_d['home'] = conf_map(cf, "HomeDirs")['centos_home']
    elif os_version == 'cab':
        prfs_d['version'] = conf_map(cf, "Version")['cab_version']
        prfs_d['home'] = conf_map(cf, "HomeDirs")['cab_home']
    elif os_version == 'ubuntu':
        prfs_d['version'] = conf_map(cf, "Version")['ubuntu_version']
        prfs_d['home'] = conf_map(cf, "HomeDirs")['ubuntu_home']
    elif os_version == 'debian':
        prfs_d['version'] = conf_map(cf, "Version")['debian_version']
        prfs_d['home'] = conf_map(cf, "HomeDirs")['debian_home']
    else:
        raise BadSettings('Operative system not chosen')

    prfs_d['fits_dir'] = conf_map(cf, "ImagesDirs")['fits_dir']
    prfs_d['fits_dir'] = '{}{}'.format(prfs_d['version'], prfs_d['fits_dir'])
    prfs_d['fpas_dir'] = conf_map(cf, "ImagesDirs")['fpas_dir']
    prfs_d['fpas_dir'] = '{}{}'.format(prfs_d['version'], prfs_d['fpas_dir'])

    # todo - comment!
    prfs_d['output_cats'] = conf_map(cf, "CatsDirs")['output_cats']
    prfs_d['output_cats'] = prfs_d['version'] + prfs_d['output_cats']
    # todo - comment!
    prfs_d['references'] = conf_map(cf, "CatsDirs")['references']
    prfs_d['references'] = prfs_d['version'] + prfs_d['references']
    # todo - comment!
    prfs_d['filtered'] = conf_map(cf, "CatsDirs")['filtered']
    prfs_d['filtered'] = prfs_d['version'] + prfs_d['filtered']

    prfs_d['time_1'] = conf_map(cf, "ImagesTime")['time_1']  # 1st dither time
    prfs_d['time_2'] = conf_map(cf, "ImagesTime")['time_2']  # 2nd dither time
    prfs_d['time_3'] = conf_map(cf, "ImagesTime")['time_3']  # 3nd dither time
    prfs_d['time_4'] = conf_map(cf, "ImagesTime")['time_4']  # 4th dither time

    outputdirs_list = ['conf_scamp', 'conf_sex', 'params_sex', 'neural_sex',
                       'params_cat', 'logger_config']
    for conf_ in outputdirs_list:
        prfs_d[conf_] = conf_map(cf, "ConfigDirs")[conf_]
        prfs_d[conf_] = prfs_d['home'] + prfs_d[conf_]

    prfs_d['detections'] = int(conf_map(cf, "Misc")['detections'])
    prfs_d['pm_low'] = float(conf_map(cf, "Misc")['pm_low'])
    prfs_d['pm_up'] = float(conf_map(cf, "Misc")['pm_up'])
    prfs_d['pm_sn'] = float(conf_map(cf, "Misc")['pm_sn'])
    pms = conf_map(cf, "Misc")['pms']
    pms = pms.replace(",", " ")
    prfs_d['pms'] = [float(x) for x in pms.split()]
    prfs_d['r_fit'] = conf_map(cf, "Misc")['r_fit']
    prfs_d['cores_number'] = conf_map(cf, "Misc")['cores_number']
    if prfs_d['cores_number'] == '0':
        prfs_d['cores_number'] = int(str(cpu_count()))
        # TODO should leave free at least 20% of processors
    else:
        prfs_d['cores_number'] = int(prfs_d['cores_number'])
    prfs_d['tolerance'] = float(conf_map(cf, "Misc")['tolerance'])

    return prfs_d


def check_source(df, i_alpha, i_delta, keys):
    """ Checks if there is any source with a right ascension and declination
    coordinates in a given dataframe

    :param df: The given dataframe.
    :param i_alpha: Right ascension of the source to be located.
    :param i_delta: Declination of the source to be located.
    :param keys: A list with the names of the columns of right ascension and
    declination.
    :return: A dataframe called df populated with the found source, if there
    is no source the dataframe will be empty.
    """
    prfs_d = extract_settings_elvis()
    # How many times the cross search radius will be increased
    tolerance_factor = 2

    df = df[df[keys[0]] + prfs_d['tolerance']*tolerance_factor > i_alpha]
    df = df[i_alpha > df[keys[0]] - prfs_d['tolerance']*tolerance_factor]
    df = df[df[keys[1]] + prfs_d['tolerance']*tolerance_factor > i_delta]
    df = df[i_delta > df[keys[1]] - prfs_d['tolerance']*tolerance_factor]

    return df


def check_distance(o_df, alpha, delta):
    """

    :param o_df:
    :param alpha:
    :param delta:
    :return:
    """
    distance_l = []
    for ix, row in o_df.iterrows():
        distance = hypot(row.ALPHA_J2000 - alpha, row.DELTA_J2000 - delta)
        distance_l.append(distance)

    index = distance_l.index(min(distance_l))

    return index
