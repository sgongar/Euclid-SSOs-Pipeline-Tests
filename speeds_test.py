#!/usr/bin/python
# -*- coding: utf-8 -*-

"""

Todo:
    * Improve log messages

"""

from collections import Counter
from filecmp import cmp
from multiprocessing import Process
from timeit import default_timer as timer

from astropy.io import fits
from astropy.table import Table
from numpy import array, sqrt
from pandas import concat, Series

__author__ = "Samuel Gongora-Garcia"
__copyright__ = "Copyright 2017"
__credits__ = ["Samuel Gongora-Garcia"]
"""
__license__ = "GPL"
"""
__version__ = "0.1"
__maintainer__ = "Samuel Gongora-Garcia"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def method_a(full_db, pm, pme, pmalpha, pmdelta, pmealpha, pmedelta):
    """ mide la velocidad de creacion del nuevo dataframe utilizando for loops

    @return full_db:
    """
    step_5_1 = timer()

    full_db = full_db.loc[~full_db['CATALOG_NUMBER'].isin([0])]

    # Filtering by detections
    full_db = concat(g for _, g in full_db.groupby("SOURCE_NUMBER")
                     if len(g) >= 3)

    for idx in set(full_db['SOURCE_NUMBER']):
        full_db.loc[full_db['SOURCE_NUMBER'] == idx,
                    'PM'] = pm.loc[idx - 1]
        full_db.loc[full_db['SOURCE_NUMBER'] == idx,
                    'PMERR'] = pme.loc[idx - 1]
        full_db.loc[full_db['SOURCE_NUMBER'] == idx,
                    'PMALPHA'] = pmalpha.loc[idx - 1]
        full_db.loc[full_db['SOURCE_NUMBER'] == idx,
                    'PMDELTA'] = pmdelta.loc[idx - 1]
        full_db.loc[full_db['SOURCE_NUMBER'] == idx,
                    'PMALPHAERR'] = pmealpha.loc[idx - 1]
        full_db.loc[full_db['SOURCE_NUMBER'] == idx,
                    'PMDELTAERR'] = pmedelta.loc[idx - 1]

    step_5_2 = timer()

    print 'elapsed time for a {}'.format(step_5_2 - step_5_1)

    full_db.to_csv('test_a.csv')


def method_b(full_db, pm, pme, pmalpha, pmdelta, pmealpha, pmedelta):
    """ mide la creacion del nuevo script creando listas a partir de las 
        frecuencias de cada valor

    @return full_db: a new Dataframe, named as the old one, populated with 
    the proper motion values
    """
    step_5_1 = timer()

    full_db = full_db.loc[~full_db['CATALOG_NUMBER'].isin([0])]

    freq_values = []

    freq_dict = Counter(full_db['SOURCE_NUMBER'].tolist())
    sources_list = freq_dict.keys()

    for source_ in sources_list:
        freq_values.append(freq_dict[source_])

    pm_list = []
    pme_list = []
    pmalpha_list = []
    pmdelta_list = []
    pmealpha_list = []
    pmedelta_list = []

    for idx, freq in enumerate(freq_values):
        tmp = [pm.iloc[idx]] * freq
        pm_list.append(tmp)
        tmp = [pme.iloc[idx]] * freq
        pme_list.append(tmp)
        tmp = [pmalpha.iloc[idx]] * freq
        pmalpha_list.append(tmp)
        tmp = [pmdelta.iloc[idx]] * freq
        pmdelta_list.append(tmp)
        tmp = [pmealpha.iloc[idx]] * freq
        pmealpha_list.append(tmp)
        tmp = [pmedelta.iloc[idx]] * freq
        pmedelta_list.append(tmp)

    pm_list = [item for sublist in pm_list for item in sublist]
    pme_list = [item for sublist in pme_list for item in sublist]
    pmalpha_list = [item for sublist in pmalpha_list for item in sublist]
    pmdelta_list = [item for sublist in pmdelta_list for item in sublist]
    pmealpha_list = [item for sublist in pmealpha_list for item in sublist]
    pmedelta_list = [item for sublist in pmedelta_list for item in sublist]

    full_db['PM'] = Series(data=pm_list, index=full_db.index)
    full_db['PMERR'] = Series(data=pme_list, index=full_db.index)
    full_db['PMALPHA'] = Series(data=pmalpha_list, index=full_db.index)
    full_db['PMDELTA'] = Series(data=pmdelta_list, index=full_db.index)
    full_db['PMALPHAERR'] = Series(data=pmealpha_list, index=full_db.index)
    full_db['PMDELTAERR'] = Series(data=pmedelta_list, index=full_db.index)

    full_db.to_csv('test_b.csv')

    step_5_2 = timer()

    print 'elapsed time for b {}'.format(step_5_2 - step_5_1)

    return full_db


def main():
    print 'loading full catalog...'
    step_1 = timer()

    fll_n = 'results/full_120_1.2_0.5_0.083_20-21_1.cat'
    full_cat = fits.open(fll_n)
    full_db = Table(full_cat[2].data)
    full_db = full_db.to_pandas()

    print 'loading merged catalog...'
    step_2 = timer()
    print 'elapsed time {}'.format(step_2 - step_1)

    mrgd_n = 'results/merged_120_1.2_0.5_0.083_20-21_1.cat'
    merged_cat = fits.open(mrgd_n)
    merged_db = Table(merged_cat[2].data)

    print 'getting "/h from miliarc/year'
    step_3 = timer()
    print 'elapsed time {}'.format(step_3 - step_2)

    pmalpha = Series(merged_db.field('PMALPHA_J2000') / 8.75e6)  # 8.75e6
    pmdelta = Series(merged_db.field('PMDELTA_J2000') / 8.75e6)
    pmealpha = Series(merged_db.field('PMALPHAERR_J2000') / 8.75e6)
    pmedelta = Series(merged_db.field('PMDELTAERR_J2000') / 8.75e6)

    print 'getting proper motions'
    step_4 = timer()
    print 'elapsed time {}'.format(step_4 - step_3)

    pm = Series(sqrt(array(pmalpha**2 + pmdelta**2), dtype=float))
    pme = Series(sqrt(array(pmealpha**2 + pmedelta**2), dtype=float))

    print 'creating new DataFrame'
    step_5 = timer()
    print 'elapsed time {}'.format(step_5 - step_4)

    method_j = []
    method_a_p = Process(target=method_a,
                         args=(full_db, pm, pme, pmalpha, pmdelta,
                               pmealpha, pmedelta,))
    method_j.append(method_a_p)
    method_a_p.start()

    method_b_p = Process(target=method_b,
                         args=(full_db, pm, pme, pmalpha, pmdelta,
                               pmealpha, pmedelta,))
    method_j.append(method_b_p)
    method_b_p.start()

    active_method = list([job.is_alive() for job in method_j])
    while True in active_method:
        active_method = list([job.is_alive() for job in method_j])
        pass

    step_6 = timer()
    print 'elapsed time for a and b {}'.format(step_6 - step_5)


if __name__ == '__main__':
    main()

    if not cmp('test_a.csv', 'test_b.csv'):
        raise Exception
    else:
        print 'test passed'
