[![Build Status](https://travis-ci.org/sgongar/Euclid-SSOs-Pipeline-Tests.svg?branch=master)](https://travis-ci.org/sgongar/Euclid-SSOs-Pipeline-Tests) [![Coverage Status](https://coveralls.io/repos/github/sgongar/Euclid-SSOs-Pipeline-Tests/badge.svg?branch=master)](https://coveralls.io/github/sgongar/Euclid-SSOs-Pipeline-Tests?branch=master) [![Code Health](https://landscape.io/github/sgongar/Euclid-SSOs-Pipeline-Tests/master/landscape.svg?style=flat)](https://landscape.io/github/sgongar/Euclid-SSOs-Pipeline-Tests/master)

### Performance module
Utilities created to check pipeline behaviour. Using these scripts pipeline performance can be validated through a collection of output catalogues and regions. These catalogues show us the objects detected (or not), organised by source type. For each catalogue there is a region file intended to be use in DS9.

##### To-do
- [ ] Implement not detected catalogue for SSOs.
- [ ] Implement not detected catalogue for stars.
- [ ] Implement not detected catalogue for galaxies.
- [x] Implement detected catalogue for SSOs.
- [ ] Increment coverage to ninety percent.

##### Extracted data
All the following scripts were created to get the data obtained by sextractor from real sources. These catalogues, stored in CSV format, can be easily read by any office suite.

###### extracted_catalog_elvis_galaxies.py
Creates a catalogue populated by galaxies.

###### extracted_catalog_elvis_stars.py
Creates a catalogue populated by stars.

###### extracted_catalog_elvis_stars.py
Creates a catalogue populated by SSOs.

###### extracted_regions_elvis_galaxies.py
Creates a regions file populated by galaxies.

###### extracted_regions_elvis_stars.py
Creates a regions file populated by stars.

###### extracted_regions_elvis_stars.py
Creates a regions file populated by SSOs.
