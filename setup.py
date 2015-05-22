#!/usr/bin/env python
# This file is part of PSAMM.
#
# PSAMM is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PSAMM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PSAMM.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 2014-2015  Jon Lund Steffensen <jon_steffensen@uri.edu>

from setuptools import setup, find_packages

# Read long description
with open('README.rst') as f:
    long_description = f.read()

setup(
    name='psamm',
    version='0.10.2',
    description='PSAMM metabolic modeling tools',
    maintainer='Jon Lund Steffensen',
    maintainer_email='jon_steffensen@uri.edu',
    url='https://github.com/zhanglab/psamm',
    license='GNU GPLv3+',

    long_description=long_description,

    classifiers = [
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Programming Language :: Python :: 2.7'
    ],

    packages=find_packages(),
    entry_points = {
        'console_scripts': [
            'psamm-model = psamm.command:main'
        ]
    },

    test_suite='psamm.tests',

    install_requires=['PyYAML>=3.11,<4.0'],
    extras_require={
        'docs': ['sphinx', 'mock']
    })
