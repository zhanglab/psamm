#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name='psamm',
    version='0.9',
    description='PSAMM metabolic modeling tools',
    maintainer='Jon Lund Steffensen',
    maintainer_email='jon_steffensen@uri.edu',
    url='https://github.com/zhanglab/model_script',
    license='GNU GPLv3+',

    classifiers = [
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)'
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
