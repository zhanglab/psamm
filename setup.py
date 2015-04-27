#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name='psamm',
    version='0.7',
    description='PSAMM metabolic modelrm -ing tools',
    maintainer='Jon Lund Steffensen',
    maintainer_email='jon_steffensen@uri.edu',
    url='https://github.com/zhanglab/model_script',

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
