#!/usr/bin/env python

from setuptools import setup, find_packages

setup(name='model_script',
      version='0.3',
      description='Metabolic modelling tools',
      maintainer='Jon Lund Steffensen',
      maintainer_email='jon_steffensen@uri.edu',
      url='https://github.com/zhanglab/model_script',
      packages=find_packages(),
      scripts=['model.py'],
      test_suite='metnet.tests',
      install_requires=['PyYAML>=3.11,<4.0'])
