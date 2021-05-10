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
# Copyright 2014-2017  Jon Lund Steffensen <jon_steffensen@uri.edu>
# Copyright 2015-2020  Keith Dufault-Thompson <keitht547@my.uri.edu>

from __future__ import print_function

import sys
from setuptools import setup, find_packages
import pkg_resources

# Read long description
with open('README.rst') as f:
    long_description = f.read()

# Test whether psamm-import is currently installed. Since the psamm-import
# functionality was moved to this package (except Excel importers), only newer
# versions of psamm-import are compatible with recent versions of PSAMM.
try:
    pkg_resources.get_distribution('psamm-import <= 0.15.2')
except (pkg_resources.DistributionNotFound,
        pkg_resources.VersionConflict):
    pass
else:
    msg = (
        'Please upgrade or uninstall psamm-import before upgrading psamm:\n'
        '$ pip install --upgrade psamm-import\n'
        ' OR\n'
        '$ pip uninstall psamm-import'
        '\n\n'
        ' The functionality of the psamm-import package has been moved into'
        ' the psamm package, and the psamm-import package now only contains'
        ' the model-specific Excel importers.')
    print(msg, file=sys.stderr)
    sys.exit(1)


setup(
    name='psamm',
    version='1.1.2',
    description='PSAMM metabolic modeling tools',
    maintainer='Jon Lund Steffensen',
    maintainer_email='jon_steffensen@uri.edu',
    url='https://github.com/zhanglab/psamm',
    license='GNU GPLv3+',

    long_description=long_description,

    classifiers=[
        'Development Status :: 4 - Beta',
        (
            'License :: OSI Approved :: '
            'GNU General Public License v3 or later (GPLv3+)'),
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],

    packages=find_packages(),

    entry_points='''
        [console_scripts]
        psamm-model = psamm.command:main
        psamm-sbml-model = psamm.command:main_sbml
        psamm-list-lpsolvers = psamm.lpsolver.generic:list_solvers
        psamm-import = psamm.importer:main
        psamm-import-bigg = psamm.importer:main_bigg

        [psamm.commands]
        chargecheck = psamm.commands.chargecheck:ChargeBalanceCommand
        console = psamm.commands.console:ConsoleCommand
        dupcheck = psamm.commands.duplicatescheck:DuplicatesCheck
        excelexport = psamm.commands.excelexport:ExcelExportCommand
        fastgapfill = psamm.commands.fastgapfill:FastGapFillCommand
        fba = psamm.commands.fba:FluxBalanceCommand
        fluxcheck = psamm.commands.fluxcheck:FluxConsistencyCommand
        fluxcoupling = psamm.commands.fluxcoupling:FluxCouplingCommand
        formulacheck = psamm.commands.formulacheck:FormulaBalanceCommand
        fva = psamm.commands.fva:FluxVariabilityCommand
        gapcheck = psamm.commands.gapcheck:GapCheckCommand
        gapfill = psamm.commands.gapfill:GapFillCommand
        genedelete = psamm.commands.genedelete:GeneDeletionCommand
        gimme = psamm.commands.gimme:GimmeCommand
        masscheck = psamm.commands.masscheck:MassConsistencyCommand
        primarypairs = psamm.commands.primarypairs:PrimaryPairsCommand
        randomsparse = psamm.commands.randomsparse:RandomSparseNetworkCommand
        robustness = psamm.commands.robustness:RobustnessCommand
        sbmlexport = psamm.commands.sbmlexport:SBMLExport
        search = psamm.commands.search:SearchCommand
        tableexport = psamm.commands.tableexport:ExportTableCommand
        psammotate = psamm.commands.psammotate:PsammotateCommand
        modelmapping = psamm.commands.model_mapping:ModelMappingCommand
        vis = psamm.commands.vis:VisualizationCommand
        tmfa = psamm.commands.tmfa:TMFACommand

        [psamm.importer]
        JSON = psamm.importers.cobrajson:Importer
        SBML = psamm.importers.sbml:NonstrictImporter
        SBML-strict = psamm.importers.sbml:StrictImporter
        MATLAB = psamm.importers.matlab:Importer
    ''',

    test_suite='psamm.tests',

    install_requires=[
        'pyyaml>=4.2b1',
        'six',
        'xlsxwriter',
        'numpy',
        'scipy',
        'future',
        'pandas'
    ],
    extras_require={
        'docs': ['sphinx', 'sphinx_rtd_theme', 'mock'],
        ':python_version=="2.7"': ['enum34'],
        ':python_version=="3.3"': ['enum34']
    })
