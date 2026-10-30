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

from setuptools import setup, find_packages

# Read long description
with open('README.rst') as f:
    long_description = f.read()


setup(
    name='psamm',
    version='1.2.1',
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
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.13',
        'Programming Language :: Python :: 3.14',
        'Programming Language :: Python :: 3.15',
    ],

    python_requires='>=3.11',

    package_data={'psamm': ['external-data/chebi_pH7_3_mapping.tsv',
                            'external-data/biomass_compound_descriptions.tsv',
                            'external-data/biomass_reaction_descriptions.tsv',
                            'external-data/tcdb_substrates.tsv',
                            'external-data/tcdb_families.tsv',
                            'external-data/biomass_compound_descriptions.tsv',
                            'external-data/biomass_reaction_descriptions.tsv']},

    packages=find_packages(),

    entry_points='''
        [console_scripts]
        psamm-model = psamm.command:main
        psamm-sbml-model = psamm.command:main_sbml
        psamm-list-lpsolvers = psamm.lpsolver.generic:list_solvers
        psamm-import = psamm.importer:main
        psamm-import-bigg = psamm.importer:main_bigg
        psamm-generate-model = psamm.generate_model:main

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

        [psamm.generate_model]
        generate-database = psamm.generate_model:main_databaseCommand
        generate-biomass = psamm.generate_biomass:main
        generate-transporters = psamm.generate_model:main_transporterCommand
    ''',

    test_suite='psamm.tests',

    install_requires=[
        'pyyaml>=4.2b1',
        'six',
        'xlsxwriter',
        'numpy',
        'scipy',
        'future',
        'pandas>=2.0',
    ],
    extras_require={
        'docs': ['sphinx', 'sphinx_rtd_theme', 'mock']
    })
