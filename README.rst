PSAMM metabolic modeling tools
==============================

.. image:: https://travis-ci.org/zhanglab/psamm.svg?branch=develop
   :alt: Build Status
   :target: https://travis-ci.org/zhanglab/psamm
.. image:: https://readthedocs.org/projects/psamm/badge/?version=latest
   :alt: Documentation Status
   :target: https://readthedocs.org/projects/psamm/?badge=latest
.. image:: https://badge.fury.io/py/psamm.svg
   :alt: Package Status
   :target: https://pypi.python.org/pypi/psamm

PSAMM is an open source software that is designed for the curation and analysis
of metabolic models. It supports model version tracking, model annotation, data
integration, data parsing and formatting, consistency checking, automatic gap
filling, and model simulations.

See NEWS_ for information on recent changes. The ``master`` branch
tracks the latest release while the ``develop`` branch is the latest version in
development. Please apply any pull requests to the ``develop`` branch when
creating the pull request.

.. _NEWS: NEWS.md

Install
-------

Use ``pip`` to install (it is recommended to use a Virtualenv_):

.. code:: shell

    $ pip install psamm

The ``psamm-import`` tool is developed in `a separate repository`_. After
installing PSAMM the ``psamm-import`` tool can be installed using:

.. code:: shell

    $ pip install git+https://github.com/zhanglab/psamm-import.git

.. _Virtualenv: https://virtualenv.pypa.io/
.. _a separate repository: https://github.com/zhanglab/psamm-import

Documentation
-------------

The documentation for PSAMM is available at `Read the Docs`_.

.. _Read the Docs: https://psamm.readthedocs.org/

Software license
----------------

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

See LICENSE_.

.. _LICENSE: LICENSE
