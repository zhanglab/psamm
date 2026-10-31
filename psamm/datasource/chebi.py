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
# Copyright 2026  Jason Vailionis <jason_vailionis@uri.edu>

"""Module related to loading ChEBI compound data.

The charge and formula of ChEBI compounds are read from the ChEBI flat files
(``chemical_data.tsv.gz`` and ``secondary_ids.tsv.gz``) which are downloaded
from the EBI FTP site and cached locally.
"""

import csv
import gzip
import logging
import os
import shutil
import tempfile
import time
from urllib.request import urlopen

logger = logging.getLogger(__name__)

CHEBI_FLAT_FILE_URL = 'https://ftp.ebi.ac.uk/pub/databases/chebi/flat_files/'
CHEMICAL_DATA_FILE = 'chemical_data.tsv.gz'
SECONDARY_IDS_FILE = 'secondary_ids.tsv.gz'
CHEBI_FILES = (CHEMICAL_DATA_FILE, SECONDARY_IDS_FILE)


def default_cache_dir():
    """Return the default directory for caching the ChEBI flat files.

    The directory can be set with the environment variable ``PSAMM_CHEBI_DIR``
    and otherwise defaults to ``psamm/chebi`` in the user cache directory.
    """
    if 'PSAMM_CHEBI_DIR' in os.environ:
        return os.environ['PSAMM_CHEBI_DIR']
    cache_home = os.environ.get('XDG_CACHE_HOME') or os.path.join(
        os.path.expanduser('~'), '.cache')
    return os.path.join(cache_home, 'psamm', 'chebi')


def download_files(directory, update=False, url=None):
    """Download the ChEBI flat files to directory.

    Files that already exist are only downloaded again if update is True.
    The files are downloaded from url (default is CHEBI_FLAT_FILE_URL).
    """
    if url is None:
        url = CHEBI_FLAT_FILE_URL
    os.makedirs(directory, exist_ok=True)
    for filename in CHEBI_FILES:
        path = os.path.join(directory, filename)
        if os.path.exists(path) and not update:
            continue
        logger.info('Downloading {}{}'.format(url, filename))
        # Download to a temporary file first so that an interrupted download
        # does not leave a partial file in the cache.
        fd, tmp_path = tempfile.mkstemp(dir=directory, suffix='.part')
        try:
            with os.fdopen(fd, 'wb') as f, urlopen(url + filename) as r:
                shutil.copyfileobj(r, f)
            # mkstemp creates the file readable only by the owner; use the
            # normal permissions so that the cache can be shared.
            umask = os.umask(0)
            os.umask(umask)
            os.chmod(tmp_path, 0o666 & ~umask)
            os.replace(tmp_path, path)
        except BaseException:
            os.remove(tmp_path)
            raise


def load(directory=None, update=False):
    """Return a :class:`ChebiDb` from the cached ChEBI flat files.

    The files are downloaded if they are not already present in directory
    (default is :func:`default_cache_dir`) or if update is True.
    """
    if directory is None:
        directory = default_cache_dir()
    download_files(directory, update=update)
    chemical_data = os.path.join(directory, CHEMICAL_DATA_FILE)
    logger.info('Using ChEBI data from {} (downloaded {})'.format(
        directory, time.strftime(
            '%Y-%m-%d', time.localtime(os.path.getmtime(chemical_data)))))
    return ChebiDb(chemical_data,
                   os.path.join(directory, SECONDARY_IDS_FILE))


def _normalize_id(chebi_id):
    """Return ChEBI ID as a string without the ``CHEBI:`` prefix."""
    chebi_id = str(chebi_id).strip()
    if chebi_id.upper().startswith('CHEBI:'):
        chebi_id = chebi_id[len('CHEBI:'):]
    return chebi_id


def _read_tsv(path):
    with gzip.open(path, 'rt', encoding='utf-8', newline='') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            yield row


class ChebiDb(object):
    """Charge and formula of ChEBI compounds.

    Loaded from the ChEBI flat files ``chemical_data.tsv.gz`` and (optionally)
    ``secondary_ids.tsv.gz``. Secondary IDs are resolved to the primary ID.
    IDs can be given with or without the ``CHEBI:`` prefix.
    """

    def __init__(self, chemical_data_path, secondary_ids_path=None):
        self._chemical_data = {}
        for row in _read_tsv(chemical_data_path):
            self._chemical_data[row['compound_id']] = (
                row['formula'] or None, row['charge'] or None)

        self._primary_ids = {}
        if secondary_ids_path is not None:
            for row in _read_tsv(secondary_ids_path):
                self._primary_ids[row['secondary_id']] = row['compound_id']

    def __contains__(self, chebi_id):
        return self.primary_id(chebi_id) in self._chemical_data

    def __len__(self):
        return len(self._chemical_data)

    def primary_id(self, chebi_id):
        """Return the primary ChEBI ID (without prefix) of an ID."""
        chebi_id = _normalize_id(chebi_id)
        return self._primary_ids.get(chebi_id, chebi_id)

    def get_formula(self, chebi_id):
        """Return formula string of compound or None if unknown."""
        return self._chemical_data.get(
            self.primary_id(chebi_id), (None, None))[0]

    def get_charge(self, chebi_id):
        """Return charge of compound as an int or None if unknown."""
        charge = self._chemical_data.get(
            self.primary_id(chebi_id), (None, None))[1]
        return int(charge) if charge is not None else None
