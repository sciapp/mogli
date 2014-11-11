# coding: utf-8
"""
Tests for mogli.read()
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import mogli
import pytest


def test_read_xyz_file():
    """
    Read an xyz file and check the number of atoms.
    """
    molecules = mogli.read('examples/dna.xyz')
    assert len(molecules) == 1
    assert len(molecules[0].positions) == 486


def test_read_xyz_file_and_get_atomic_radii():
    """
    Read an xyz file, check the number of atomic_radii and their memoization.
    """
    molecule = mogli.read('examples/dna.xyz')[0]
    atomic_radii = molecule.atomic_radii
    assert len(atomic_radii) == 486
    assert atomic_radii is molecule.atomic_radii


def test_read_unknown_file():
    """
    Read an unknown file and check the expected failure.
    """
    with pytest.raises(mogli.UnknownFileFormatException):
        mogli.read('unknown')[0]
