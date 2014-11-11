# coding: utf-8
"""
Tests for mogli.Molecule.calculate_bonds()
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import mogli
import pytest


def test_calculate_bonds_with_default_method():
    """
    Calculate bonds with no parameters.
    """
    molecule = mogli.read('examples/dna.xyz')[0]
    assert molecule.bonds is None
    molecule.calculate_bonds()
    assert molecule.bonds is not None
    assert molecule.bonds.method == 'radii'


def test_calculate_bonds_with_wrong_method():
    """
    Calculate bonds with an invalid method name and check the expected failure.
    """
    molecule = mogli.read('examples/dna.xyz')[0]
    with pytest.raises(ValueError):
        molecule.calculate_bonds(method='invalid_method_name')


def test_calculate_bonds_with_radii_method():
    """
    Calculate bonds with the radii method.
    """
    molecule = mogli.read('examples/dna.xyz')[0]
    molecule.calculate_bonds(method='radii')
    assert molecule.bonds is not None


def test_calculate_bonds_with_constant_delta_method_with_default_param():
    """
    Calculate bonds with the constant_delta method without providing a value
    for param and check the expected failure.
    """
    molecule = mogli.read('examples/dna.xyz')[0]
    with pytest.raises(ValueError):
        molecule.calculate_bonds(method='constant_delta')
    assert molecule.bonds is None


def test_calculate_bonds_with_constant_delta_method_with_explicit_param():
    """
    Calculate bonds with the constant_delta method.
    """
    molecule = mogli.read('examples/dna.xyz')[0]
    molecule.calculate_bonds(method='constant_delta', param=3.0)
    assert molecule.bonds is not None


def test_calculate_bonds_empty():
    """
    Calculate bonds in a way that should not produce any bonds.
    """
    molecule = mogli.read('examples/dna.xyz')[0]
    molecule.calculate_bonds(method='constant_delta', param=0.0)
    assert molecule.bonds is not None
    assert len(molecule.bonds) == 0


def test_calculate_bonds_memoization_with_other_method_or_param():
    """
    Calculate bonds and check the function memoization with varying parameters.
    """
    molecule = mogli.read('examples/dna.xyz')[0]
    assert molecule.bonds is None
    molecule.calculate_bonds(method='constant_delta', param=0.0)
    bonds = molecule.bonds
    molecule.calculate_bonds(method='radii', param=0.0)
    assert molecule.bonds is not bonds
    bonds = molecule.bonds
    molecule.calculate_bonds(method='radii')
    assert molecule.bonds is not bonds
    bonds = molecule.bonds
    molecule.calculate_bonds(method='radii', param=0.0)
    assert molecule.bonds is not bonds


def test_calculate_bonds_memoization_with_default_param():
    """
    Calculate bonds and check the function memoization for default paremeters.
    """
    molecule = mogli.read('examples/dna.xyz')[0]
    assert molecule.bonds is None
    molecule.calculate_bonds()
    bonds = molecule.bonds
    molecule.calculate_bonds()
    assert molecule.bonds is bonds


def test_calculate_bonds_memoization_with_explicit_param():
    """
    Calculate bonds and check the function memoization for given paremeters.
    """
    molecule = mogli.read('examples/dna.xyz')[0]
    assert molecule.bonds is None
    molecule.calculate_bonds(param=1.0)
    bonds = molecule.bonds
    molecule.calculate_bonds(param=1.0)
    assert molecule.bonds is bonds
