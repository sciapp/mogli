# coding: utf-8
"""
Tests for mogli.export()

These tests only check if the function does not fail, they do not check the
results.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os.path
import mogli


def test_export():
    """
    Export a molecule to an HTML file.
    """
    molecule = mogli.read('examples/dna.xyz')[0]
    mogli.export(molecule, 'test.html')
    assert os.path.exists('test.html')
