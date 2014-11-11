# coding: utf-8
"""
Tests for mogli.draw()

These tests only check if the function does not fail, they do not check the
results.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import gr
import mogli


def test_draw():
    """
    Draw a molecule using GR
    """
    molecule = mogli.read('examples/dna.xyz')[0]
    # Draw the molecule (with the current GR workstation type)
    mogli.draw(molecule)
    # Close GR again
    gr.emergencyclosegks()
