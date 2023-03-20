# coding: utf-8
"""
Tests for mogli.show()

These tests only check if the function does not fail, they do not check the
results.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import signal
import mogli
import glfw


def test_show():
    """
    Show a molecule and close the window after 1 second.
    """

    def alarm_handler(signal_number, stack_frame):
        """
        Close the current GLFW window when called due to a SIGALRM.
        """
        assert signal_number == signal.SIGALRM
        window = glfw.get_current_context()
        assert window is not None
        glfw.set_window_should_close(window, True)

    molecule = mogli.read("examples/dna.xyz")[0]
    # Set the signal handler
    signal.signal(signal.SIGALRM, alarm_handler)
    # Set the alarm to 1 second (and check no previous alarm had been set)
    assert signal.alarm(1) == 0
    # Show the molecule (and open the GLFW window closed by the alarm handler)
    mogli.show(molecule)


def test_show_custom_colors():
    """
    Show a molecule and close the window after 1 second.
    """

    def alarm_handler(signal_number, stack_frame):
        """
        Close the current GLFW window when called due to a SIGALRM.
        """
        assert signal_number == signal.SIGALRM
        window = glfw.get_current_context()
        assert window is not None
        glfw.set_window_should_close(window, True)

    molecule: mogli.Molecule = mogli.read("examples/dna.xyz")[0]

    # Set a list of custom colors
    molecule.atom_colors = [
        (i / (molecule.atom_count - 1), 0, 0) for i in range(molecule.atom_count)
    ]

    # Set the signal handler
    signal.signal(signal.SIGALRM, alarm_handler)
    # Set the alarm to 1 second (and check no previous alarm had been set)
    assert signal.alarm(1) == 0
    # Show the molecule (and open the GLFW window closed by the alarm handler)
    mogli.show(molecule)


if __name__ == "__main__":
    test_show()
