#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Utilities/helpers for testing"""

import os
import sys


class HidePrints:
    """ 
    Context manager that hides all standard output within its body.
    The optional hide_prints argument can be used to quickly enable printing during debugging of tests.


    Use (This will result in all console output of noisy_function to be suppressed):
    from tests.utils import HidePrints
    with HidePrints():
        noise_function()
    """

    @staticmethod
    def hide(func, *args, **kwargs):
        """Decorator that hides all prints of the decorated function.

        Use:
        from tests.utils import HidePrints
        @HidePrints.hide
        def noise()
            print("NOISE")
        """

        def silencer(*args, **kwargs):
            with HidePrints():
                func(*args, **kwargs)
        return silencer

    def __init__(self, hide_prints=True):
        """Create object with hide_prints=False to disable print hiding"""
        self.hide = hide_prints

    def __enter__(self):
        """Redirect stdout to /dev/null, save old stdout"""
        if self.hide:
            self._original_stdout = sys.stdout
            devnull = open(os.devnull, "w")
            sys.stdout = devnull

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Reinstate old stdout"""
        if self.hide:
            sys.stdout = self._original_stdout