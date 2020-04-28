#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions to deal with console logging and debugging.

@author: tangir
"""

import os
import numpy as np
import inspect
from logging import DEBUG, ERROR, INFO, WARN
import logging

import pdb

# init the logger
logging.basicConfig(format='(%(levelname)s) %(message)s', level=DEBUG)
# characters used in log messages
_INDENT_CHAR = '>'
_INDENT_CHAR2 = ' '
# separator lines in logs
_LINE_CHAR = '-'
_LINE_WIDTH = 40


def getLevel():
    """
    Return the current logging level: DEBUG, ERROR, INFO or WARN.

    Returns
    -------
    logging.getLogger().level() : int/enum
        Logging level: DEBUG, ERROR, INFO or WARN
    """
    return(logging.getLogger().level)


def setLevel(l):
    """
    Set the current logging level: DEBUG, ERROR, INFO or WARN.

    Parameters
    ----------
    l : int/enum
        Logging level: DEBUG, ERROR, INFO or WARN
    """
    logging.getLogger().setLevel(l)
    # to avoid DEBUG messages from matplotlib!
    logging.getLogger('matplotlib.font_manager').disabled = True


def _wrap_message(s, cl):
    """
    Take the logging string and add the function name, some symbols and indentation to illustrate the current stack.

    Parameters
    ----------
    s : string
        Logging string
    cl : int/enum
        Current logging level: DEBUG, ERROR, INFO or WARN

    Returns
    -------
    s : string
        Logging string with
    """
    # find class' and method's name
    stack = inspect.stack()
    if("self" in list(stack[2][0].f_locals.keys())):
        the_base = stack[2][0].f_locals["self"].__class__.__name__
    else:
        the_base = os.path.basename(stack[2][0].f_code.co_filename)

    the_method = stack[2][0].f_code.co_name

    # final pretty logging string
    s = the_base + "." + the_method + ": " + s
    # add 1 more space if info
    if(cl == INFO):
        s = ' ' + s

    return(s)


def info(s):
    """
    Call logger for some INFO output.

    Parameters
    ----------
    s : string
        Logging string
    """
    s = _wrap_message(s, INFO)
    logging.info(s)


def info_line________________________():
    """Log a line as INFO."""
    logging.info(_LINE_CHAR * _LINE_WIDTH)


def info_line_break():
    """Log a line break as INFO."""
    logging.info("")


def debug(s):
    """
    Call logger for some DEBUG output.

    Parameters
    ----------
    s : string
        Logging string
    """
    s = _wrap_message(s, DEBUG)
    logging.debug(s)


def debug_line________________________():
    """Log a line as DEBUG."""
    logging.debug(_LINE_CHAR * _LINE_WIDTH)


def debug_line_break():
    """Log a line break as DEBUG."""
    logging.info("")


def warning(s):
    """
    Call logger for some WARNING output. Raise a regular python warning.

    Parameters
    ----------
    s : string
        Logging string
    """
    s = _wrap_message(s, WARN)
    logging.warning(s)
    Warning(s)


def error(s):
    """
    Call logger for some ERROR output. Raise a regular python exception.

    Parameters
    ----------
    s : string
        Logging string
    """
    s = _wrap_message(s, ERROR)
    logging.error(s)
    raise Exception(s)


class progressbar:
    """A class that deals with console progressbars to show the progression of a process."""

    def __init__(self, s, number_steps, backspace_mode=False):
        """
        Initialize a progressbar object.

        Parameters
        ----------
        s : string
            Logging string
        number_steps: int
            Size of progress bar / step value for which the bar should be full
        backspace_mode : boolean
            Refresh the whole progress bar (True) or just add characters (False)
        """
        # current bar progression
        self._progress_bar_counter = 0
        # max bar progression
        self._progress_bar_counter_max = number_steps
        # size of bar (number of characters in console)
        self._progressbar_length = 30
        # empty character
        self._empty_char = '░'
        # full character
        self._full_char = '█'
        # print mode
        self._backspace_mode = backspace_mode

        # wrap logging string
        bar_str = s + ": "
        if(self._backspace_mode):
            # add empty an empty progressbar
            bar_str += self._empty_char * self._progressbar_length
        # only print progress bar if we are logging at INFO level or more
        if(logging.getLogger().getEffectiveLevel() <= INFO):
            print(bar_str, end="", flush=True)

    def update(self, current_step):
        """
        Print a small and pretty ASCII progress bar in the terminal. There is maybe some n / n + 1 bug but I don't care.

        Parameters
        ----------
        current_step : int
            Integer describing current status or number of achieved steps
        """
        current_bar_counter = int(np.ceil(current_step / self._progress_bar_counter_max * self._progressbar_length))
        # growing bar?
        if(current_bar_counter > self._progress_bar_counter):
            if(self._backspace_mode):
                # remove backspace whole bar
                bar_str = '\b' * self._progressbar_length
                # fill bar to current progression
                bar_str += self._full_char * current_bar_counter
                # rest of bar is empty
                bar_str += self._empty_char * (self._progressbar_length - current_bar_counter)
            else:
                bar_str = self._full_char * (current_bar_counter - self._progress_bar_counter)
            # only print progress bar if we are logging at INFO level or more
            if(logging.getLogger().getEffectiveLevel() <= INFO):
                print(bar_str, end="", flush=True)

            # remember progression
            self._progress_bar_counter = current_bar_counter

    def finish(self, s):
        """
        Force fill of progressbar with a message.

        Parameters
        ----------
        s : string
            Finish string
        """
        if(self._backspace_mode):
            # remove backspace whole bar
            bar_str = '\b' * self._progressbar_length
            # fill bar to current progression
            bar_str += self._full_char * self._progressbar_length
        else:
            bar_str = ""
        # add finish message
        bar_str += " " + s + "."
        # only print progress bar if we are logging at INFO level or more
        if(logging.getLogger().getEffectiveLevel() <= INFO):
            print(bar_str, flush=True)
