# -*- coding: utf-8 -*-
"""
Functions used in the Matplotlib plotting scripts.
"""

from collections import OrderedDict
import os

import matplotlib.pyplot as plt

def save_figure(fig, out_dir, basename, png = True, eps = False, \
        close = False, return_fig = False):
    """
    "Shortcut" to run fig.savefig, optionally returns fig.gca().
    fig: matplotlib figure object
    out_dir: directory for image output
    basename: basename for images
    png: save a png
    eps: save an eps
    close: close fig when done
    """
    # path is given in components
    path = os.path.join(out_dir, basename)

    # make sure destination dir exists
    if not os.path.exists(out_dir):
        try:
            os.makedirs(out_dir)
        except OSError:
            if not os.path.exists(out_dir):
                raise

    if png:
        fig.savefig('%s.png' % (path))
    if eps:
        fig.savefig('%s.eps' % (path))

    if close:
        plt.close(fig)

    # XXX: VERY VERY BAD, LEGACY, DEPRECATED.
    if return_fig:
        return fig, fig.gca()

def show_legend(centre=False, extra_labels={}):
    """
    TODO: description
    """
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    by_label.update(extra_labels)

    if centre:
        plt.legend(by_label.values(), by_label.keys(), loc='upper center', \
                bbox_to_anchor=(0.5, 1.05), mode='expand', \
                ncol=3, fancybox=True, shadow=True)
    else:
        plt.legend(by_label.values(), by_label.keys(), loc='best', fontsize=9)

def convert_strings_to_floats(string_list):
    """
    Questionable logic to map the float function onto a list of strings.
    List of comma-separated floats for each string are returned
    BUT if there was one value then the type returned is the single float.
    """

    def extract_numbers(number_string):
        """
        Returns all floats in comma seperated number string.
        """
        floats = map(float, number_string.split(','))
        if len(floats) == 1:
            return floats[0]
        return floats

    return map(extract_number, string_list)

def is_virtual_station(station_name):
    """
    Checks if all restraints on virtual station names are met.
    """
    # 7 characters long
    if len(station_name) != 7:
        return False

    # no capitals
    if sum(map(str.isupper, tuple(station_name))):
        return False

    # valid hex string
    try:
        if not isinstance(s, int):
            int(station_name, 16)
    except (ValueError, TypeError):
        return False

    # all tests passed
    return True
