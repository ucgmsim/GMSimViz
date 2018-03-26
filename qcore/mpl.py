# -*- coding: utf-8 -*-
"""
A collection of common functions used in the plotting scripts.

"""

from collections import OrderedDict
import os

import matplotlib.pyplot as plt

def save_figure(fig, filename, directory = None, png = True, eps = False):
    """Function to save .png and .eps and print info to stdout.       """
    if not os.path.exists(directory):
        os.makedirs(directory)
    path = os.path.abspath(os.path.join(directory, filename))
    if png:
        png_filename = path + '.png'
        fig.savefig(png_filename)
        print png_filename + ' saved to disk'
    if eps:
        eps_filename = path + '.eps'
        fig.savefig(eps_filename)
        print eps_filename + ' saved to disk'

def save_and_end_figure(fig, directory, filename, return_figure_handle=False):
    save_figure(fig, directory, filename)
    if return_figure_handle:
        ax = fig.gca()
        return fig, ax
    else:
        plt.close(fig)
        return None

def show_legend(centre=False, extra_labels={}):
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    by_label.update(extra_labels)
    if centre:
        plt.legend(by_label.values(), by_label.keys(), loc='upper center', bbox_to_anchor=(0.5, 1.05), mode='expand',
                    ncol=3, fancybox=True, shadow=True)
    else:
        plt.legend(by_label.values(), by_label.keys(), loc='best', fontsize=9)

def extract_number(input):
    result = input.split(',')
    if len(result) <= 1:
        result = float(result[0])
    else:
        result = map(float, result)
    return result

def convert_string_to_floats(in_text):
    return [CommonPlot.extract_number(ratio) for ratio in in_text]

def is_hex(s):
    try:
        int(s, 16)
        return True
    except ValueError:
        return False

def is_non_uniform(name):
    n_caps = sum(1 for c in name if c.isupper())
    return len(name) == 7 and n_caps == 0 and CommonPlot.is_hex(name)
