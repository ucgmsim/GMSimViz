# -*- coding: utf-8 -*-
"""
Functions used in the Matplotlib plotting scripts.
"""

from collections import OrderedDict
import os
import matplotlib.pyplot as plt


def save_figure(fig, out_dir, basename, png=True, eps=False, close=False):
    """
    "Shortcut" to run fig.savefig
    :param fig: matplotlib figure object
    :param out_dir: directory for image output
    :param basename: basename for images
    :param png: save a png
    :param eps: save an eps
    :param close: close fig when done
    :return:
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


def show_legend(centre=False, extra_labels={}):
    """
    Collapses the legend to remove mulitple entries of the same name
    
    extra_labels: (dict) Optionally adds extra labels to the legend
    centre: optional arguement to place the legend in the top-centre of the plot
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
    Map the float function onto a list of strings.
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

    return map(extract_numbers, string_list)
