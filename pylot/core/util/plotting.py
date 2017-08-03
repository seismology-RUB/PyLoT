#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt


def create_bin_list(l_boundary, u_boundary, nbins=100):
    """
    takes two boundaries and a number of bins and creates a list of bins for
    histogram plotting
    :param l_boundary: Any number.
    :type  l_boundary: float
    :param u_boundary: Any number that is greater than l_boundary.
    :type  u_boundary: float
    :param nbins: Any positive integer.
    :type  nbins: int
    :return: A list of equidistant bins.
    """
    if u_boundary <= l_boundary:
        raise ValueError('Upper boundary must be greather than lower!')
    elif nbins <= 0:
        raise ValueError('Number of bins is not valid.')
    binlist = []
    for i in range(nbins):
        binlist.append(l_boundary + i * (u_boundary - l_boundary) / nbins)
    return binlist


def histplot(array, binlist, xlab='Values',
             ylab='Frequency', title=None, fnout=None):
    """
    function to quickly show some distribution of data. Takes array like data,
    and a list of bins. Editing detail and inserting a legend is not possible.
    :param array: List of values.
    :type  array: Array like
    :param binlist: List of bins.
    :type  binlist: list
    :param xlab: A label for the x-axes.
    :type  xlab: str
    :param ylab: A label for the y-axes.
    :type  ylab: str
    :param title: A title for the Plot.
    :type  title: str
    :param fnout: A path to save the plot instead of showing.
                  Has to contain filename and type. Like: 'path/to/file.png'
    :type  fnout. str
    :return: -
    """

    plt.hist(array, bins=binlist)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    if title:
        plt.title(title)
    if fnout:
        plt.savefig(fnout)
    else:
        plt.show()
