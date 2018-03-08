#!/usr/bin/env python2
"""
Required functions for sosfiltfilt (from scipy 0.18+) to work with scipy 0.16+.
Functions copied from scipy git 13/07/2016 9AM NZST.
"""

from numpy import atleast_2d
from scipy.signal import sosfilt, sosfilt_zi
from scipy.signal._arraytools import axis_slice, axis_reverse, odd_ext, even_ext, const_ext

def _validate_pad(padtype, padlen, x, axis, ntaps):
    """Helper to validate padding for filtfilt"""
    if padtype not in ['even', 'odd', 'constant', None]:
        raise ValueError(("Unknown value '%s' given to padtype.  padtype "
                          "must be 'even', 'odd', 'constant', or None.") %
                         padtype)

    if padtype is None:
        padlen = 0

    if padlen is None:
        # Original padding; preserved for backwards compatibility.
        edge = ntaps * 3
    else:
        edge = padlen

    # x's 'axis' dimension must be bigger than edge.
    if x.shape[axis] <= edge:
        raise ValueError("The length of the input vector x must be at least "
                         "padlen, which is %d." % edge)

    if padtype is not None and edge > 0:
        # Make an extension of length `edge` at each
        # end of the input array.
        if padtype == 'even':
            ext = even_ext(x, edge, axis=axis)
        elif padtype == 'odd':
            ext = odd_ext(x, edge, axis=axis)
        else:
            ext = const_ext(x, edge, axis=axis)
    else:
        ext = x
    return edge, ext

def sosfiltfilt(sos, x, axis=-1, padtype='odd', padlen=None):
    """
    A forward-backward filter using cascaded second-order sections.
    See `filtfilt` for more complete information about this method.
    Parameters
    ----------
    sos : array_like
        Array of second-order filter coefficients, must have shape
        ``(n_sections, 6)``. Each row corresponds to a second-order
        section, with the first three columns providing the numerator
        coefficients and the last three providing the denominator
        coefficients.
    x : array_like
        The array of data to be filtered.
    axis : int, optional
        The axis of `x` to which the filter is applied.
        Default is -1.
    padtype : str or None, optional
        Must be 'odd', 'even', 'constant', or None.  This determines the
        type of extension to use for the padded signal to which the filter
        is applied.  If `padtype` is None, no padding is used.  The default
        is 'odd'.
    padlen : int or None, optional
        The number of elements by which to extend `x` at both ends of
        `axis` before applying the filter.  This value must be less than
        ``x.shape[axis] - 1``.  ``padlen=0`` implies no padding.
        The default value is::
            3 * (2 * len(sos) + 1 - min((sos[:, 2] == 0).sum(),
                                        (sos[:, 5] == 0).sum()))
        The extra subtraction at the end attempts to compensate for poles
        and zeros at the origin (e.g. for odd-order filters) to yield
        equivalent estimates of `padlen` to those of `filtfilt` for
        second-order section filters built with `scipy.signal` functions.
    Returns
    -------
    y : ndarray
        The filtered output with the same shape as `x`.
    See Also
    --------
    filtfilt, sosfilt, sosfilt_zi
    Notes
    -----
    .. versionadded:: 0.18.0
    """
    sos, n_sections = _validate_sos(sos)

    # `method` is "pad"...
    ntaps = 2 * n_sections + 1
    ntaps -= min((sos[:, 2] == 0).sum(), (sos[:, 5] == 0).sum())
    edge, ext = _validate_pad(padtype, padlen, x, axis,
                              ntaps=ntaps)

    # These steps follow the same form as filtfilt with modifications
    zi = sosfilt_zi(sos)  # shape (n_sections, 2) --> (n_sections, ..., 2, ...)
    zi_shape = [1] * x.ndim
    zi_shape[axis] = 2
    zi.shape = [n_sections] + zi_shape
    x_0 = axis_slice(ext, stop=1, axis=axis)
    (y, zf) = sosfilt(sos, ext, axis=axis, zi=zi * x_0)
    y_0 = axis_slice(y, start=-1, axis=axis)
    (y, zf) = sosfilt(sos, axis_reverse(y, axis=axis), axis=axis, zi=zi * y_0)
    y = axis_reverse(y, axis=axis)
    if edge > 0:
        y = axis_slice(y, start=edge, stop=-edge, axis=axis)
    return y

def _validate_sos(sos):
    """Helper to validate a SOS input"""
    sos = atleast_2d(sos)
    if sos.ndim != 2:
        raise ValueError('sos array must be 2D')
    n_sections, m = sos.shape
    if m != 6:
        raise ValueError('sos array must be shape (n_sections, 6)')
    if not (sos[:, 3] == 1).all():
        raise ValueError('sos[:, 3] should be all ones')
    return sos, n_sections
