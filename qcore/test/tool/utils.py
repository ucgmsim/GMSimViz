import numpy as np
import filecmp
import os


ERROR_LIMIT = 0.001

def compare_np_array(array1, array2, error_limit=ERROR_LIMIT):
    """first test if the two numpy arrays are of the same shape
       if pass, then test if the relative error between the two arrays are <= the preset error limit
       if any of the test fails, assertion error will be raised;
       array1: a numpy array from sample output, will be used as the denominator,
       array2: a numpy array from test output, makes part of the numerator.
       error_limit: preset error_limit to be compared with the relative error (array1-array2)/array1
    """
    # print("array1",array1.shape)
    # print("array2",array2.shape)
    # print("diff",array1 - array2)

    assert array1.shape == array2.shape
    relative_error = np.divide((array1 - array2), array1)
    max_relative_error = np.nanmax(np.abs(relative_error))
    assert max_relative_error <= error_limit


def compare_binary_files(file1, file2):
    """ Compare if the two binary files are equal
        if not eqaul, raise assertion error
        file1: abs path to file1
        file2: abs path to file2
    """
    try:
        filecmp.cmp(file1, file2)
    except (OSError,IOError) as e:
        raise AssertionError


def remove_file(abs_path):
    """ remove test output file if successfully passed the test
        abs_path: abs path of the file to be removed
    """
    try:
        os.remove(abs_path)
    except (IOError, OSError):
        raise AssertionError