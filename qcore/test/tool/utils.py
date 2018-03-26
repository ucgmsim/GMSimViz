import numpy as np
import filecmp
import os


ERROR_LIMIT = 0.001

def compare_np_array(array1, array2, error_limit=ERROR_LIMIT):
    """array1: a numpy array from sample output, will be used as the denominator,
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
    try:
        filecmp.cmp(file1, file2)
    except (OSError,IOError) as e:
        print (e)


def remove_file(abs_path):
    try:
        os.remove(abs_path)
    except (IOError, OSError):
        raise