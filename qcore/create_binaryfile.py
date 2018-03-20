""" Creates 2 benchmark binary files for testing of the functions srf.srf2llv(srffile) and srf.srf2llv_py(srffile)
First create these binary files and put in the sample output folder, (if you are testing it for a new SRF) """

from qcore import srf
import numpy as np


SRFFILE = "/home/aas105/qcore/qcore/test/test_srf/sample2/input/Tuakana13_HYP01-01_S1244.srf"


def create_file(srffile, islist=False):
    """srffile: abs path to input srf file
       islist: srffile function return a numpy array or a list(of numpy array(s))
       islist is true for srf2llv_py; false for srfllv. default false
    """
    filename = "out_array_srf2llv.bin"
    if not islist:
        out_array = srf.srf2llv(srffile)
    else:
        parts = filename.split('.')
        filename = parts[0] + '_py.' + parts[1]
        out_array_list = srf.srf2llv_py(srffile)
        out_array = out_array_list[0]
        for array in out_array_list[1:]:
                out_array = np.concatenate([out_array, array])
                print("first out array", out_array)
    try:
        out_array.astype(np.float32).tofile(filename)
    except Exception as (e):
        print e


create_file(SRFFILE)
create_file(SRFFILE, islist=True)
