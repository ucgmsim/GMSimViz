""" Creates 2 benchmark binary files for testing of the functions srf.srf2llv(srffile) and srf.srf2llv_py(srffile)
First create these binary files and put in the sample output folder, (if you are testing it for a new SRF) 
command to run:
python create_binaryfile.py test/test_srf/sample2/input/Tuakana13_HYP01-01_S1244.srf """

from qcore import srf
import numpy as np
import sys

# output for Python and C implementations
OUT_C = "out_array_srf2llv.bin"
OUT_PY = "out_array_srf2llv_py.bin"

def srf2bin(srf_file, out_file, py_method=False):
    """srffile: path to input srf file
       islist: srffile function return a numpy array or a list(of numpy array(s))
       islist is true for srf2llv_py; false for srfllv. default false
    """

    # load data
    if py_method:
        out_array_list = srf.srf2llv_py(srf_file)
        out_array = out_array_list[0]
        for array in out_array_list[1:]:
            out_array = np.concatenate([out_array, array])
    else:
        out_array = srf.srf2llv(srf_file)

    # save data
    try:
        out_array.astype(np.float32).tofile(out_file)
    except Exception as e:
        sys.exit(e)

# retrieve input srf from command line
try:
    srf_file = sys.argv[1]
except IndexError:
    sys.exit("Please provide the path to the SRF file as an argument.")

# produce 2 versions using the different methods
srf2bin(srf_file, OUT_C)
srf2bin(srf_file, OUT_PY, py_method=True)
