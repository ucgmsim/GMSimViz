#!/usr/bin/env python

from glob import glob
from multiprocessing import cpu_count
import os
from shutil import rmtree
from subprocess import call, Popen, PIPE

tests_dir = os.path.abspath(os.path.dirname(__file__))
sample_data = os.path.abspath(os.path.join(tests_dir, os.pardir, "sample_data"))
ref_dir = os.path.join(tests_dir, "reference_frames")
out_dir = os.path.join(tests_dir, "test_frames")
nproc = cpu_count() - 1


def img_diff(img_a, img_b):
    cmd = ["magick", "compare", "-metric", "AE", img_a, img_b, "NULL:"]
    p = Popen(cmd, stderr=PIPE)
    d = p.communicate()[1]
    # proportion of resolution
    try:
        return (float(d) * 100.0) / (1536 * 864)
    except ValueError:
        # could be file does not exist
        return 100.0


if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

call(
    [
        "gmsimviz",
        os.path.join(sample_data, "fault.srf"),
        "-a",
        "--gm-cut",
        "10.0",
        "--crude",
        "-n{}".format(nproc),
        "--title",
        "Test Animation",
        "--dpi",
        "96",
        "-f",
        "10",
        "--downscale",
        "1",
        "-x",
        os.path.join(sample_data, "xyts.e3d"),
        "--liquefaction-s",
        os.path.join(sample_data, "liquefaction_s.hdf5"),
        "--liquefaction-p",
        os.path.join(sample_data, "liquefaction_p.hdf5"),
        "--landslide-s",
        os.path.join(sample_data, "landslide_s.hdf5"),
        "--landslide-p",
        os.path.join(sample_data, "landslide_p.hdf5"),
        "--paths",
        os.path.join(sample_data, "transport"),
        "--temp",
        out_dir,
        "-k",
    ]
)

ref_imgs = glob(os.path.join(ref_dir, "*.png"))
diff_total = 0.0
for ref_img in ref_imgs:
    d = img_diff(ref_img, os.path.join(out_dir, os.path.basename(ref_img)))
    print("{:6.2f}% {}".format(d, os.path.basename(ref_img)))
    diff_total += d

print("average pixels changed: {:.2f}%".format(diff_total / len(ref_imgs)))
rmtree(out_dir)
