from setuptools import setup, find_packages
import os
import sys
import tarfile
from time import sleep
from urllib.request import urlretrieve

PACKAGE_NAME = "gmsimviz"
PACKAGE_URL = "https://github.com/ucgmsim/GMSimViz"
DATA_VERSION = "1.4"
DATA_NAME = "GMSimViz_resources.tar.xz"
DATA_URL = f"{PACKAGE_URL}/releases/download/{DATA_VERSION}/{DATA_NAME}"


def extract_data(archive, destination):
    with tarfile.open(archive) as xz:
        xz.extractall(destination)


def get_version(version_path):
    if os.path.isfile(version_path):
        return open(version_path).read().strip()
    else:
        return None


def prepare_data():
    loc_version = os.path.join(PACKAGE_NAME, "data", "version")
    # extract existing archive
    have_ver = get_version(loc_version)
    if str(have_ver) != DATA_VERSION and os.path.isfile(DATA_NAME):
        # extract available archive
        print("checking available archive version...")
        extract_data(DATA_NAME, PACKAGE_NAME)
    # download missing archive
    have_ver = get_version(loc_version)
    print(have_ver)
    if str(have_ver) != DATA_VERSION:
        print("data package missing or incorrect version")
        print("downloading...")
        urlretrieve(DATA_URL, DATA_NAME)
        extract_data(DATA_NAME, PACKAGE_NAME)
    # final check
    have_ver = get_version(loc_version)
    if str(have_ver) != DATA_VERSION:
        sys.exit("data package issue, please contact repository maintainer")


prepare_data()
setup(
    name="GMSimViz",
    version="1.4",
    packages=[PACKAGE_NAME],
    url=PACKAGE_URL,
    description="Ground Motion Visualisation",
    package_data={PACKAGE_NAME: ["data/*", "data/*/*", "data/*/*/*"]},
    include_package_data=True,
    install_requires=["numpy", "mpi4py", "h5py"],
    scripts=["bin/gmsimviz"],
)
