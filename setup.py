from setuptools import setup, find_packages
import os, sys

error_msg = """Please download and extract the resources folder into GMSimViz.
It is available at: https://goo.gl/mYFCQn"""

script_dir = os.path.abspath(os.path.dirname(__file__))
resources = os.path.join(script_dir, "resources")
if not os.path.isdir(resources):
    sys.exit(error_msg)
with open(os.path.join(script_dir, "gmsimviz", "config.json"), "w") as cfg:
    cfg.write('{\n"GMT_DATA" : "%s"\n}\n' % (resources))

setup(
    name="GMSimViz",
    version="1.1",
    packages=["gmsimviz"],
    url="https://github.com/ucgmsim/GMSimViz",
    description="Ground Motion Visualisation",
    package_data={"gmsimviz": ["*.json", "*.png"]},
    install_requires=["numpy", "scipy>=0.16"],
    scripts=["bin/gmsimviz"],
)
