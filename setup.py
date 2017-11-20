from setuptools import setup, find_packages

pkgs = find_packages()
print pkgs
setup(
    name="qcore",
    version="1.0",
    packages=pkgs,
    url="https://github.com/ucgmsim/qcore",
    description="QuakeCoRE Library",
    package_data={'': ['*.json', '*.sh', 'cpt/*.cpt']},
    install_requires=['numpy', 'mpi4py', 'scipy'],
)
