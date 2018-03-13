from setuptools import setup, find_packages

setup(
    name='qcore',
    version='1.1',
    packages=['qcore'],
    url='https://github.com/ucgmsim/qcore',
    description='QuakeCoRE Library',
    package_data={'qcore': ['*.json']},
    install_requires=['numpy', 'scipy>=0.16'],
)
