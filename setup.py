#!/usr/bin/env python
"""Setup script for turbomoleio."""
import sys
import os

from glob import glob
from setuptools import find_packages

# changing Python version requirements.
if sys.version[0:3] < '3.6':
    sys.stderr.write("ERROR: turbomoleio requires Python Version 3.6 or above. "
                     "You are using version {}. Exiting.".format(sys.version))
    sys.exit(1)


def file_doesnt_end_with(test, endings):
    """
    A little utility we'll need below, since glob() does NOT allow you to do exclusion on multiple endings!

    Return true if test is a file and its name does NOT end with any
    of the strings listed in endings.
    """
    if not os.path.isfile(test):
        return False
    for e in endings:
        if test.endswith(e):
            return False
    return True

#---------------------------------------------------------------------------
# Basic project information
#---------------------------------------------------------------------------


def find_package_data():
    """Find package_data."""
    package_data = {"turbomoleio": ["input/templates/*.yaml"]}
    return package_data


def find_scripts():
    """Find scripts."""
    scripts = []
    #
    # All python files in turbomoleio/scripts
    # scripts for more specific technical tasks should be put in devscripts these are not standard installed
    pyfiles = glob(os.path.join('turbomoleio', 'scripts', "*.py") )
    scripts.extend(pyfiles)
    return scripts


def cleanup():
    """Clean up the junk left around by the build process."""
    if "develop" not in sys.argv:
        import shutil
        try:
            shutil.rmtree('turbomoleio.egg-info')
        except:
            try:
                os.unlink('turbomoleio.egg-info')
            except:
                pass


# List of external packages we rely on.
# Note setup install will download them from Pypi if they are not available.
install_requires = [
    "numpy",
    "matplotlib",
    "pandas",
    "pymatgen>=2019.2.28",
    "monty>=2.0.4",
    "pexpect>=4.7",
    "cerberus",
]

#---------------------------------------------------------------------------
# Find all the packages, package data, and data_files
#---------------------------------------------------------------------------

# Get the set of packages to be included.
my_packages = find_packages(exclude=())

my_scripts = find_scripts()

my_package_data = find_package_data()

about = {}
with open(os.path.join(os.path.dirname(__file__), "turbomoleio", "__version__.py"), "r") as f:
    exec(f.read(), about)

# Create a dict with the basic information
# This dict is eventually passed to setup after additional keys are added.
setup_args = dict(name="turbomoleio",
                  version=about["__version__"],
                  description=about['__description__'],
                  long_description="  ",
                  author="Michiel van Setten, David Waroquiers, Guido Petretto, Jan Kloppenburg",
                  author_email="mjvansetten@gmail.com, david.waroquiers@matgenix.com, guido.petretto@matgenix.com, "
                               "jan.kloppenburg@uclouvain.be",
                  #url=url,
                  #download_url=download_url,
                  license="license",
                  platforms="platforms",
                  keywords="keywords",
                  install_requires=install_requires,
                  packages=my_packages,
                  package_data=my_package_data,
                  scripts=my_scripts,
                  extras_require={
                      "testing": [
                          "pytest >= 3.8.0",
                          "pytest-cov >= 2.6.0",
                      ],
                      "docs": [
                          "sphinx >= 1.8.1",
                      ]
                  },
                  )


if __name__ == "__main__":
    from setuptools import setup
    setup(**setup_args)
    cleanup()
