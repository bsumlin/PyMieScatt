# -*- coding: utf-8 -*-
"""
Documentation is hosted at http://pymiescatt.readthedocs.io/en/latest/index.html
PyMieScatt - the Python Mie Scattering Package
Written by Benjamin Sumlin, Washington University in St. Louis
Aerosol Impacts and Research Laboratory
Department of Energy, Environmental, and Chemical Engineering
Special thanks to Dr. Rajan Chakrabarty, Dr. William Heinson, Claire Fortenberry, and Apoorva Pandey
"""
from setuptools import setup
import re
VERSIONFILE="PyMieScatt/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
  verstr = mo.group(1)
else:
  raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))
    
setup(name='PyMieScatt',
      version=verstr,
      description="A collection of forward and inverse Mie solving routines based on Bohren and Huffman's Mie Theory derivations.",
      long_description=verstr + ' - Added error bounds as an option for the graphical inversion method. Added new automatic inversion methods Inversion() and Inversion_SD(). Significantly improved the iterative methods.\nDocs are hosted at `ReadTheDocs <http://pymiescatt.readthedocs.io/>`_.',
      url='http://air.eece.wustl.edu/people/ben-sumlin/',
      author='Benjamin Sumlin',
      author_email='bsumlin@wustl.edu',
      license='GPL',
      packages=['PyMieScatt'],
      keywords=['Mie Rayleigh scattering absorption extinction light refraction'],
      classifiers = ['Development Status :: 5 - Production/Stable','Intended Audience :: Science/Research','Programming Language :: Python :: 3 :: Only','Topic :: Scientific/Engineering :: Physics'],
      install_requires=['numpy','scipy','matplotlib','shapely'],
      zip_safe=False)