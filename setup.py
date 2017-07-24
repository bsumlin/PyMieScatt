# -*- coding: utf-8 -*-

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
      long_description=verstr + ' - considerably refactored and commented code to be more human-readable. Added a graphical inversion method.',
      url='http://air.eece.wustl.edu/people/ben-sumlin/',
      author='Benjamin Sumlin',
      author_email='bsumlin@wustl.edu',
      license='GPL',
      packages=['PyMieScatt'],
      keywords=['Mie Rayleigh scattering absorption extinction light refraction'],
      classifiers = ['Development Status :: 5 - Production/Stable','Intended Audience :: Science/Research','Programming Language :: Python :: 3 :: Only','Topic :: Scientific/Engineering :: Physics'],
      install_requires=['numpy','scipy','matplotlib'],
      zip_safe=False)