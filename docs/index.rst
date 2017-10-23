.. PyMieScatt documentation master file, created by
   sphinx-quickstart on Sat Jun 24 18:09:26 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
   

Online user's guide for the Python Mie Scattering package (PyMieScatt)
======================================================================

AAAR 2017: thanks everyone for the feedback and comments! To continue the dialog, `e-mail me <http://pymiescatt.readthedocs.io/en/latest/#author-contact-information>`_.

NEWS: PyMieScatt is being ported to Julia! I did a few timing and comparison tests and got some basic Mie functions to run 2-6 times faster on Julia than Python. This is a side project and could take some time to see complete functionality. The inverse algorithms are the hard part. I may port the numerical algorithms first, and then the graphical ones later once I figure out how plotting objects in Julia work.

Documentation is currently under development. A manuscript communicating the development and use of this package was submitted to the `Journal of Quantative Spectroscopy and Radiative Transfer <http://www.sciencedirect.com/science/journal/00224073>`_. It has been accepted and is currently being proofed for publication. The JQSRT preprint is `available now <https://doi.org/10.1016/j.jqsrt.2017.10.012>`_, and the arXiv preprint can be found `here <https://arxiv.org/abs/1710.05288>`_.


Install PyMieScatt
------------------

You can install PyMieScatt from `The Python Package Index (PyPI) <https://pypi.python.org/pypi/PyMieScatt>`_ with ::

   $ pip install PyMieScatt


or from `GitHub <https://github.com/bsumlin/PyMieScatt>`_. Clone the repository and then run ::

   $ python setup.py install
   
The current version is 1.3.3.
   

.. toctree::
   :maxdepth: 2
   :caption: Table of Contents
   
   Documentation Home <index>
   Forward Functions for Homogeneous Spheres <forward>
   Forward Functions for Coated Spheres <forwardCS>
   Inverse Mie Functions for Homogeneous Spheres <inverse>
   General Usage tips and Example Scripts (constantly updating) <examples>
   
   

Revision Notes - version 1.3.3
------------------------------

- Fixed a big that caused SF_SD() to throw errors when a custom angle range was specified.
- Added MieS1S2() and MiePiTau() to __init__.py. Dunno why they weren't always there.

Revision History
----------------

- 1.3.2
  - Renamed GraphicalInversion() and GraphicalInversion_SD() to ContourIntersection() and ContourIntersection_SD(), respectively.
- 1.3.1
  - Optimizations to the resolution of the survey-intersection inversion method.

Revisions in Progress
---------------------

- Trying to speed up the scattering function of a size distribution.
- Would like to be able to pass array objects directly to all functions (within reason).

Documentation To-Do List
------------------------

- More example scripts, I guess?

PyMieScatt To-Do List
---------------------

- Upload package to Anaconda cloud.

Publications Using PyMieScatt
-----------------------------

If you use PyMieScatt in your research, please let me know and I'll link the publications here.

- Sumlin BJ, Pandey A, Walker MJ, Pattison RS, Williams BJ, Chakrabarty RK. Atmospheric Photooxidation Diminishes Light Absorption by Primary Brown Carbon Aerosol from Biomass Burning. Environ Sci Tech Let. 2017. DOI: `10.1021/acs.estlett.7b00393 <http://doi.org/10.1021/acs.estlett.7b00393>`_

- Sumlin BJ, Heinson WR, Chakrabarty RK. Retrieving the Aerosol Complex Refractive Index using PyMieScatt: A Mie Computational Package with Visualization Capabilities. J. Quant. Spectros. Rad. Trans. 2017. DOI: `10.1016/j.jqsrt.2017.10.012 <https://doi.org/10.1016/j.jqsrt.2017.10.012>`_

Author Contact Information
--------------------------
PyMieScatt was written by `Benjamin Sumlin <https://air.eece.wustl.edu/people/ben-sumlin/>`_. Special thanks to Dr. William Heinson, Dr. Rajan Chakrabarty, Claire Fortenberry, and Apoorva Pandey for their insights and support.

Email: `bsumlin@wustl.edu <mailto:bsumlin@wustl.edu?subject=PyMieScatt>`_


