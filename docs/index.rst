.. PyMieScatt documentation master file, created by
   sphinx-quickstart on Sat Jun 24 18:09:26 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
   

Online user's guide for the Python Mie Scattering package (PyMieScatt)
======================================================================

Welcome AAAR 2017 Attendees! You can find Ben at posters 2CA.18, 2CA.21 (Tuesday) and 8IM.44 (Thursday). I'm happy to answer any questions, and you can always `e-mail me <http://pymiescatt.readthedocs.io/en/latest/#author-contact-information>`_.

Documentation is currently under development. This documentation includes a brief discussion of Mie theory and the development of the functions included in the package. A manuscript communicating the development and use of this package was submitted to the `Journal of Quantative Spectroscopy and Radiative Transfer <http://www.sciencedirect.com/science/journal/00224073>`_. It has been accepted and is currently being proofed for publication. The JQSRT preprint is `available now <https://doi.org/10.1016/j.jqsrt.2017.10.012>`_.


Install PyMieScatt
------------------

You can install PyMieScatt from `The Python Package Index (PyPI) <https://pypi.python.org/pypi/PyMieScatt>`_ with ::

   $ pip install PyMieScatt


or from `GitHub <https://github.com/bsumlin/PyMieScatt>`_. Clone the repository and then run ::

   $ python setup.py install
   
The current version is 1.3.2.
   

.. toctree::
   :maxdepth: 2
   :caption: Table of Contents
   
   Documentation Home <index>
   Forward Functions for Homogeneous Spheres <forward>
   Forward Functions for Coated Spheres <forwardCS>
   Inverse Mie Functions for Homogeneous Spheres (incomplete) <inverse>
   General Usage tips and Example Scripts (incomplete) <examples>
   
   

Revision Notes - version 1.3.2
------------------------------

- Renamed GraphicalInversion() and GraphicalInversion_SD() to ContourIntersection() and ContourIntersection_SD(), respectively.

Revision History
----------------

- 1.3.1
  - Optimizations to the resolution of the survey-intersection inversion method

Revisions in Progress
---------------------

- Trying to speed up the scattering function of a size distribution

Documentation To-Do List
------------------------

- More example scripts

PyMieScatt To-Do List
---------------------

- Upload package to Anaconda cloud

Publications Using PyMieScatt
-----------------------------

- Sumlin BJ, Pandey A, Walker MJ, Pattison RS, Williams BJ, Chakrabarty RK. Atmospheric Photooxidation Diminishes Light Absorption by Primary Brown Carbon Aerosol from Biomass Burning. Environ Sci Tech Let. 2017. DOI: `10.1021/acs.estlett.7b00393 <http://doi.org/10.1021/acs.estlett.7b00393>`_

- Sumlin BJ, Heinson WR, Chakrabarty RK. Retrieving the Aerosol Complex Refractive Index using PyMieScatt: A Mie Computational Package with Visualization Capabilities. J. Quant. Spectros. Rad. Trans. 2017. DOI: `10.1016/j.jqsrt.2017.10.012 <https://doi.org/10.1016/j.jqsrt.2017.10.012>`_

Author Contact Information
--------------------------
PyMieScatt was written by `Benjamin Sumlin <https://air.eece.wustl.edu/people/ben-sumlin/>`_. Special thanks to Dr. William Heinson, Dr. Rajan Chakrabarty, Claire Fortenberry, and Apoorva Pandey for their insights and support.

Email: `bsumlin@wustl.edu <mailto:bsumlin@wustl.edu?subject=PyMieScatt>`_


