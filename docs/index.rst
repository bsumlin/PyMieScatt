.. PyMieScatt documentation master file, created by
   sphinx-quickstart on Sat Jun 24 18:09:26 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
   

Online user's guide for the Python Mie Scattering package (PyMieScatt)
======================================================================

Documentation is always under development, but mostly complete. A manuscript communicating the development of the inverse Mie algorithms was published by the `Journal of Quantative Spectroscopy and Radiative Transfer <http://www.sciencedirect.com/science/journal/00224073>`_. The JQSRT article is `available here <https://doi.org/10.1016/j.jqsrt.2017.10.012>`_.

**NOTE TO USERS:** When using PyMieScatt, pay close attention to the units of the your inputs and outputs. Wavelength and particle diameters are always in nanometers, efficiencies are unitless, cross-sections are in nm\ :sup:`2`, coefficients are in Mm\ :sup:`-1`, and size distribution concentration is always in cm\ :sup:`-3`. If you use other units, your outputs may not make sense.

**NOTE TO THOSE WITH MIEPLOT EXPERIENCE:** The functions in PyMieScatt take particle *diameter*. MiePlot's default is to take the particle *radius* in micrometers. Make sure all your particle dimensions, whether for a single particle or for a distribution, are for the diamaters, in nanometers.


Install PyMieScatt
------------------

NOTE: You must install `Shapely <https://shapely.readthedocs.io/>`_ first, preferably from GitHub. Users have reported difficulty installing it with pip. Conda works, too.

The current version is 1.7.5. You can install PyMieScatt from `The Python Package Index (PyPI) <https://pypi.python.org/pypi/PyMieScatt>`_ with ::

   $ pip install PyMieScatt


or from `GitHub <https://github.com/bsumlin/PyMieScatt>`_. Clone the repository and then run ::

   $ python setup.py install

Revision Notes - version 1.7.5 (23 February, 2020)
------------------------------------------------------------------------------

  - Fixed :py:func:`AutoMieQ` per discussions with Gerard van Ewijk. In the case of nMedium!=1, :py:func:`AutoMieQ` was calculating effective n and wavelength, and then passing those parameters to the relevant Mie function. Those functions then re-calculated the effective n and wavelength, leading to errors.
  - Fixed :py:func:`ContourIntersection_SD` per discussions with Hans Moosmuller. The inputs should now correctly scale for units of Mm-1.

Revision History
----------------

- 1.7.4 (6 May, 2019)

  - Fixed :py:func:`ScatteringFunction` per discussions with @zcm73400 on GitHub. View the pull request for more info.

- 1.7.3 (23 August, 2018) - 1.7.2 was skipped ¯\\_(ツ)_/¯

  - Added :py:func:`CoreShellS1S2` to __init__.py. Also added :py:func:`CoreShellMatrixElements` to the documentation. Thanks Jonathan Taylor for the heads up!

- 1.7.1 (12 April, 2018)

  - Fixed a bug in :py:func:`MieQ_withWavelengthRange` where the inputs would be affected by in-place math performed within the function. This bug was also present in :py:func:`MieQ_withSizeParameterRange` and has been fixed.

- 1.7.0 (5 April, 2018)

  - Updated most of the forward homogeneous sphere functions with a new optional parameter **nMedium**, which allows for Mie calculations in media other than vacuum/air. Please see documentation.

- 1.6.0b0 (23 March, 2018)

  - Updated :py:func:`ContourIntersection` and :py:func:`ContourIntersection_SD` to take optional constraint parameters of an assumed *n* or *k*. Please see the documentation for more information.

- 1.5.2 (9 March, 2018)

  - Fixed a bug in :py:func:`ContourIntersection` and :py:func:`ContourIntersection_SD` that would occasionally cause a single solution from two optical measurements to not be reported (thanks to Miriam Elser for pointing this bug out).

- 1.5.1 (7 March, 2018)

  - Added the option to report single-particle Mie efficiencies as optical cross-sections. This affects :py:func:`MieQ`, :py:func:`RayleighMieQ`, :py:func:`AutoMieQ`, :py:func:`LowFrequencyMieQ`, and :py:func:`MieQCoreShell`. The results carry units of nm\ :sup:`2`.

- 1.4.3 (21 February, 2018)

  - Fixed a small bug in :py:func:`ContourIntersection` and :py:func:`ContourIntersection_SD` that would produce an error if no intersections were detected. Now it just throws a warning. I'll update soon to have better reporting.

- 1.4.2 (25 January, 2018)

  - Very minor adjustment to :py:func:`AutoMieQ`; changed the crossover from Rayleigh to Mie to x=0.01 (previously 0.5). Thanks to `John Kendrick <https://github.com/JohnKendrick/PDielec>`_ for the suggestion.

- 1.4.1 (25 January, 2018)

  - Added `Shapely <https://shapely.readthedocs.io/>`_ support! Shapely is a geometric manipulation and analysis package. I wrote it in as a slightly faster, more robust way to look for intersections in n-k space when doing inversions. It also makes the code more readable and makes it clearer how the intersection method works, especially when including backscatter to find a unique solution. There is no change to the user experience, other than slight speedups.

- 1.3.7

  - Fixed a major bug in :py:func:`ContourIntersection` and :py:func:`ContourIntersection_SD` that prevented them from using the actual input values to derive solutions.

- 1.3.6

  - Added new normalization options to :py:func:`ScatteringFunction` and :py:func:`SF_SD`. Docs for those functions have details.
  
- 1.3.5

  - Fixed a bug that prevented SF_SD from properly scaling with the number of particles.
  
- 1.3.4.1

  - Added a new sub-version delimiter. 1.x.y.z will be for minor revisions including some optimizations I've been working on that don't merit a full 1.x.y release.
  - Added a new AutoMie_ab() function that uses LowFrequencyMie_ab() for *x = πd/λ* < 0.5 and Mie_ab() otherwise.
  - Sped up the MieS1S2() function by using the new AutoMie_ab() function.
  - Sped up the SF_SD() function by about 33% (on average) when the MieS1S2() optimizations are considered.
  - Added Mie_cd() to __init__.py.

- 1.3.4

  - Fixed a really dumb bug introduced in 1.3.3.
  
- 1.3.3

  - Fixed a big that caused SF_SD() to throw errors when a custom angle range was specified.
  - Added MieS1S2() and MiePiTau() to __init__.py. Dunno why they weren't always there.
  
- 1.3.2

  - Renamed GraphicalInversion() and GraphicalInversion_SD() to ContourIntersection() and ContourIntersection_SD(), respectively.
  
- 1.3.1

  - Optimizations to the resolution of the survey-intersection inversion method.

Revisions in Progress
---------------------

- Would like to re-write the inversion functions to be as general as possible, i.e., if I pass scattering, absorption, particle size, and refractive index, it would solve for the wavelength.
- Ablility to pass array objects directly to all functions (within reason).
- Auto-graphing capabilities for sacttering functions.

Documentation To-Do List
------------------------

- More example scripts, I guess?
- As a few function names and parameter names get updated, there may be some typos in old examples. I'll catch those as they crop up.

PyMieScatt To-Do List
---------------------

- Upload package to Anaconda cloud.

Publications Using PyMieScatt
-----------------------------

If you use PyMieScatt in your research, please let me know and I'll link the publications here.

- Sumlin BJ, Pandey A, Walker MJ, Pattison RS, Williams BJ, Chakrabarty RK. Atmospheric Photooxidation Diminishes Light Absorption by Primary Brown Carbon Aerosol from Biomass Burning. Environ Sci Tech Let. 2017 4 (12) 540-545 (`Cover article <http://pubs.acs.org/toc/estlcu/4/12>`_). DOI: `10.1021/acs.estlett.7b00393 <http://doi.org/10.1021/acs.estlett.7b00393>`_

- Sumlin BJ, Heinson WR, Chakrabarty RK. Retrieving the Aerosol Complex Refractive Index using PyMieScatt: A Mie Computational Package with Visualization Capabilities. J. Quant. Spectros. Rad. Trans. 2018 (205) 127-134. DOI: `10.1016/j.jqsrt.2017.10.012 <https://doi.org/10.1016/j.jqsrt.2017.10.012>`_

- Sumlin BJ, Heinson YW, Shetty N, Pandey A, Pattison RS, Baker S, Hao WM, Chakrabarty RK. UV-Vis-IR Spectral Complex Refractive Indices and Optical Properties of Brown Carbon Aerosol fro Biomass Burning. J. Quant. Spectros. Rad. Trans. 2018 (206) 392-398 DOI: `10.1016/j.jqsrt.2017.12.009 <https://doi.org/10.1016/j.jqsrt.2017.12.009>`_

- Sumlin BJ, Oxford C, Seo B, Pattison R, Williams B, Chakrabarty RK. Density and Homogeneous Internal Composition of Primary Brown Carbon Aerosol. Environ. Sci. Tech., In press. DOI: `10.1021/acs.est.8b00093 <https://doi.org/10.1021/acs.est.8b00093>`_

Author Contact Information
--------------------------
PyMieScatt was written by `Benjamin Sumlin <https://air.eece.wustl.edu/people/ben-sumlin/>`_. Special thanks to Dr. William Heinson, Dr. Rajan Chakrabarty, Claire Fortenberry, and Apoorva Pandey for their insights and support.

Email: `bsumlin@wustl.edu <mailto:bsumlin@wustl.edu?subject=PyMieScatt>`_


.. toctree::
   :maxdepth: 2
   :caption: Table of Contents
   
   Documentation Home <index>
   Forward Functions for Homogeneous Spheres <forward>
   Forward Functions for Coated Spheres <forwardCS>
   Inverse Mie Functions for Homogeneous Spheres <inverse>
   General Usage tips and Example Scripts (constantly updating) <examples>
