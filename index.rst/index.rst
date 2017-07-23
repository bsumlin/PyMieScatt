.. PyMieScatt documentation master file, created by
   sphinx-quickstart on Sat Jun 24 18:09:26 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the online user's guide for PyMieScatt!
==================================================

Documentation is currently under development and will be available Monday, July 24th.

Install PyMieScatt
------------------

You can install PyMieScatt from `PyPI <https://pypi.python.org/pypi/PyMieScatt>`_ with ::

   $ pip install PyMieScatt


.. toctree::
   :maxdepth: 2
   :caption: Contents:
   Introduction <Introduction>
   Theory <Theory>
   Functions <Function>
   Examples <Examples>

.. py:function:: MieQ(m, wavelength, diameter[, asDict=False])

   Return Mie efficencies of a spherical particle with a given refractive index *m*, *wavelength*, and *diameter*. Optionally returns the parameters as a dict when *asDict* is specified and set to True.
   
.. py:Function:: Mie_ab(m,x)

   Returns external field coefficients :math:`a_n` and :math:`b_n` based on inputs of *m* and :math:`x=\pi d/\lambda`. Typically not available as a top level call but can be specifically imported via ::

   $ from PyMieScatt import Mie_ab

.. py:Function:: Mie_cd(m,x)

   Returns internal field coefficients :math:`c_n` and :math:`d_n` based on inputs of *m* and :math:`x=\pi d/\lambda`. Typically not available as a top level call but can be specifically imported via ::

   $ from PyMieScatt import Mie_cd

.. py:Function:: RayleighMieQ(m, wavelength, diameter[, asDict=False])

   Returns Mie efficencies of a spherical particle in the Rayleigh regime (:math:`x=\pi d/\lambda << 1`) given refractive index *m*, *wavelength*, and *diameter*. Optionally returns the parameters as a dict when *asDict* is specified and set to True. Uses :py:func:`RayleighMie_ab` to calculate :math:`a_n` and :math:`b_n`.
   
.. py:Function::`LowFrequencyMieQ(m, wavelength, diameter[, asDict=False])

   Returns Mie efficencies of a spherical particle in the low-frequency regime given refractive index *m*, *wavelength*, and *diameter*. Optionally returns the parameters as a dict when *asDict* is specified and set to True. Uses :py:func:`LowFrequencyMie_ab` to calculate :math:`a_n` and :math:`b_n`.

