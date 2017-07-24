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

   Return Mie efficencies of a spherical particle with a given refractive index *m*, *wavelength*, and *diameter*. Optionally returns the parameters as a dict when *asDict* is specified and set to True. Uses :py:func:`Mie_ab` to calculate :math:`a_n` and :math:`b_n`, and then calculates :math:`Q_i` via:
   
		:math:`Q_{ext}=\frac{2}{x^2}\sum_{n=1}^{n_{max}}(2n+1)\:\text{Re}\left\{a_n+b_n\right\}`
		
		:math:`Q_{sca}=\frac{2}{x^2}\sum_{n=1}^{n_{max}}(2n+1)(|a_n|^2+|b_n|^2)`
		
		:math:`Q_{abs}=Q_{ext}-Q_{sca}`
		
		:math:`Q_{back}=\frac{1}{x^2}\left|\sum_{n=1}^{n_{max}}(2n+1)(-1)^n(a_n-b_n)\right|^2`
		
		
   
.. py:Function:: Mie_ab(m,x)

   Returns external field coefficients :math:`a_n` and :math:`b_n` based on inputs of *m* and :math:`x=\pi d/\lambda`. Typically not available as a top level call but can be specifically imported via ::

   $ from PyMieScatt import Mie_ab

.. py:Function:: Mie_cd(m,x)

   Returns internal field coefficients :math:`c_n` and :math:`d_n` based on inputs of *m* and :math:`x=\pi d/\lambda`. Typically not available as a top level call but can be specifically imported via ::

   $ from PyMieScatt import Mie_cd

.. py:Function:: RayleighMieQ(m, wavelength, diameter[, asDict=False])

   Returns Mie efficencies of a spherical particle in the Rayleigh regime (:math:`x=\pi d/\lambda << 1`) given refractive index *m*, *wavelength*, and *diameter*. Optionally returns the parameters as a dict when *asDict* is specified and set to True. Uses Rayleigh-regime approximations:
   
		:math:`Q_{sca}=\frac{8x^4}{3}\left|{\frac{m^2-1}{m^2+2}}\right|^2`
   
		:math:`Q_{abs}=4x\:\text{Im}\left\{\frac{m^2-1}{m^2+2}\right\}`
   
		:math:`Q_{ext}=Q_{sca}+Q_{abs}`
   
		:math:`Q_{back}=\frac{3Q_{sca}}{2}`
   
		:math:`Q_{ratio}=1.5`
   
		:math:`Q_{pr}=Q_{ext}`      
   
.. py:Function::`LowFrequencyMieQ(m, wavelength, diameter[, asDict=False])

   Returns Mie efficencies of a spherical particle in the low-frequency regime given refractive index *m*, *wavelength*, and *diameter*. Optionally returns the parameters as a dict when *asDict* is specified and set to True. Uses :py:func:`LowFrequencyMie_ab` to calculate :math:`a_n` and :math:`b_n`.

.. py:Function:: LowFrequencyMie_ab(m,x)

   Returns external field coefficients :math:`a_n` and :math:`b_n` based on inputs of *m* and :math:`x=\pi d/\lambda` by limiting the expansion of :math:`a_n` and :math:`b_n` to second order:
   
		:math:`a_1=-\frac{i2x^3}{3}\frac{(m^2-1)}{m^2+2}`
   
		:math:`a_2=-\frac{ix^5}{15}\frac{(m^2-1)}{2m^2+3}`
   
		:math:`b_1=-\frac{ix^5}{45}(m^2-1)`
   
		:math:`b_2=0`
