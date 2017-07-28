.. PyMieScatt documentation master file, created by
   sphinx-quickstart on Sat Jun 24 18:09:26 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the online user's guide for PyMieScatt!
==================================================

Documentation is currently under development. Migrating to readthedocs broke things in new and exciting ways. Documentation is scheduled to be complete on Monday, July 31st.

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
   
Functions for homogeneous spheres
---------------------------------

.. py:function:: MieQ(m, wavelength, diameter[, asDict=False])

   Compute Mie efficencies *Q* and asymmetry parameter *g* of a single, homogeneous particle. Uses :py:func:`Mie_ab` to calculate :math:`a_n` and :math:`b_n`, and then calculates :math:`Q_i` via:
   
		:math:`Q_{ext}=\frac{2}{x^2}\sum_{n=1}^{n_{max}}(2n+1)\:\text{Re}\left\{a_n+b_n\right\}`
		
		:math:`Q_{sca}=\frac{2}{x^2}\sum_{n=1}^{n_{max}}(2n+1)(|a_n|^2+|b_n|^2)`
		
		:math:`Q_{abs}=Q_{ext}-Q_{sca}`
		
		:math:`Q_{back}=\frac{1}{x^2}\left|\sum_{n=1}^{n_{max}}(2n+1)(-1)^n(a_n-b_n)\right|^2`
		
		:math:`Q_{ratio}=\frac{Q_{back}}{Q_{sca}}`
		
		:math:`g=\frac{4}{Q_{sca}x^2}\left[\sum\limits_{n=1}^{n_{max}}\frac{n(n+2)}{n+1}\text{Re}\left\{a_n a_{n+1}^*+b_n b_{n+1}^*\right\}+\sum\limits_{n=1}^{n_{max}}\frac{2n+1}{n(n+1)}\text{Re}\left\{a_n b_n^*\right\}\right]`
		
		:math:`Q_{pr}=Q_{ext}-gQ_{sca}`
		
   where asterisks denote the complex conjugates.
   
   **Parameters**
   
   m : complex
	The complex refractive index, with the convention :math:`m=n+ik`.
   wavelength : float
	The wavelength of incident light, in nanometers.
   diameter : float
	The diameter of the particle, in nanometers.
   asDict : bool, optional
	If specified and set to True, returns the results as a dict.
	
   **Returns**
   
   qext, qsca, qabs, g, qpr, qback, qratio : float
	The Mie efficencies described above.
   q : dict
	If asDict==True, :py:func:`MieQ` returns a dict of the above values with appropriate keys.
   
   For example, compute the Mie efficencies of a particle 300 nm in diameter with m=1.77+0.63i, illuminated by :math:`\lambda` = 375 nm: ::
   
		>>> import PyMieScatt as ps
		>>> ps.MieQ(1.77+0.63j,375,300,asDict=True)
		{'Qext': 2.8584971991564112,
		 'Qsca': 1.3149276685170939,
		 'Qabs': 1.5435695306393173,
		 'g': 0.7251162362148782,
		 'Qpr': 1.9050217972664911,
		 'Qback': 0.20145510481352547,
		 'Qratio': 0.15320622543498222}
   
.. py:Function:: Mie_ab(m,x)

   Returns external field coefficients :math:`a_n` and :math:`b_n` based on inputs of *m* and :math:`x=\frac{\pi\,d_p}{\lambda}`. Typically not available as a top level call but can be specifically imported via ::

   $ from PyMieScatt import Mie_ab

.. py:Function:: Mie_cd(m,x)

   Returns internal field coefficients :math:`c_n` and :math:`d_n` based on inputs of *m* and :math:`x=\frac{\pi\,d_p}{\lambda}`. Typically not available as a top level call but can be specifically imported via ::

   $ from PyMieScatt import Mie_cd

.. py:Function:: RayleighMieQ(m, wavelength, diameter[, asDict=False])

   Returns Mie efficencies of a spherical particle in the Rayleigh regime (:math:`x=\frac{\pi\,d_p}{\lambda} \ll 1`) given refractive index *m*, *wavelength*, and *diameter*. Optionally returns the parameters as a dict when *asDict* is specified and set to True. Uses Rayleigh-regime approximations:
   
		:math:`Q_{sca}=\frac{8x^4}{3}\left|{\frac{m^2-1}{m^2+2}}\right|^2`
   
		:math:`Q_{abs}=4x\:\text{Im}\left\{\frac{m^2-1}{m^2+2}\right\}`
   
		:math:`Q_{ext}=Q_{sca}+Q_{abs}`
   
		:math:`Q_{back}=\frac{3Q_{sca}}{2}`
   
		:math:`Q_{ratio}=1.5`
   
		:math:`Q_{pr}=Q_{ext}`      
   
.. py:Function::`LowFrequencyMieQ(m, wavelength, diameter[, asDict=False])

   Returns Mie efficencies of a spherical particle in the low-frequency regime given refractive index *m*, *wavelength*, and *diameter*. Optionally returns the parameters as a dict when *asDict* is specified and set to True. Uses :py:func:`LowFrequencyMie_ab` to calculate :math:`a_n` and :math:`b_n`.

.. py:Function:: LowFrequencyMie_ab(m,x)

   Returns external field coefficients :math:`a_n` and :math:`b_n` based on inputs of *m* and :math:`x=\frac{\pi\,d_p}{\lambda}` by limiting the expansion of :math:`a_n` and :math:`b_n` to second order:
   
		:math:`a_1=-\frac{i2x^3}{3}\frac{(m^2-1)}{m^2+2}`
   
		:math:`a_2=-\frac{ix^5}{15}\frac{(m^2-1)}{2m^2+3}`
   
		:math:`b_1=-\frac{ix^5}{45}(m^2-1)`
   
		:math:`b_2=0`

Functions for polydisperse size distributions of homogeneous spheres
--------------------------------------------------------------------

When an efficency *Q* is integrated over a size distribution :math:`n_d(d_p)`, the result is the *coefficient* :math:`\beta`, which carries units of inverse length. The general form is:

		:math:`\beta=\int\limits_{0}^{\infty}\frac{\pi d_p^2}{4}Q(m,\lambda,d_p)n_d(d_p)(10^{-6})dd_p`
		
where :math:`d_p` is the diameter of the particle (in nm), :math:`n_d(d_p)` is the number of particles of diameter :math:`d_p` (per cubic centimeter), and the factor :math:`10^{-6}` is used to cast the result in units of :math:`\text{Mm}^{-1}`.

The bulk asymmetry parameter *G* is calculated by:

		:math:`G=\frac{\int g(d_p)\beta_{sca}(d_p)dd_p}{\int \beta_{sca}(d_p)dd_p}`

.. py:Function::`MieQ_withSizeDistribution(m, wavelength, sizeDistributionDiameterBins, sizeDistribution[, asDict=False])

   Returns Mie coefficients :math:`\beta_{ext}`, :math:`\beta_{sca}`, :math:`\beta_{abs}`, :math:`G`, :math:`\beta_{pr}`, :math:`\beta_{back}`,  and :math:`\beta_{ratio}`.
