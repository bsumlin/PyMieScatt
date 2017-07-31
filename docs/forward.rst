Functions for Forward Mie Calculations
======================================

Functions for homogeneous spheres
---------------------------------

.. py:function:: MieQ(m, wavelength, diameter[, asDict=False])

   Computes Mie efficencies *Q* and asymmetry parameter *g* of a single, homogeneous particle. Uses :py:func:`Mie_ab` to calculate :math:`a_n` and :math:`b_n`, and then calculates Q\ :sub:`i` via:
   
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

   Computes external field coefficients :math:`a_n` and :math:`b_n` based on inputs of *m* and :math:`x=\pi\,d_p/\lambda`. Typically not available as a top level call but can be specifically imported via ::

   $ from PyMieScatt.Mie import Mie_ab
   
   **Parameters**
   
   
   m : complex
	The complex refractive index with the convention :math:`m=n+ik`.
   x : float
	The size parameter :math:`x=\pi\,d_p/\lambda`.
	
	**Returns**
	
	
   :math:`a_n`, :math:`b_n` : numpy.ndarray
	Arrays of size n\ :sub:`max` = 2+x+4x\ :sup:`1/3`

.. py:Function:: Mie_cd(m,x)

   Computes internal field coefficients :math:`c_n` and :math:`d_n` based on inputs of *m* and :math:`x=\pi\,d_p/\lambda`. Typically not available as a top level call but can be specifically imported via ::

   $ from PyMieScatt.Mie import Mie_cd
   
  **Parameters**
   
   
   m : complex
	The complex refractive index with the convention :math:`m=n+ik`.
   x : float
	The size parameter :math:`x=\pi\,d_p/\lambda`.
	
	**Returns**
	
	
   :math:`c_n`, :math:`d_n` : numpy.ndarray
	Arrays of size n\ :sub:`max` = 2+x+4x\ :sup:`1/3`

.. py:Function:: RayleighMieQ(m, wavelength, diameter[, asDict=False])

   Computes Mie efficencies of a spherical particle in the Rayleigh regime (:math:`x=\pi\,d_p/\lambda \ll 1`) given refractive index *m*, *wavelength*, and *diameter*. Optionally returns the parameters as a dict when *asDict* is specified and set to True. Uses Rayleigh-regime approximations:
   
		:math:`Q_{sca}=\frac{8x^4}{3}\left|{\frac{m^2-1}{m^2+2}}\right|^2`
   
		:math:`Q_{abs}=4x\:\text{Im}\left\{\frac{m^2-1}{m^2+2}\right\}`
   
		:math:`Q_{ext}=Q_{sca}+Q_{abs}`
   
		:math:`Q_{back}=\frac{3Q_{sca}}{2}`
   
		:math:`Q_{ratio}=1.5`
   
		:math:`Q_{pr}=Q_{ext}`      
		
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
	If asDict==True, :py:func:`RayleighMieQ` returns a dict of the above values with appropriate keys.
   
   For example, compute the Mie efficencies of a particle 50 nm in diameter with m = 1.33+0.01i, illuminated by λ = 870 nm: ::
   
		>>> import PyMieScatt as ps
		>>> ps.MieQ(1.33+0.01j,870,50,asDict=True)
		{'Qabs': 0.004057286640269908,
		 'Qback': 0.00017708468873118297,
		 'Qext': 0.0041753430994240295,
		 'Qpr': 0.0041753430994240295,
		 'Qratio': 1.5,
		 'Qsca': 0.00011805645915412197,
		 'g': 0}
   
.. py:Function:: LowFrequencyMieQ(m, wavelength, diameter[, asDict=False])

   Returns Mie efficencies of a spherical particle in the low-frequency regime (:math:`x=\pi\,d_p/\lambda \ll 1`) given refractive index *m*, *wavelength*, and *diameter*. Optionally returns the parameters as a dict when *asDict* is specified and set to True. Uses :py:func:`LowFrequencyMie_ab` to calculate :math:`a_n` and :math:`b_n`, and follows the same math as :py:func:'MieQ'.

.. py:Function:: LowFrequencyMie_ab(m,x)

   Returns external field coefficients :math:`a_n` and :math:`b_n` based on inputs of *m* and :math:`x=\pi\,d_p/\lambda` by limiting the expansion of :math:`a_n` and :math:`b_n` to second order:
   
		:math:`a_1=-\frac{i2x^3}{3}\frac{(m^2-1)}{m^2+2}`
   
		:math:`a_2=-\frac{ix^5}{15}\frac{(m^2-1)}{2m^2+3}`
   
		:math:`b_1=-\frac{ix^5}{45}(m^2-1)`
   
		:math:`b_2=0`
		
   **Parameters**
   
   
   m : complex
	The complex refractive index with the convention :math:`m=n+ik`.
   x : float
	The size parameter :math:`x=\pi\,d_p/\lambda`.
	
	**Returns**
   
   
   :math:`a_n`, :math:`b_n` : numpy.ndarray
	Arrays of size 2.

Functions for polydisperse size distributions of homogeneous spheres
--------------------------------------------------------------------

When an efficency *Q* is integrated over a size distribution :math:`n_d(d_p)`, the result is the *coefficient* :math:`\beta`, which carries units of inverse length. The general form is:

		:math:`\beta=\int\limits_{0}^{\infty}\frac{\pi d_p^2}{4}Q(m,\lambda,d_p)n_d(d_p)(10^{-6})dd_p`
		
where :math:`d_p` is the diameter of the particle (in nm), :math:`n_d(d_p)` is the number of particles of diameter :math:`d_p` (per cubic centimeter), and the factor :math:`10^{-6}` is used to cast the result in units of :math:`\text{Mm}^{-1}`.

The bulk asymmetry parameter *G* is calculated by:

		:math:`G=\frac{\int g(d_p)\beta_{sca}(d_p)dd_p}{\int \beta_{sca}(d_p)dd_p}`
		

.. py:Function:: MieQ_withSizeDistribution(m, wavelength, sizeDistributionDiameterBins, sizeDistribution[, asDict=False])

   Returns Mie coefficients :math:`\beta_{ext}`, :math:`\beta_{sca}`, :math:`\beta_{abs}`, :math:`G`, :math:`\beta_{pr}`, :math:`\beta_{back}`, and :math:`\beta_{ratio}`.
   
   **Parameters**
   
   
   m : complex
	The complex refractive index, with the convention :math:`m=n+ik`.
   wavelength : float
	The wavelength of incident light, in nanometers.
   sizeDistributionDiameterBins : list, tuple, or numpy.ndarray
	The diameter bin midpoints of the size distribution, in nanometers.
   sizeDistribution : list, tuple, or numpy.ndarray
	The number concentrations of the size distribution bins. Must be the same size as sizeDistributionDiameterBins.
   asDict : bool, optional
	If specified and set to True, returns the results as a dict.
	
   **Returns**
   
   
   Bext, Bsca, Babs, G, Bpr, Bback, Bratio : float
	The Mie coefficients calculated by :py:func:`MieQ`, integrated over the size distribution.
   q : dict
	If asDict==True, :py:func:`MieQ_withSizeDistribution` returns a dict of the above values with appropriate keys.

.. py:Function:: MieQ_withLognormalDistribution(m, wavelength, geoStdDev, geoMean, numberOfParticles[, numberOfBins=1000, lower=1, upper=1000, returnDistribution=False, asDict=False])

   Returns Mie coefficients :math:`\beta_{ext}`, :math:`\beta_{sca}`, :math:`\beta_{abs}`, :math:`G`, :math:`\beta_{pr}`, :math:`\beta_{back}`,  and :math:`\beta_{ratio}`.
   
   **Parameters**
   
   
   m : complex
	The complex refractive index, with the convention :math:`m=n+ik`.
   wavelength : float
	The wavelength of incident light, in nanometers.
   geoStdDev : float
	The geometric standard deviation :math:`\sigma_g`.
   geoMean : float
	The geometric mean diameter :math:`d_{pg}`, in nanometers.
   numberOfParticles : float
	The total number of particles in the distribution.
   numberOfBins : int, optional
	The number of discrete bins in the distribution. Defaults to 1000.
   lower : float, optional
	The smallest diameter bin, in nanometers. Defaults to 1 nm.
   upper : float, optional
	The largest diameter bin, in nanometers. Defaults to 1000 nm.
   returnDistribution : bool, optional
	If True, both the size distribution bins and number concentrations will be returned.
   asDict : bool, optional
	If True, returns the results as a dict.
	
   **Returns**
   
   
   Bext, Bsca, Babs, G, Bpr, Bback, Bratio : float
	The Mie coefficients calculated by :py:func:`MieQ`, integrated over the size distribution.
   diameters, nd : numpy.ndarray
	The diameter bins and number concentrations per bin, respectively.
   B : dict
	If asDict==True, :py:func:`MieQ_withLognormalDistribution` returns a dict of the above values with appropriate keys.
   
   For example, compute the Mie coefficients of a lognormal size distribution with 1000000 particles, σ\ :sub:`g` = 1.7, and d\ :sub:`pg` = 200 nm; with m = 1.60+0.08i and λ = 532 nm: ::
   
		>>> import PyMieScatt as ps
		>>> ps.MieQ_withLognormalDistribution(1.60+0.08j,532,1.7,200,1e6,asDict=True)
		{'Babs': 33537.324569179938,
		'Bback': 10188.473118449627,
		'Bext': 123051.1109783932,
		'Bpr': 62038.347528346232,
		'Bratio': 12701.828124508347,
		'Bsca': 89513.786409213266,
		'bigG': 0.6816018615403715}
		
Functions for Coated Spheres (Core-Shell Particles)
---------------------------------------------------

.. py:Function:: MieQCoreShell(mCore, mShell, wavelength, dCore, dShell[, asDict=False])

   Compute Mie efficencies *Q* and asymmetry parameter *g* of a single, coated particle. Uses :py:func:`CoreShell_ab` to calculate a\ :sub:`n` and b\ :sub:`n` , and then calculates Q\ :sub:`i` following closely from the original BHMIE.
   
   **Parameters**
   
   
   mCore : complex
	The complex refractive index of the core region, with the convention :math:`m=n+ik`.
   mShell : complex
	The complex refractive index of the shell region, with the convention :math:`m=n+ik`.
   wavelength : float
	The wavelength of incident light, in nanometers.
   dCore : float
	The diameter of the core, in nanometers.
   dShell : float
	The diameter of the shell, in nanomaters. This is equal to the total diameter of the particle.
   asDict : bool, optional
	If True, returns the results as a dict.
	
   **Returns**
   
   
   qext, qsca, qabs, g, qpr, qback, qratio : float
	The Mie efficencies described above.
   q : dict
	If asDict==True, :py:func:`MieQCoreShell` returns a dict of the above values with appropriate keys.
	
   **Considerations**
   
   
   When using this function in a script, there are three simplifying clauses that can speed up computation when considering both coated and homogeneous particles. Upon determining the size parameters of the core and the shell:
   
   - if x\ :sub:`core` == x\ :sub:`shell`, then :py:func:`MieQCoreShell` returns Mie efficencies calculated by :py:func:`MieQ(mCore,wavelength,dShell)`.
   - If x\ :sub:`core`==0, then :py:func:`MieQCoreShell` returns efficencies calculated by :py:func:`MieQ(mShell,wavelength,dShell)`.
   - If m\ :sub:`core`==m\ :sub:`shell`, returns efficencies calculated by :py:func:`MieQ(mCore,wavelength,dShell)`.
   
.. py:Function: CoreShell_ab(m, x)

   Computes external field coefficients :math:`a_n` and :math:`b_n` based on inputs of *m* and :math:`x=\pi\,d_p/\lambda`. Typically not available as a top level call but can be specifically imported via ::

   $ from PyMieScatt.CoreShell import CoreShell_ab
   
   **Parameters**
   
   
   m : complex
	The complex refractive index with the convention :math:`m=n+ik`.
   x : float
	The size parameter :math:`x=\pi\,d_p/\lambda`.
	
	**Returns**
	
	
   :math:`a_n`, :math:`b_n` : numpy.ndarray
	Arrays of size n\ :sub:`max` = 2+x+4x\ :sup:`1/3`

Angular Functions
-----------------

These functions compute the angle-dependent scattering functions and create arrays that are useful for plotting.

Homogeneous Spheres
~~~~~~~~~~~~~~~~~~~

.. py:Function: ScatteringFunction(m, wavelength, diameter[, minAngle=0, maxAngle=180, angularResolution=0.5, normed=False])

   Creates arrays for plotting the angular scattering intensity functions in theta-space with parallel, perpendicular, and unpolarized light. Uses :py:func:`MieS1S2` to compute S\ :sub:`1` and S\ :sub:`2`, then computes parallel, perpendicular, and unpolarized intensities by
   
		:math:`SR(\Theta)=|S_1|^2`
		:math:`SL(\Theta)=|S_2|^2`
		:math:`SU(\Theta)=\frac{1}{2}(SR+SL)`
   
   **Parameters**
   
   
   m : complex
	The complex refractive index with the convention :math:`m=n+ik`.
   wavelength : float
	The wavelength of incident light, in nanometers.
   diameter : float
	The diameter of the particle, in nanometers.
   minAngle : float, optional
	The minimum scattering angle (in degrees) to be calculated. Defaults to 0.
   maxAngle : float, optional
	The maximum scattering angle (in degrees) to be calculated. Defaults to 180.
   angularResolution : float, optional
	The resolution of the output. Defaults to 0.5, meaning a value will be calculated for every 0.5 degrees.
   normed : bool, optional
	If True, will normalize the output such that the maximum intensity will be 1.0. Defaults to False.
	
   **Returns**
   
   
   theta : numpy.ndarray
	An array of the angles used in calculations. Values will be spaced according to *angularResolution*, and the size of the array will be *(maxAngle-minAngle)/angularResolution*.
   SL : numpy.ndarray
	An array of the scattered intensity of left-polarized (parallel) light. Same size as the *theta* array.
   SR : numpy.ndarray
	An array of the scattered intensity of right-polarized (perpendicular) light. Same size as the *theta* array.
   SU : numpy.ndarray
	An array of the scattered intensity of unpolarized light, which is the average of SL and SR. Same size as the *theta* array.


.. py:Function: qSpaceScatteringFunction(m, wavelength, diameter[, normed=False])

   Creates arrays for plotting the angular scattering intensity functions in q-space with parallel, perpendicular, and unpolarized light.
   
   **Parameters**
   
   
   m : complex
	The complex refractive index with the convention :math:`m=n+ik`.
   wavelength : float
	The wavelength of incident light, in nanometers.
   diameter : float
	The diameter of the particle, in nanometers.
   normed : bool, optional
	If True, will normalize the output such that the maximum intensity will be 1.0. Defaults to False.
	
   **Returns**
   
   
   qR : numpy.ndarray
	An array of the q-space angles used in calculations. Size is 3600.
   SL : numpy.ndarray
	An array of the scattered intensity of left-polarized (parallel) light. Same size as the *qR* array.
   SR : numpy.ndarray
	An array of the scattered intensity of right-polarized (perpendicular) light. Same size as the *qR* array.
   SU : numpy.ndarray
	An array of the scattered intensity of unpolarized light, which is the average of SL and SR. Same size as the *qR* array.