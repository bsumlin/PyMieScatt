Functions for Forward Mie Calculations of Homogeneous Spheres
=============================================================

Functions for single particles
---------------------------------

.. py:function:: MieQ(m, wavelength, diameter[, asDict=False])

   Computes Mie efficencies *Q* and asymmetry parameter *g* of a single, homogeneous particle. Uses :py:func:`Mie_ab` to calculate :math:`a_n` and :math:`b_n`, and then calculates *Q* via:
   
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
	The complex refractive index, with the convention *m = n+ik*.
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
   
   For example, compute the Mie efficencies of a particle 300 nm in diameter with m = 1.77+0.63i, illuminated by λ = 375 nm: ::
   
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

   Computes external field coefficients a\ :sub:`n` and b\ :sub:`n` based on inputs of *m* and :math:`x=\pi\,d_p/\lambda`. Must be explicitly imported via ::

   >>> from PyMieScatt.Mie import Mie_ab
   
   **Parameters**
   
   
   m : complex
	The complex refractive index with the convention *m = n+ik*.
   x : float
	The size parameter :math:`x=\pi\,d_p/\lambda`.
	
	**Returns**
	
	
   an, bn : numpy.ndarray
	Arrays of size n\ :sub:`max` = 2+x+4x\ :sup:`1/3`

.. py:Function:: Mie_cd(m,x)

   Computes internal field coefficients c\ :sub:`n` and d\ :sub:`n` based on inputs of *m* and :math:`x=\pi\,d_p/\lambda`. Must be explicitly imported via ::

   >>> from PyMieScatt.Mie import Mie_cd
   
  **Parameters**
   
   
   m : complex
	The complex refractive index with the convention *m = n+ik*.
   x : float
	The size parameter :math:`x=\pi\,d_p/\lambda`.
	
	**Returns**
	
	
   cn, dn : numpy.ndarray
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
	The complex refractive index, with the convention *m = n+ik*.
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

   Returns Mie efficencies of a spherical particle in the low-frequency regime (:math:`x=\pi\,d_p/\lambda \ll 1`) given refractive index *m*, *wavelength*, and *diameter*. Optionally returns the parameters as a dict when *asDict* is specified and set to True. Uses :py:func:`LowFrequencyMie_ab` to calculate a\ :sub:`n` and b\ :sub:`n`, and follows the same math as :py:func:`MieQ`.

.. py:Function:: LowFrequencyMie_ab(m,x)

   Returns external field coefficients a\ :sub:`n` and b\ :sub:`n` based on inputs of *m* and :math:`x=\pi\,d_p/\lambda` by limiting the expansion of a\ :sub:`n` and b\ :sub:`n` to second order:
   
		:math:`a_1=-\frac{i2x^3}{3}\frac{(m^2-1)}{m^2+2}`
   
		:math:`a_2=-\frac{ix^5}{15}\frac{(m^2-1)}{2m^2+3}`
   
		:math:`b_1=-\frac{ix^5}{45}(m^2-1)`
   
		:math:`b_2=0`
		
   **Parameters**
   
   
   m : complex
	The complex refractive index with the convention *m = n+ik*.
   x : float
	The size parameter :math:`x=\pi\,d_p/\lambda`.
	
	**Returns**
   
   
   an, bn : numpy.ndarray
	Arrays of size 2.

Functions for single particles across various ranges
----------------------------------------------------

.. py:Function:: MieQ_withDiameterRange(m, wavelength[, diameterRange=[10,1000], nd=1000, logD=False])

   Computes the Mie efficencies of particles across a diameter range using :py:func:`MieQ`.
   
   **Parameters**
   
   
   m : complex
	The complex refractive index with the convention *m = n+ik*.
   wavelength : float
	The wavelength of incident light, in nanomaters
   diameterRange : list, optional
	The diameter range, in nanometers. Convention is [*smallest*, *largest*]. Defaults to [10, 1000].
   nd : int, optional
	The number of diameter bins in the range. Defaults to 1000.
   logD : bool, optional
	If True, will use logarithmically-spaced diameter bins. Defaults to False.
	
   **Returns**
   
   
   diameters : numpy.ndarray
	An array of the diameter bins that calculations were performed on. Size is equal to *nd*.
   qext, qsca, qabs, g, qpr, qback, qratio : numpy.ndarray
	The Mie efficencies at each diameter in *diameters*.
	
.. py:Function:: MieQ_withWavelengthRange(m, diameter[, wavelengthRange=[100,1600], nw=1000, logW=False])

   Computes the Mie efficencies of particles across a wavelength range using :py:func:`MieQ`.
   
   **Parameters**
   
   
   m : complex
	The complex refractive index with the convention *m = n+ik*.
   diameter : float
	The diameter of the particle, in nanometers.
   wavelengthRange : list, optional
	The wavelength range of incident light, in nanomaters. Convention is [*smallest*, *largest*]. Defaults to [100, 1600].
   nw : int, optional
	The number of wavelength bins in the range. Defaults to 1000.
   logW : bool, optional
	If True, will use logarithmically-spaced wavelength bins. Defaults to False.
	
   **Returns**
   
   
   wavelengths : numpy.ndarray
	An array of the wavelength bins that calculations were performed on. Size is equal to *nw*.
   qext, qsca, qabs, g, qpr, qback, qratio : numpy.ndarray
	The Mie efficencies at each wavelength in *wavelengths*.
	
.. py:Function:: MieQ_withSizeParameterRange(m[, xRange=[1,10], nx=1000, logX=False])

   Computes the Mie efficencies of particles across a size parameter range (\ :math:`x=\pi\,d_p/\lambda`\ ) using :py:func:`MieQ`.
   
   **Parameters**
   
   
   m : complex
	The complex refractive index with the convention *m = n+ik*.
   xRange : list, optional
	The size parameter range. Convention is [*smallest*, *largest*]. Defaults to [1, 10].
   nx : int, optional
	The number of size parameter bins in the range. Defaults to 1000.
   logX : bool, optional
	If True, will use logarithmically-spaced size parameter bins. Defaults to False.
	
   **Returns**
   
   
   xValues : numpy.ndarray
	An array of the size parameter bins that calculations were performed on. Size is equal to *nx*.
   qext, qsca, qabs, g, qpr, qback, qratio : numpy.ndarray
	The Mie efficencies at each size parameter in *xValues*.


Functions for polydisperse size distributions of homogeneous spheres
--------------------------------------------------------------------

When an efficiency *Q* is integrated over a size distribution n\ :sub:`d`\ (d\ :sub:`p`), the result is the *coefficient* :math:`\beta`, which carries units of inverse length. The general form is:

		:math:`\beta=\int\limits_{0}^{\infty}\frac{\pi d_p^2}{4}Q(m,\lambda,d_p)n_d(d_p)(10^{-6})dd_p`
		
where d\ :sub:`p` is the diameter of the particle (in nm), n\ :sub:`d`\ (d\ :sub:`p`) is the number of particles of diameter d\ :sub:`p` (per cubic centimeter), and the factor 10\ :sup:`-6` is used to cast the result in units of Mm\ :sup:`-1`.

The bulk asymmetry parameter *G* is calculated by:

		:math:`G=\frac{\int g(d_p)\beta_{sca}(d_p)dd_p}{\int \beta_{sca}(d_p)dd_p}`
		

.. py:Function:: MieQ_withSizeDistribution(m, wavelength, sizeDistributionDiameterBins, sizeDistribution[, asDict=False])

   Returns Mie coefficients β\ :sub:`ext`, β\ :sub:`sca`, β\ :sub:`abs`, G, β\ :sub:`pr`, β\ :sub:`back`, β\ :sub:`ratio`.
   
   **Parameters**
   
   
   m : complex
	The complex refractive index, with the convention *m = n+ik*.
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
	The complex refractive index, with the convention *m = n+ik*.
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
		


Angular Functions
-----------------

These functions compute the angle-dependent scattered field intensities, scattering matrix elements, and create arrays that are useful for plotting.

.. py:Function:: ScatteringFunction(m, wavelength, diameter[, minAngle=0, maxAngle=180, angularResolution=0.5, normed=False])

   Creates arrays for plotting the angular scattering intensity functions in theta-space with parallel, perpendicular, and unpolarized light. Uses :py:func:`MieS1S2` to compute S\ :sub:`1` and S\ :sub:`2`, then computes parallel, perpendicular, and unpolarized intensities by
   
		:math:`SR(\theta)=|S_1|^2`
		
		:math:`SL(\theta)=|S_2|^2`
		
		:math:`SU(\theta)=\frac{1}{2}(SR+SL)`
   
   **Parameters**
   
   
   m : complex
	The complex refractive index with the convention *m = n+ik*.
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


.. py:Function:: qSpaceScatteringFunction(m, wavelength, diameter[, normed=False])

   Creates arrays for plotting the angular scattering intensity functions in q-space with parallel, perpendicular, and unpolarized light. Uses :py:func:`MieS1S2` to compute S\ :sub:`1` and S\ :sub:`2`. The scattering angle variable, *qR*, is calculated by :math:`qR=(4\pi /\lambda)\,sin(\theta /2)\,(d_p /2)`.
   
   **Parameters**
   
   
   m : complex
	The complex refractive index with the convention *m = n+ik*.
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
	
	
.. py:Function:: MatrixElements(m, wavelength, diameter, mu)

   Calculate the four nonzero scattering matrix elements S\ :sub:`11`, S\ :sub:`12`, S\ :sub:`33`, and S\ :sub:`34` as functions of *μ*\ =cos(*θ*\ ), where *θ* is the scattering angle:
   
		:math:`S_{11}=\frac{1}{2}\left(|S_2|^2+|S_1|^2\right)`
		
		:math:`S_{12}=\frac{1}{2}\left(|S_2|^2-|S_1|^2\right)`
		
		:math:`S_{33}=\frac{1}{2}(S_2^*S_1^*+S_2S_1^*)`
		
		:math:`S_{34}=\frac{i}{2}(S_1S_2^*-S_2S_1^*)`
		
		
   **Parameters**
   
   
   m : complex
	The complex refractive index with the convention *m = n+ik*.
   wavelength : float
	The wavelength of incident light, in nanometers.
   diameter : float
	The diameter of the particle, in nanometers.
   mu : float
	The cosine of the scattering angle.

   **Returns**
   S11, S12, S33, S34 : float
	The matrix elements described above.
	

	
.. py:Function:: MieS1S2(m,x,mu)

   Calculates S\ :sub:`1` and S\ :sub:`2` at μ=cos(θ), where θ is the scattering angle. Must be explicitly imported via: ::
   
   >>> from PyMieScatt.Mie import MieS1S2
   
   Uses :py:func:`Mie_ab` to calculate a\ :sub:`n` and b\ :sub:`n`, and :py:func:`MiePiTau` to calculate π\ :sub:`n` and τ\ :sub:`n`. S\ :sub:`1` and S\ :sub:`2` are calculated by:
   
		:math:`S_1=\sum\limits_{n=1}^{n_{max}}\frac{2n+1}{n(n+1)}(a_n\pi_n+b_n\tau_n)`
		
		:math:`S_2=\sum\limits_{n=1}^{n_{max}}\frac{2n+1}{n(n+1)}(a_n\tau_n+b_n\pi_n)`
		
   **Parameters**
   
   
   m : complex
	The complex refractive index with the convention *m = n+ik*.
   x : float
	The size parameter :math:`x=\pi\,d_p/\lambda`.
   mu : float
	The cosine of the scattering angle.
	
   **Returns**
   
   
   S1, S2 : complex
	The S\ :sub:`1` and S\ :sub:`2` values.

.. py:Function:: MiePiTau(mu,nmax)

   Calculates π\ :sub:`n` and τ\ :sub:`n`. Must be explicitly imported via: ::
   
   >>> from PyMieScatt.Mie import MiePiTau
   
   This function uses recurrence relations to calculate π\ :sub:`n` and τ\ :sub:`n`, beginning with π\ :sub:`0` = 1, π\ :sub:`1` = 3μ (where μ is the cosine of the scattering angle), τ\ :sub:`0` = μ, and τ\ :sub:`1` = 3cos(2cos\ :sup:`-1` (μ)):
   
		:math:`\pi_n=\frac{2n-1}{n-1}\mu\pi_{n-1}-\frac{n}{n-1}\pi_{n-2}`
		
		:math:`\tau_n=n\mu\pi_n-(n+1)\pi_{n-1}`
		
   **Parameters**
   
   
   mu : float
	The cosine of the scattering angle.
   nmax : int
	The number of elements to compute, where n\ :sub:`max` = 2+x+4x\ :sup:`1/3`.
   **Returns**
   
   
   p, t : numpy.ndarray
	The π\ :sub:`n` and τ\ :sub:`n` arrays, of length *nmax*.