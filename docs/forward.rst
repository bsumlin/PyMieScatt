Functions for Forward Mie Calculations of Homogeneous Spheres
=============================================================

Functions for single particles
---------------------------------

.. py:function:: MieQ(m, wavelength, diameter[, asDict=False])

   Computes Mie efficencies *Q* and asymmetry parameter *g* of a single, homogeneous particle. Uses :py:func:`Mie_ab` to calculate :math:`a_n` and :math:`b_n`, and then calculates *Q* via:
   
		:math:`${\displaystyle Q_{ext}=\frac{2}{x^2}\sum_{n=1}^{n_{max}}(2n+1)\:\text{Re}\left\{a_n+b_n\right\}}$`
		
		:math:`${\displaystyle Q_{sca}=\frac{2}{x^2}\sum_{n=1}^{n_{max}}(2n+1)(|a_n|^2+|b_n|^2)}$`
		
		:math:`${\displaystyle Q_{abs}=Q_{ext}-Q_{sca}}$`
		
		:math:`${\displaystyle Q_{back}=\frac{1}{x^2} \left| \sum_{n=1}^{n_{max}}(2n+1)(-1)^n(a_n-b_n) \right| ^2}$`
		
		:math:`${\displaystyle Q_{ratio}=\frac{Q_{back}}{Q_{sca}}}$`
		
		:math:`${\displaystyle g=\frac{4}{Q_{sca}x^2}\left[\sum\limits_{n=1}^{n_{max}}\frac{n(n+2)}{n+1}\text{Re}\left\{a_n a_{n+1}^*+b_n b_{n+1}^*\right\}+\sum\limits_{n=1}^{n_{max}}\frac{2n+1}{n(n+1)}\text{Re}\left\{a_n b_n^*\right\}\right]}$`
		
		:math:`${\displaystyle Q_{pr}=Q_{ext}-gQ_{sca}}$`
		
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
   
		:math:`${\displaystyle Q_{sca}=\frac{8x^4}{3}\left|{\frac{m^2-1}{m^2+2}}\right|^2}$`
   
		:math:`${\displaystyle Q_{abs}=4x\:\text{Im}\left\{\frac{m^2-1}{m^2+2}\right\}}$`
   
		:math:`${\displaystyle Q_{ext}=Q_{sca}+Q_{abs}}$`
   
		:math:`${\displaystyle Q_{back}=\frac{3Q_{sca}}{2}}$`
   
		:math:`${\displaystyle Q_{ratio}=1.5}$`
   
		:math:`${\displaystyle Q_{pr}=Q_{ext}}$`
		
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
   
   
.. py:Function:: AutoMieQ(m, wavelength, diameter[, crossover=0.5, asDict=False])

   Returns Mie efficencies of a spherical particle according to either :py:func:`MieQ` or :py:func:`RayleighMieQ` depending on the magnitude of the size parameter. Good for studying parameter ranges or size distributions.
   
   **Parameters**
   
   m : complex
	The complex refractive index, with the convention *m = n+ik*.
   wavelength : float
	The wavelength of incident light, in nanometers.
   diameter : float
	The diameter of the particle, in nanometers.
   crossover : float, optional
	The size parameter that dictates where calculations switch from Rayleigh approximation to actual Mie.
   asDict : bool, optional
	If specified and set to True, returns the results as a dict.
	
   **Returns**
   
   
   qext, qsca, qabs, g, qpr, qback, qratio : float
	The Mie efficencies described above.
   q : dict
	If asDict==True, :py:func:`RayleighMieQ` returns a dict of the above values with appropriate keys.


.. py:Function:: LowFrequencyMieQ(m, wavelength, diameter[, asDict=False])

   Returns Mie efficencies of a spherical particle in the low-frequency regime (:math:`x=\pi\,d_p/\lambda \ll 1`) given refractive index **m**, **wavelength**, and **diameter**. Optionally returns the parameters as a dict when **asDict** is specified and set to True. Uses :py:func:`LowFrequencyMie_ab` to calculate a\ :sub:`n` and b\ :sub:`n`, and follows the same math as :py:func:`MieQ`.
   
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
   
   For example, compute the Mie efficencies of a particle 100 nm in diameter with m = 1.33+0.01i, illuminated by λ = 1600 nm: ::
   
		>>> import PyMieScatt as ps
		>>> ps.LowFrequencyMieQ(1.33+0.01j,1600,100,asDict=True)
		{'Qabs': 0.0044765816617916582,
		 'Qback': 0.00024275862007727458,
		 'Qext': 0.0046412326004135135,
		 'Qpr': 0.0046400675577583459,
		 'Qratio': 1.4743834569616665,
		 'Qsca': 0.00016465093862185558,
		 'g': 0.0070758336692078412}

.. py:Function:: LowFrequencyMie_ab(m,x)

   Returns external field coefficients a\ :sub:`n` and b\ :sub:`n` based on inputs of **m** and :math:`x=\pi\,d_p/\lambda` by limiting the expansion of a\ :sub:`n` and b\ :sub:`n` to second order:
   
		:math:`${\displaystyle a_1=\frac{m^2-1}{m^2+2} \left[ -\frac{i2x^3}{3}-\frac{2ix^5}{5}\left( \frac{m^2-2}{m^2+2}\right) +\frac{4x^6}{9}\left( \frac{m^2-1}{m^2+2} \right) \right]}$`
   
		:math:`${\displaystyle a_2=-\frac{ix^5}{15}\frac{(m^2-1)}{2m^2+3}}$`
   
		:math:`${\displaystyle b_1=-\frac{ix^5}{45}(m^2-1)}$`
   
		:math:`${\displaystyle b_2=0}$`
		
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

.. py:Function:: MieQ_withDiameterRange(m, wavelength[, diameterRange=(10,1000), nd=1000, logD=False])

   Computes the Mie efficencies of particles across a diameter range using :py:func:`AutoMieQ`.
   
   **Parameters**
   
   
   m : complex
	The complex refractive index with the convention *m = n+ik*.
   wavelength : float
	The wavelength of incident light, in nanomaters
   diameterRange : tuple or list, optional
	The diameter range, in nanometers. Convention is (smallest, largest). Defaults to (10, 1000).
   nd : int, optional
	The number of diameter bins in the range. Defaults to 1000.
   logD : bool, optional
	If True, will use logarithmically-spaced diameter bins. Defaults to False.
	
   **Returns**
   
   
   diameters : numpy.ndarray
	An array of the diameter bins that calculations were performed on. Size is equal to **nd**.
   qext, qsca, qabs, g, qpr, qback, qratio : numpy.ndarray
	The Mie efficencies at each diameter in **diameters**.
	
.. py:Function:: MieQ_withWavelengthRange(m, diameter[, wavelengthRange=(100,1600), nw=1000, logW=False])

   Computes the Mie efficencies of particles across a wavelength range using :py:func:`AutoMieQ`. This function can optionally take a list, tuple, or numpy.ndarray for **m**. If your particles have a wavelength-dependent refractive index, you can study it by specifying **m** as list-like. When doing so, **m** must be the same size as **wavelengthRange**, which is also specified as list-like in this situation. Otherwise, the function will construct a range from **wavelengthRange[0]** to **wavelengthRange[1]** with **nw** entries.
   
   **Parameters**
   
   
   m : complex or list-like
	The complex refractive index with the convention *m = n+ik*. If dealing with a dispersive material, then len(**m**) must be equal to len(**wavelengthRange**).
   diameter : float
	The diameter of the particle, in nanometers.
   wavelengthRange : tuple or list, optional
	The wavelength range of incident light, in nanomaters. Convention is (smallest, largest). Defaults to (100, 1600). When **m** is list-like, len(**wavelengthRange**) must be equal to len(**m**).
   nw : int, optional
	The number of wavelength bins in the range. Defaults to 1000. This parameter is ignored if **m** is list-like.
   logW : bool, optional
	If True, will use logarithmically-spaced wavelength bins. Defaults to False. This parameter is ignored if **m** is list-like.
	
   **Returns**
   
   
   wavelengths : numpy.ndarray
	An array of the wavelength bins that calculations were performed on. Size is equal to **nw**, unless **m** was list-like. Then **wavelengths** = **wavelengthRange**.
   qext, qsca, qabs, g, qpr, qback, qratio : numpy.ndarray
	The Mie efficencies at each wavelength in **wavelengths**.
	
.. py:Function:: MieQ_withSizeParameterRange(m[, xRange=(1,10), nx=1000, logX=False])

   Computes the Mie efficencies of particles across a size parameter range (\ :math:`x=\pi\,d_p/\lambda`\ ) using :py:func:`AutoMieQ`.
   
   **Parameters**
   
   
   m : complex
	The complex refractive index with the convention *m = n+ik*.
   xRange : tuple or list, optional
	The size parameter range. Convention is (smallest, largest). Defaults to (1, 10).
   nx : int, optional
	The number of size parameter bins in the range. Defaults to 1000.
   logX : bool, optional
	If True, will use logarithmically-spaced size parameter bins. Defaults to False.
	
   **Returns**
   
   
   xValues : numpy.ndarray
	An array of the size parameter bins that calculations were performed on. Size is equal to **nx**.
   qext, qsca, qabs, g, qpr, qback, qratio : numpy.ndarray
	The Mie efficencies at each size parameter in **xValues**.


Functions for polydisperse size distributions of homogeneous spheres
--------------------------------------------------------------------

When an efficiency *Q* is integrated over a size distribution n\ :sub:`d`\ (d\ :sub:`p`), the result is the *coefficient* :math:`\beta`, which carries units of inverse length. The general form is:

		:math:`${\displaystyle \beta=10^{-6} \int\limits_{0}^{\infty}\frac{\pi d_p^2}{4}Q(m,\lambda,d_p)n(d_p)dd_p}$`
		
where d\ :sub:`p` is the diameter of the particle (in nm), n(d\ :sub:`p`) is the number of particles of diameter d\ :sub:`p` (per cubic centimeter), and the factor 10\ :sup:`-6` is used to cast the result in units of Mm\ :sup:`-1`. 

The bulk asymmetry parameter *G* is calculated by:

		:math:`${\displaystyle G=\frac{\int g(d_p)\beta_{sca}(d_p)dd_p}{\int \beta_{sca}(d_p)dd_p}}$`
		

.. py:Function:: MieQ_withSizeDistribution(m, wavelength, sizeDistributionDiameterBins, sizeDistribution[, asDict=False])

   Returns Mie coefficients β\ :sub:`ext`, β\ :sub:`sca`, β\ :sub:`abs`, G, β\ :sub:`pr`, β\ :sub:`back`, β\ :sub:`ratio`. Uses `scipy.integrate.trapz <https://docs.scipy.org/doc/scipy-0.10.1/reference/generated/scipy.integrate.trapz.html>`_ to compute the integral, which can introduce errors if your distribution is too sparse. Best used with a continuous, compactly-supported distribution.
   
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
	The Mie coefficients calculated by :py:func:`AutoMieQ`, integrated over the size distribution.
   q : dict
	If asDict==True, :py:func:`MieQ_withSizeDistribution` returns a dict of the above values with appropriate keys.

.. py:Function:: Mie_Lognormal(m, wavelength, geoStdDev, geoMean, numberOfParticles[, numberOfBins=1000, lower=1, upper=1000, gamma=[1], returnDistribution=False, decomposeMultimodal=False, asDict=False])

   Returns Mie coefficients :math:`\beta_{ext}`, :math:`\beta_{sca}`, :math:`\beta_{abs}`, :math:`G`, :math:`\beta_{pr}`, :math:`\beta_{back}`,  and :math:`\beta_{ratio}`, integrated over a mathematically-generated k-modal lognormal particle number distribution. Uses `scipy.integrate.trapz <https://docs.scipy.org/doc/scipy-0.10.1/reference/generated/scipy.integrate.trapz.html>`_ to compute the integral.
   
   The general form of a k-modal lognormal distribution is given by:
   
		:math:`${\displaystyle n(d_p)=\frac{N_\infty}{\sqrt{2\pi}} \sum_{i}^{k}\frac{\gamma_i}{d_p\ln\sigma_{g_i}}\exp\left\{ \frac{-(\ln d_p-\ln d_{pg_i})^2}{2 \ln^2\sigma_{g_i}}\right\}}$`
		
   where :math:`d_{p}` is the diameter of the particle (in nm), :math:`n(d_{p})` is the number of particles of diameter :math:`d_{p}` (per cubic centimeter), :math:`N_\infty` is the total number of particles in the distribution, :math:`\sigma_{g_i}` is the geometric standard deviation of mode :math:`i`, and :math:`d_{pg_i}` is the geometric mean diameter (in nm) of the *i*\ :sup:`th` moment. :math:`\gamma_i` is a porportionality constant that determines the fraction of total particles in the *i*\ :sup:`th` moment.
   
   This function is essentially a wrapper for :py:func:`Mie_withSizeDistribution`. A warning will be raised if the distribution is not compactly-supported on the interval specified by **lower** and **upper**.
   
   
   **Parameters**
   
   
   m : complex
	The complex refractive index, with the convention *m = n+ik*.
   wavelength : float
	The wavelength of incident light, in nanometers.
   geoStdDev : float or list-like
	The geometric standard deviation(s) :math:`\sigma_g` or :math:`\sigma_{g_i}` if list-like.
   geoMean : float or list-like
	The geometric mean diameter(s) :math:`d_{pg}` or :math:`d_{pg_i}` if list-like, in nanometers.
   numberOfParticles : float
	The total number of particles in the distribution.
   numberOfBins : int, optional
	The number of discrete bins in the distribution. Defaults to 1000.
   lower : float, optional
	The smallest diameter bin, in nanometers. Defaults to 1 nm.
   upper : float, optional
	The largest diameter bin, in nanometers. Defaults to 1000 nm.
   gamma : list-like, optional
	The porportionality coefficients for dividing total particles among modes.
   returnDistribution : bool, optional
	If True, both the size distribution bins and number concentrations will be returned.
   decomposeMultimodal: bool, optional
	If True (and returnDistribution==True), then the function returns an additional parameter containing the individual modes of the distribution.
   asDict : bool, optional
	If True, returns the results as a dict.
	
   **Returns**
   
   
   Bext, Bsca, Babs, G, Bpr, Bback, Bratio : float
	The Mie coefficients calculated by :py:func:`MieQ`, integrated over the size distribution.
   diameters, nd : numpy.ndarray
	The diameter bins and number concentrations per bin, respectively. Only if returnDistribution is True.
   ndi : list of numpy.ndarray objects
	A list whose entries are the individual modes that created the multimodal distribution. Only returned if both returnDistribution and decomposeMultimodal are True.
   B : dict
	If asDict==True, :py:func:`MieQ_withLognormalDistribution` returns a dict of the above values with appropriate keys.
   
   For example, compute the Mie coefficients of a lognormal size distribution with 1000000 particles, σ\ :sub:`g` = 1.7, and d\ :sub:`pg` = 200 nm; with m = 1.60+0.08i and λ = 532 nm: ::
   
		>>> import PyMieScatt as ps
		>>> ps.MieQ_Lognormal(1.60+0.08j,532,1.7,200,1e6,asDict=True)
		{'Babs': 33537.324569179938,
		'Bback': 10188.473118449627,
		'Bext': 123051.1109783932,
		'Bpr': 62038.347528346232,
		'Bratio': 12701.828124508347,
		'Bsca': 89513.786409213266,
		'bigG': 0.6816018615403715}
		
.. py:Function:: Mie_OtherDistribution(m, wavelength, numberOfParticles[, distribution='powerlaw', params=[numberOfBins=1000, lower=1, upper=1000, returnDistribution=False, decomposeMultimodal=False, asDict=False])

   Returns Mie coefficients :math:`\beta_{ext}`, :math:`\beta_{sca}`, :math:`\beta_{abs}`, :math:`G`, :math:`\beta_{pr}`, :math:`\beta_{back}`,  and :math:`\beta_{ratio}`, integrated over a mathematically-generated unimodal particle number distribution of a specified type. Uses `scipy.integrate.trapz <https://docs.scipy.org/doc/scipy-0.10.1/reference/generated/scipy.integrate.trapz.html>`_ to compute the integral.
   
   Not yet implemented.

Angular Functions
-----------------

These functions compute the angle-dependent scattered field intensities and scattering matrix elements. They return arrays that are useful for plotting.

.. py:Function:: ScatteringFunction(m, wavelength, diameter[, minAngle=0, maxAngle=180, angularResolution=0.5, space='theta', angleMeasure='radians', normed=False])

   Creates arrays for plotting the angular scattering intensity functions in theta-space with parallel, perpendicular, and unpolarized light. Also includes an array of the angles for each step. This angle can be in either degrees, radians, or gradians for some reason. The angles can either be geometrical angle or the qR vector (see `Sorensen, M. Q-space analysis of scattering by particles: a review. J. Quant. Spectrosc. Radiat. Transfer 2013, 131, 3-12 <http://www.sciencedirect.com/science/article/pii/S0022407313000083>`_). Uses :py:func:`MieS1S2` to compute S\ :sub:`1` and S\ :sub:`2`, then computes parallel, perpendicular, and unpolarized intensities by
   
		:math:`${\displaystyle SR(\theta)=|S_1|^2}$`
		
		:math:`${\displaystyle SL(\theta)=|S_2|^2}$`
		
		:math:`${\displaystyle SU(\theta)=\frac{1}{2}(SR+SL)}$`
   
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
   space : str, optional
	The measure of scattering angle. Can be 'theta' or 'qspace'. Defaults to 'theta'.
   angleMeasure : str, optional
	The units for the scattering angle
   normed : bool, optional
	If True, will normalize the output such that the maximum intensity will be 1.0. Defaults to False.
	
   **Returns**
   
   
   theta : numpy.ndarray
	An array of the angles used in calculations. Values will be spaced according to **angularResolution**, and the size of the array will be *(maxAngle-minAngle)/angularResolution*.
   SL : numpy.ndarray
	An array of the scattered intensity of left-polarized (parallel) light. Same size as the **theta** array.
   SR : numpy.ndarray
	An array of the scattered intensity of right-polarized (perpendicular) light. Same size as the **theta** array.
   SU : numpy.ndarray
	An array of the scattered intensity of unpolarized light, which is the average of SL and SR. Same size as the **theta** array.

	
	
.. py:Function:: MatrixElements(m, wavelength, diameter, mu)

   Calculates the four nonzero scattering matrix elements S\ :sub:`11`, S\ :sub:`12`, S\ :sub:`33`, and S\ :sub:`34` as functions of *μ*\ =cos(*θ*\ ), where *θ* is the scattering angle:
   
		:math:`${\displaystyle S_{11}=\frac{1}{2}\left(|S_2|^2+|S_1|^2\right)}$`
		
		:math:`${\displaystyle S_{12}=\frac{1}{2}\left(|S_2|^2-|S_1|^2\right)}$`
		
		:math:`${\displaystyle S_{33}=\frac{1}{2}(S_2^*S_1^*+S_2S_1^*)}$`
		
		:math:`${\displaystyle S_{34}=\frac{i}{2}(S_1S_2^*-S_2S_1^*)}$`
		
		
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
   
		:math:`${\displaystyle S_1=\sum\limits_{n=1}^{n_{max}}\frac{2n+1}{n(n+1)}(a_n\pi_n+b_n\tau_n)}$`
		
		:math:`${\displaystyle S_2=\sum\limits_{n=1}^{n_{max}}\frac{2n+1}{n(n+1)}(a_n\tau_n+b_n\pi_n)}$`
		
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
   
		:math:`${\displaystyle \pi_n=\frac{2n-1}{n-1}\mu\pi_{n-1}-\frac{n}{n-1}\pi_{n-2}}$`
		
		:math:`${\displaystyle \tau_n=n\mu\pi_n-(n+1)\pi_{n-1}}$`
		
   **Parameters**
   
   
   mu : float
	The cosine of the scattering angle.
   nmax : int
	The number of elements to compute. Typically, n\ :sub:`max` = floor(2+x+4x\ :sup:`1/3`\ ), but can be given any integer.
	
   **Returns**
   
   
   p, t : numpy.ndarray
	The π\ :sub:`n` and τ\ :sub:`n` arrays, of length **nmax**.