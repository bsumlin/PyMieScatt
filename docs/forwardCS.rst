Functions for Coated Spheres (Core-Shell Particles)
===================================================

.. py:Function:: MieQCoreShell(mCore, mShell, wavelength, dCore, dShell[, asDict=False, asCrossSection=False])

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
   asCrossSection : bool, optional
	If specified and set to True, returns the results as optical cross-sections with units of nm\ :sup:`2`.
	
   **Returns**
   
   
   qext, qsca, qabs, g, qpr, qback, qratio : float
	The Mie efficencies described above.
   cext, csca, cabs, g, cpr, cback, cratio : float
	If asCrossSection==True, :py:func:`MieQCoreShell` returns optical cross-sections.
   q : dict
	If asDict==True, :py:func:`MieQCoreShell` returns a dict of the above efficiencies with appropriate keys.
   c : dict
	If asDict==True and asCrossSection==True, returns a dict of the above cross-sections with appropriate keys.
	
   **Considerations**
   
   
   When using this function in a script, there are three simplifying clauses that can speed up computation when considering both coated and homogeneous particles. Upon determining the size parameters of the core and the shell:
   
   - if x\ :sub:`core` == x\ :sub:`shell`, then :py:func:`MieQCoreShell` returns Mie efficencies calculated by MieQ(mCore,wavelength,dShell).
   - If x\ :sub:`core` == 0, then :py:func:`MieQCoreShell` returns efficencies calculated by MieQ(mShell,wavelength,dShell).
   - If m\ :sub:`core` == m\ :sub:`shell`, then :py:func:`MieQCoreShell` returns efficencies calculated by MieQ(mCore,wavelength,dShell).
   
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

.. py:Function:: CoreShellScatteringFunction(mCore, mShell, wavelength, dCore, dShell[, minAngle=0, maxAngle=180, angularResolution=0.5, normed=False])
   
   Computes the angle-dependent scattering intensity of a coated sphere.
   
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
   thetaSteps : int
	The number of points between 0 and 180 degrees to use in calculations.

   **Returns**
   
   
   theta : numpy.ndarray
	An array of the angles used in calculations. Values will be spaced according to *angularResolution*, and the size of the array will be *(maxAngle-minAngle)/angularResolution*.
   SL : numpy.ndarray
	An array of the scattered intensity of left-polarized (parallel) light. Same size as the *theta* array.
   SR : numpy.ndarray
	An array of the scattered intensity of right-polarized (perpendicular) light. Same size as the *theta* array.
   SU : numpy.ndarray
	An array of the scattered intensity of unpolarized light, which is the average of SL and SR. Same size as the *theta* array.
   
.. py:Function:: CoreShellS1S2(mCore, mShell, xCore, xShell, mu)

   Computes S1 and S2 of a coated sphere as a function of mu, the cosine of the scattering angle.
   
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
   mu : float
	The cosine of the scattering angle.
	
   **Returns**
   
   
   S1, S2 : complex
	The S\ :sub:`1` and S\ :sub:`2` values.
	
	
.. py:Function:: CoreShellMatrixElements(mCore, mShell, xCore, xShell, mu)

   Calculates the four nonzero scattering matrix elements S\ :sub:`11`, S\ :sub:`12`, S\ :sub:`33`, and S\ :sub:`34` as functions of *μ*\ =cos(*θ*\ ), where *θ* is the scattering angle:
   
		:math:`S_{11}=\frac{1}{2}\left(|S_2|^2+|S_1|^2\right)`
		
		:math:`S_{12}=\frac{1}{2}\left(|S_2|^2-|S_1|^2\right)`
		
		:math:`S_{33}=\frac{1}{2}(S_2^*S_1^*+S_2S_1^*)`
		
		:math:`S_{34}=\frac{i}{2}(S_1S_2^*-S_2S_1^*)`
		
		
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
   mu : float
	The cosine of the scattering angle.

   **Returns**
   S11, S12, S33, S34 : float
	The matrix elements described above.