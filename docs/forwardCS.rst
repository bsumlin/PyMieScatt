Functions for Coated Spheres (Core-Shell Particles)
===================================================

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