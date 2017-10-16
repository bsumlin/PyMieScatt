Inverse Mie Theory Functions
============================

Contour Intersection Inversion Functions
----------------------------------------

For more details on the contour intersection inversion method, please see Sumlin BJ, Heinson WR, Chakrabarty RK. *Retrieving the Aerosol Complex Refractive Index using PyMieScatt: A Mie Computational Package with Visualization Capabilities*. J. Quant. Spectros. Rad. Trans. 2017. DOI: `10.1016/j.jqsrt.2017.10.012 <https://doi.org/10.1016/j.jqsrt.2017.10.012>`_

.. py:function:: ContourIntersection(Qsca, Qabs, wavelength, diameter[, nMin=1, nMax=3, kMin=0.00001, kMax=1, Qback=None, gridPoints=100, interpolationFactor=2, maxError=0.005, fig=None, ax=None, axisOption=0])

   Computes complex *m = n+ik* from a particle diameter (in nm), incident wavelength (in nm), and scattering and absorption efficiencies. Optionally, backscatter efficiency may be specified to constrain the problem to produce a unique solution.
   
   **Parameters**
   
   
   Qsca : float or list-like
	The scattering efficiency, or optionally, a list, tuple, or numpy.ndarray of scattering efficiency and its associated error.
   Qabs : float or list-like
	The absorption efficiency, or optionally, a list, tuple, or numpy.ndarray of absorption efficiency and its associated error..
   wavelength : float
	The wavelength of incident light, in nm.
   diameter : float
	The diameter of the particle, in nm.
   nMin : float, optional
	The minimum value of *n* to search.
   nMax : float, optional
	The maximum value of *n* to search.
   kMin : float, optional
	The minimum value of *k* to search.
   kMax : float, optional
	The maximum value of *k* to search.
   Qback : float or list-like, optional
	The backscatter efficiency, or optionally, a list, tuple, or numpy.ndarray of backscatter efficiency and its associated error.
   gridPoints : int, optional
	The number of gridpoints for the search mesh. Defaults to 200. Increase for better resolution but longer run times.
   interpolationFactor : int, optional
	The interpolation to apply to the search fields, artificially increasing their resolutions. This is applied after calculations, so some features may be lost if **interpolationFactor** is too high and **gridPoints** is too low.
   maxError : float, optional
	The allowed error in forward calculations of the retrived *m*.
   fig : matplotlib.figure object, optional (but recommended)
	The figure object to send to the geometric inversion routine. If unspecified, one will be created.
   ax : matplotlib.axes object, optional (but recommended)
	The axes object to send to the geometric inversion routine. If unspecified, one will be created.

   **Returns**
   
   
   solutionSet : list
	A list of all valid solutions
   ForwardCalculations : list
	A list of scattering and absorption efficencies produced by forward Mie calculations using the derived refractive indices
   solutionErrors : list
	The relative errors of the efficencies in **ForwardCalculations**.
   fig : matplotlib.figure object
	The figure object now associated with the inversion calculations.
   ax : matplotlib.axes object
	The axes object now associated with the inversion calculations.
   graphElements : dict
	A dict of all artists necessary to fully manipulate the appearance of the output. The keys will depend on the options passed to the inversion function itself (i.e., errors specified, backscatter specified). Maximally, it will contain:
	'Qsca', 'Qabs', 'Qback' : the major contours;
	'QscaErrFill', 'QscaErrOutline1', 'QscaErrOutline2' : the error bound outlines;
	'QabsErrFill', 'QabsErrOutline1', 'QabsErrOutline2' : the error bound fills;
	'SolMark', 'SolFill' : the circly thingies at each solution;
	'CrosshairsH', 'CrosshairsV' : solution crosshairs;
	'LeftSpine', 'RightSpine', 'BottomSpine', 'TopSpine' : graph spines;
	'XAxis', 'YAxis' : the individual matplotlib axis objects.
	