.. highlight:: python

Example Scripts
===============

Mie Efficiencies of a Single Homogeneous Particle
-------------------------------------------------

To calculate the efficencies of a single homogeneous particle, use the :py:func:`MieQ` function. ::

	>>> import PyMieScatt as ps
	>>> ps.MieQ(1.5+0.5j,532,200,asDict=True)
	{'Qabs': 1.2206932456722366,
	 'Qback': 0.2557593071989655,
	 'Qext': 1.6932375984850729,
	 'Qpr': 1.5442174328606282,
	 'Qratio': 0.5412387338385265,
	 'Qsca': 0.47254435281283641,
	 'g': 0.3153569918620277}


Mie Efficencies of a Weibull Distribution
-----------------------------------------

Consider the 405 nm Mie coefficients of 105 particles/cm3, with m = 1.5+0.5i, in a `Weibull distribution <https://en.wikipedia.org/wiki/Weibull_distribution>`_ with shape parameter sh = 5 and scale parameter sc = 200. This would take seven lines in Python, three of which are module import commands: ::

	>>> import PyMieScatt as ps
	>>> import numpy as np
	>>> import matplotlib.pyplot as plt
	>>> dp = np.linspace(10,1000,1000)
	>>> N,sh,sc = 1e5,5,200
	>>> w=[N*((sh/sc)*(d/sc)**(sh-1))*np.exp(-(d/sc)**sh) for d in dp]
	>>> ps.Mie_withSizeDistribution(1.5+0.5j,405,dp,w,asDict=True)
	{'Babs': 3762.0479602613427,
	 'Bback': 286.65698999981691,
	 'Bext': 5747.4466502095638,
	 'Bpr': 4662.181554274106,
	 'Bratio': 550.87163111634698,
	 'Bsca': 1985.3986899482211,
	 'G': 0.54662325578736115}
