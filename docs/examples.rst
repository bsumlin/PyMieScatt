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