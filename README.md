sage-foliations
===============

GENERAL DESCRIPTION:

	A package whose primary function is the construction of pseudo-Anosov maps from measured foliations on 
	surfaces using a method similar to that of Arnoux and Yoccoz. It also handles canonical train tracks carrying
	the measured foliations, and maps between train tracks. Given a pseudo-Anosov map, the Teichmuller 
	polynomial of the corresponding fibered face of its mapping torus can be calculated (even in the 
	non-orintable case), thus helping the construction of even more pseudo-Anosovs by other fibrations of 
	the 3-manifold.

	The idea of the pseudo-Anosov construction is as follows: 
	1. Pick a random measured foliation on any (possibly non-orientable or punctured) surface with negative
	Euler characteristic, defined by an interval exchange map on a simple closed curve on the surface.
	2. Find another simple closed curve transverse to the foliation.
	3. Calculate the canonical train tracks carrying the foliation with respect to the two curves. 
	The first train track always carries the second one, so a train track map can be calculated between the two.
	4. If the two train tracks are isomorphic and the train track map is "mixing" enough, then the eigendata 
	of the transition map defines a measure on the train tracks which in turn define a "corrected" measured 
	foliation (with the same combinatorial data, but possibly different length parameters) which is invariant 
	under a pseudo-Anosov map.

	The Teichmuller polynomial then can be calculated from train track map as described by McMullen.

USE:

	Right now the easiest way to learn how to use the package is to contact me (strenner@math.wisc.edu) and 
	have me explain what is possible to do with it. I am planning to make it more user-friendly over time, 
	but it is a slow process.

TODO:

	1. A large part of the code needs polishing and documentation. Especially for adding this code to Sage 
	one day, as planned originally.
	2. Other useful features could be added, like the calculation of the Alexander polynomial.
	3. One major, but very interesting project would be to make the code capable of finding the other 
	invariant foliation of pseudo-Anosovs as well, thus constructing flat surfaces. The motivation for this 
	is that one could hopefully construct exotic flat surfaces, like the Arnoux-Yoccoz surfaces that doesn't 
	have parabolic elements in the Veech group.

