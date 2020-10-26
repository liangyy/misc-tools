Here I provide a light weight variance component solver on the basis of [EMMA](http://mouse.cs.ucla.edu/emma/).
Currently, I implemented both MLE and REML based estimations. 
The solver uses Python numpy as backend.

It provides the following minimal functionality:

* Load GRM in GRM.GZ format (the gcta format).
* Perform EVD on the GRM and cache the result.
* Perform variance component analysis using GRM EVD and a matrix of observed y's.


I provide a prototype script in R (`emma_solver.R`) as the side note of the implementation.
