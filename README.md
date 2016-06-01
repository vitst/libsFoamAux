# libsOpenFOAMauxiliary
Auxiliary external libraries for OpenFOAM.

* danckwerts: User library.


* dissolMeshRelax: User library.


* nonLinearFv: User library.


* steadyStateControl: User library. Based on the simpleControl class
($FOAM_SRC/finiteVolume/cfdTools/general/solutionControl/simpleControl/)

steadyStateControl iterates to convergence without updating the timestamp or writing intermediate data files. It is useful in simulations where a field (porosity for example) is evolving slowly and flow (calculated by the SIMPLE method) is assumed to be steady at each instance of the field. The runTime database then refers to the slowly varying field, rather than the flow field.

steadyStateControl requires dataDictPatch to be installed, which can be found in https://github.com/vitst/patches-OpenFOAM-2.x.x-unofficial.git

Files: steadyStateControl.H
       steadyStateControl.C

