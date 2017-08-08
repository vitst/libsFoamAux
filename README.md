# libsFoamAux: Auxiliary libraries for OpenFOAM.

Notes:  Patches to OF source no longer needed (after commit c084134)
        Source common to official and extended code streams
        Forks to accomodate extended source will be in developPlus
        

* steadyStateControl: User library. Based on the simpleControl class
($FOAM_SRC/finiteVolume/cfdTools/general/solutionControl/simpleControl/)

steadyStateControl iterates to convergence without updating the timestamp or writing intermediate data files. It is useful in simulations where a field (porosity for example) is evolving slowly and flow (calculated by the SIMPLE method) is assumed to be steady at each instance of the field. The runTime database then refers to the slowly varying field, rather than the flow field.

Files: steadyStateControl.H
       steadyStateControl.C

* dissolMeshRelax: User library.


* nonLinearFv: User library. Replaced by coded bc (template in boundary conditions)


* danckwerts: User library. To be replaced by coded bc

