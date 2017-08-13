# libsFoamAux: Auxiliary libraries for OpenFOAM.

Notes:  Patches to OF source no longer needed (after commit c084134)  
        Source common to official and extended code streams  
        Forks to accomodate extended source will be in developPlus  
        Script mkInclude makes symbolic links  
        

steadyStateControl: User library. Based on the simpleControl class
($FOAM_SRC/finiteVolume/cfdTools/general/solutionControl/simpleControl/)

steadyStateControl iterates to convergence without updating the timestamp or writing intermediate data files. It is useful in simulations where a field (porosity for example) is evolving slowly and flow (calculated by the SIMPLE method) is assumed to be steady at each instance of the field. The runTime database then refers to the slowly varying field, rather than the flow field.

Files: steadyStateControl.H
       steadyStateControl.C

dissolMeshRelax: User library for mesh relaxation

Files: dissolMotionPointPatchVectorField.H
       dissolMotionPointPatchVectorField.C

coupledPatchInterpolation: User library for face to point interpolation
over a coupled patch

Files: coupledPatchInterpolation.H
       CoupledPatchInterpolation.H
       CoupledPatchInterpolation.C

OFstreamMod: User library adds support for appending to an existing file

Files: OFstreamMod.H
       OFstreamMod.C

boundaryConditions: Templates for coded boundary conditions

Files: nonlinear.H
       danckwerts.H


