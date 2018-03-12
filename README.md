# libsFoamAux: Auxiliary libraries for OpenFOAM.

Notes:  Source codes built with OpenFOAM-v1706; v1712 should also work
        Script mkInclude makes symbolic links  
        

steadyStateControl: User library. Based on the simpleControl class,
steadyStateControl iterates to convergence without updating the timestamp or writing intermediate data files. It is useful in simulations where a field (porosity for example) is evolving slowly and flow (calculated by the SIMPLE method) is assumed to be steady at each instance of the field. The runTime database then refers to the slowly varying field, rather than the flow field.

Files: steadyStateControl.H
       steadyStateControl.C

normalMotionSlip: User library for mesh relaxation. Templated from OF slip class with added normal velocity.

Files: normalMotionSlipFvPatchField/normalMotionSlipFvPatchField.H
       normalMotionSlipFvPatchField/normalMotionSlipFvPatchField.C
       normalMotionSlipPointPatchField/normalMotionSlipPointPatchField.H
       normalMotionSlipPointPatchField/normalMotionSlipPointPatchField.C

boundaryConditions: Templates for coded boundary conditions

Files: linear.H
       nonlinear.H
       danckwerts.H

coupledPatchInterpolation: User library for face to point interpolation
over a coupled patch

Files: coupledPatchInterpolation.H
       CoupledPatchInterpolation.H
       CoupledPatchInterpolation.C

OFstreamMod: User library adds support for appending to an existing file; used by dissolCalc

Files: OFstreamMod.H
       OFstreamMod.C

dissolMotion: User library for mesh relaxation. Superseded by normalMotionSlip, but still used by surfRoughGen.

Files: dissolMotionPointPatchVectorField.H
       dissolMotionPointPatchVectorField.C


