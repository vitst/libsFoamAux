# libsFoamAux: Auxiliary libraries for OpenFOAM.

Notes:  Source codes built with OpenFOAM-v1706; v1712 should also work
        Script mkInclude makes symbolic links  
        The user must set the path to the codeTemplates directory (in libsFoamAux): e.g.
        export FOAM_CODE_TEMPLATES=path-to-libsFoamAux/codeTemplates
        

steadyStateControl: User library. Based on the simpleControl class,
steadyStateControl iterates to convergence without updating the timestamp or writing intermediate data files. It is useful in simulations where a field (porosity for example) is evolving slowly and flow (calculated by the SIMPLE method) is assumed to be steady at each instance of the field. The runTime database then refers to the slowly varying field, rather than the flow field.

Files: steadyStateControl.H
       steadyStateControl.C

normalMotionSlip: User library for mesh relaxation. Templated from OF slip class with added normal velocity. Needs a coded include file (roughMotion.H or gradCMotion.H)

Files: normalMotionSlipFvPatchField/normalMotionSlipFvPatchField.H
       normalMotionSlipFvPatchField/normalMotionSlipFvPatchField.C
       normalMotionSlipPointPatchField/normalMotionSlipPointPatchField.H
       normalMotionSlipPointPatchField/normalMotionSlipPointPatchField.C

boundaryConditions: Templates for coded boundary conditions

Files: linear.H                     // linear kinetics boundary condition
       nonlinear.H                  // nonlinear kinetics boundary condition
       danckwerts.H                 // Danckwerts boundary condition
       roughMotion.H                // pointMotionU boundary condition for rough surface
       gradCMotion.H                // pointMotionU boundary condition for surface concentration gradient

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

