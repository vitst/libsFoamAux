# libsFoamAux: Auxiliary libraries for OpenFOAM.

This is a branch of the v2012 to make it compartable with OF-dev (July 2021)

## Notes:  
* Source codes built with OpenFOAM-v2012(R)
* Script `mkInclude` makes symbolic links  
* The user must set the path to the `codeTemplates` directory (in `libsFoamAux`): e.g. `export FOAM_CODE_TEMPLATES=path-to-libsFoamAux/codeTemplates`
        
## steadyStateControl: 
User library. Based on the simpleControl class, `steadyStateControl` iterates to convergence without updating the timestamp or writing intermediate data files. It is useful in simulations where a field (porosity for example) is evolving slowly and flow (calculated by the `SIMPLE` method) is assumed to be steady at each instance of the field. The `runTime` database then refers to the slowly varying field, rather than the flow field.

## normalMotionSlip: 
User library for mesh relaxation. Templated from OF slip class with added normal velocity. Needs a coded include file (roughMotion.H or gradCMotion.H).

## boundaryConditions: 
Templates for coded boundary conditions.

## coupledPatchInterpolation:
User library for face to point interpolation over a coupled patch.

## OFstreamMod:
User library adds support for appending to an existing file.

## dissolMotion:
User library for mesh relaxation. Superseded by normalMotionSlip, but still used by surfRoughGen utility.

## velocityDeltatLaplacianFvMotionSolver:
This library is a modified `velocityLaplacianFvMotionSolver` and is used for dynamic mesh motion. The difference include:
* point displacement is dt*`pointMotionU`
* laplacian solver iterates if `iterateLaplacianSolve` in `dynamicMeshDict` is set to true
* laplacian solver iterates until `tolerance` specified in `dynamicMeshDict` is reached

## functionObjects:
This directory containes functionObjects. Thet mostly are postprocessing utilities which can be called using `postProcess` tool.
