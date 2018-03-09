/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "velocityDeltatLaplacianFvMotionSolver.H"
#include "motionInterpolation.H"
#include "motionDiffusivity.H"
#include "fvmLaplacian.H"
#include "fvcLaplacian.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(velocityDeltatLaplacianFvMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        velocityDeltatLaplacianFvMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::velocityDeltatLaplacianFvMotionSolver::velocityDeltatLaplacianFvMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    velocityMotionSolver(mesh, dict, typeName),
    fvMotionSolver(mesh),
    cellMotionU_
    (
        IOobject
        (
            "cellMotionU",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvMesh_,
        dimensionedVector
        (
            "cellMotionU",
            pointMotionU_.dimensions(),
            Zero
        ),
        cellMotionBoundaryTypes<vector>(pointMotionU_.boundaryField())
    ),
    interpolationPtr_
    (
        coeffDict().found("interpolation")
      ? motionInterpolation::New(fvMesh_, coeffDict().lookup("interpolation"))
      : motionInterpolation::New(fvMesh_)
    ),
    diffusivityPtr_
    (
        motionDiffusivity::New(fvMesh_, coeffDict().lookup("diffusivity"))
    )
{
  if( !coeffDict().readIfPresent<scalar>("tolerance", tolerance) ){
    SeriousErrorIn("velocityDeltatLaplacianFvMotionSolver.C")
        << "No 'tolerance' parameter in "<< coeffDict().name()
        << exit(FatalError);
  }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::velocityDeltatLaplacianFvMotionSolver::~velocityDeltatLaplacianFvMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::velocityDeltatLaplacianFvMotionSolver::curPoints() const
{
    interpolationPtr_->interpolate
    (
        cellMotionU_,
        pointMotionU_
    );
    
    tmp<pointField> tcurPoints
    (
        fvMesh_.points()
      + 
        pointMotionU_.primitiveField()
    );

    twoDCorrectPoints(tcurPoints.ref());

    return tcurPoints;
}


void Foam::velocityDeltatLaplacianFvMotionSolver::solve()
{
    // The points have moved so before interpolation update
    // the fvMotionSolver accordingly
    movePoints(fvMesh_.points());

    diffusivityPtr_->correct();
    pointMotionU_.boundaryFieldRef().updateCoeffs();
    
    int iter = 0;
    while ( true )
    {

      fvVectorMatrix UEqn
      (
          fvm::laplacian
          (
              diffusivityPtr_->operator()(),
              cellMotionU_,
              "laplacian(diffusivity,cellMotionD)"
          )
      );
      
      SolverPerformance<vector> sp 
              = UEqn.solveSegregatedOrCoupled(UEqn.solverDict());
      scalar residual = cmptMax(sp.initialResidual());

      if( residual < tolerance )
      {
          Info << "velocity laplacian: Converged in " 
               << iter << " steps.  Residual = "
               << residual << nl << endl;
          break;
      }
      if(debug)
      {
          Info << " Step " << iter << token::TAB
               << " residual: "<< residual << " > " << tolerance << endl;
      }

      iter++;
    }
}


void Foam::velocityDeltatLaplacianFvMotionSolver::updateMesh
(
    const mapPolyMesh& mpm
)
{
    velocityMotionSolver::updateMesh(mpm);

    // Update diffusivity. Note two stage to make sure old one is de-registered
    // before creating/registering new one.
    diffusivityPtr_.reset(nullptr);
    diffusivityPtr_ = motionDiffusivity::New
    (
        fvMesh_,
        coeffDict().lookup("diffusivity")
    );
}


// ************************************************************************* //
