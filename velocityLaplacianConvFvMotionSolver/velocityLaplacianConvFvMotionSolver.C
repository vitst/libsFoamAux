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

#include "velocityLaplacianConvFvMotionSolver.H"
#include "motionInterpolation.H"
#include "motionDiffusivity.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(velocityLaplacianConvFvMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        velocityLaplacianConvFvMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::velocityLaplacianConvFvMotionSolver::velocityLaplacianConvFvMotionSolver
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
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::velocityLaplacianConvFvMotionSolver::~velocityLaplacianConvFvMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::velocityLaplacianConvFvMotionSolver::curPoints() const
{
    interpolationPtr_->interpolate
    (
        cellMotionU_,
        pointMotionU_
    );

    tmp<pointField> tcurPoints
    (
        fvMesh_.points()
      + fvMesh_.time().deltaTValue()*pointMotionU_.primitiveField()
    );

    twoDCorrectPoints(tcurPoints.ref());

    return tcurPoints;
}


void Foam::velocityLaplacianConvFvMotionSolver::solve()
{
    // The points have moved so before interpolation update
    // the fvMotionSolver accordingly
    movePoints(fvMesh_.points());

    diffusivityPtr_->correct();
    pointMotionU_.boundaryFieldRef().updateCoeffs();
    
    scalar tolerance = 0.001;
    int iter = 0;
    while ( true )
    {
      iter++;

      fvVectorMatrix UEqn
      (
          fvm::laplacian
          (
              diffusivityPtr_->operator()(),
              cellMotionU_,
              "laplacian(diffusivity,cellMotionU)"
          )
      );
      
      //UEqn.relax();
      
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
      else{
        if((iter-1)%100==0)
        {
          Info << " Step " << iter
               << " residual: "<< residual << " > " << tolerance << endl;
          //Info << " iniR " << sp.initialResidual() << nl;
        }
        
      }
      if( iter>1000 )
      {
        Info << nl << "velocityLaplacianConvFvMotionSolver WARNING:"
             << "Laplacian solver did not converge." << nl
             << "Maximum number of iterations"
             << "  iter: "<< iter << endl;
        break;
      }
    }
}


void Foam::velocityLaplacianConvFvMotionSolver::updateMesh
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
