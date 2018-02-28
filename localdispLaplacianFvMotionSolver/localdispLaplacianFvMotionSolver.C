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

#include "localdispLaplacianFvMotionSolver.H"
#include "motionInterpolation.H"
#include "motionDiffusivity.H"
#include "fvmLaplacian.H"
#include "fvcLaplacian.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(localdispLaplacianFvMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        localdispLaplacianFvMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::localdispLaplacianFvMotionSolver::localdispLaplacianFvMotionSolver
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
    W_
    (
        IOobject
        (
            "inverseDiffusivity0",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvMesh_,
        dimensionedScalar("1.0", dimless, 1.0)
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
    W_ = diffusivityPtr_->operator()();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::localdispLaplacianFvMotionSolver::~localdispLaplacianFvMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::localdispLaplacianFvMotionSolver::curPoints() const
{
    interpolationPtr_->interpolate
    (
        cellMotionU_,
        pointMotionU_
    );
    
    tmp<pointField> tcurPoints
    (
        fvMesh_.points()
      + pointMotionU_.primitiveField()
      //+ fvMesh_.time().deltaTValue()*pointMotionU_.primitiveField()
    );

    twoDCorrectPoints(tcurPoints.ref());

    return tcurPoints;
}


void Foam::localdispLaplacianFvMotionSolver::solve()
{
    // The points have moved so before interpolation update
    // the fvMotionSolver accordingly
    movePoints(fvMesh_.points());

    diffusivityPtr_->correct();
    pointMotionU_.boundaryFieldRef().updateCoeffs();
    
    //scalarField W = fvMesh_.V();
    //W = 1/W;
    //Info<<"diffusivity:  "<<diffusivityPtr_->operator()()<<nl;
    //Info<<"diffusivity:  "<<W_<<nl;
    
    
    /*
    volVectorField Y
    (
        IOobject
        (
            "Y",
            fvMesh_.time().timeName(),
            fvMesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        fvMesh_,
        dimensionedVector("Y", dimensionSet(0,1,0,0,0,0,0), vector::zero),
        fixedValueFvPatchVectorField::typeName
        //zeroGradientFvPatchVectorField::typeName
    );
    Y.primitiveFieldRef() = fvMesh_.C(); // fvMesh_.time().deltaTValue();
    //Y.correctBoundaryConditions();

    forAll(Y.boundaryField(), I)
    {
       vectorField& bvField = Y.boundaryFieldRef()[I];
       bvField = fvMesh_.boundaryMesh()[I].faceCentres(); // / fvMesh_.time().deltaTValue();
    } 
    */
    
    
    /*
    Info<<"fvc::laplacian    "
            <<
          fvc::laplacian
          (
              W_,
              Y,
              "laplacian(diffusivity,Y)"
          )
          << nl;
    //std::exit(0);
     */
    
    //Info<<"max pm BEF: "<< max(pointMotionU_.primitiveField())<<nl;
    //pointMotionU_ = pointMotionU_ * fvMesh_.time().deltaTValue();
    //Info<<"max pm AFT: "<< max(pointMotionU_.primitiveField())<<nl;
    
    //Info<<"W size   "<<W_.size()<<nl;
    //Info<<"cellMotionU_   "<<cellMotionU_.size()<<nl;
    
    scalar tolerance = 0.0001;
    int iter = 0;
    while ( true )
    {
      iter++;

      fvVectorMatrix UEqn
      (
          fvm::laplacian
          (
              diffusivityPtr_->operator()(),
//              W_,
              cellMotionU_,
              "laplacian(diffusivity,cellMotionD)"
          )
      /*
      +
          fvc::laplacian
          (
              W_,
              Y,
              "laplacian(diffusivity,Y)"
          )
       */
      );
      
/*
*/
      
      //UEqn.relax();
      
      //Y = Y + cellMotionU_*fvMesh_.time().deltaTValue();
      
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
          Info << " Step " << iter << token::TAB
               << " residual: "<< residual << " > " << tolerance << endl;
          //Info << " iniR " << sp.initialResidual() << nl;
        }
        
      }
      if( iter>1000 )
      {
        Info << nl << "localdispLaplacianFvMotionSolver WARNING:"
             << "Laplacian solver did not converge." << nl
             << "Maximum number of iterations"
             << "  iter: "<< iter << endl;
        break;
      }
    }
    
    //Info<<"max cm BEF: "<< max(cellMotionU_.primitiveField())<<nl;
    //pointMotionU_ = pointMotionU_ / fvMesh_.time().deltaTValue();
    //cellMotionU_ = cellMotionU_ / fvMesh_.time().deltaTValue();
    //Info<<"max cm AFT: "<< max(cellMotionU_.primitiveField())<<nl;
}


void Foam::localdispLaplacianFvMotionSolver::updateMesh
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
