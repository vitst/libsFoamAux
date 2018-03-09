/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
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

Class
    Foam::velocityDeltatLaplacianFvMotionSolver

Group
    grpMeshMotionSolvers 

Description
    Mesh motion solver for an fvMesh.  Based on solving the cell-centre
    Laplacian for the motion velocity.

SourceFiles
    velocityDeltatLaplacianFvMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef velocityDeltatLaplacianFvMotionSolver_H
#define velocityDeltatLaplacianFvMotionSolver_H

#include "velocityMotionSolver.H"
#include "fvMotionSolver.H"
#include "motionInterpolation.H"
#include "surfaceFields.H"
#include "surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class motionInterpolation;
class motionDiffusivity;

/*---------------------------------------------------------------------------*\
              Class velocityDeltatLaplacianFvMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class velocityDeltatLaplacianFvMotionSolver
:
    public velocityMotionSolver,
    public fvMotionSolver
{
    // Private data
  
        // tolerance for Laplace eq. solver
        scalar tolerance;

        //- Cell-centre motion field
        mutable volVectorField cellMotionU_;

        //- Interpolation used to transfer cell displacement to the points
        autoPtr<motionInterpolation> interpolationPtr_;

        //- Diffusivity used to control the motion
        autoPtr<motionDiffusivity> diffusivityPtr_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        velocityDeltatLaplacianFvMotionSolver
        (
            const velocityDeltatLaplacianFvMotionSolver&
        );

        //- Disallow default bitwise assignment
        void operator=(const velocityDeltatLaplacianFvMotionSolver&);


public:

    //- Runtime type information
    TypeName("velocityDeltatLaplacian");


    // Constructors

        //- Construct from polyMesh and IOdictionary
        velocityDeltatLaplacianFvMotionSolver
        (
            const polyMesh&,
            const IOdictionary&
        );


    //- Destructor
    ~velocityDeltatLaplacianFvMotionSolver();


    // Member Functions

        //- Return reference to the cell displacement field
        volVectorField& cellMotionU()
        {
            return cellMotionU_;
        }

        //- Return const reference to the cell displacement velocity field
        const volVectorField& cellMotionU() const
        {
            return cellMotionU_;
        }
        
        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Update topology
        virtual void updateMesh(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
