/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::steadyStateControl

Description
    Modification of simpleControl class to supply convergence
    information/checks for a steady state solution, without
    incrementing the time counter.
    Vitaliy Starchenko (Dec 2017)


\*---------------------------------------------------------------------------*/

#ifndef steadyStateControl_H
#define steadyStateControl_H

#include "solutionControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class steadyStateControl Declaration
\*---------------------------------------------------------------------------*/

class steadyStateControl
:
    public solutionControl
{

protected:

    // Protected Data

        //- Initialised flag
        bool initialised_;


    // Protected Member Functions

        //- Read controls from fvSolution dictionary
        void read();

        //- Return true if all convergence checks are satisfied
        bool criteriaSatisfied();

        //- Disallow default bitwise copy construct
        steadyStateControl(const steadyStateControl&);

        //- Disallow default bitwise assignment
        void operator=(const steadyStateControl&);

    	int iter_counter;

        template<class Type>
        void maxTypeResidual
        (
            const word& fieldName,
            ITstream& data,
            scalar& firstRes,
            scalar& lastRes
        ) const;

        scalar maxResidual
        (
            const word& fieldName,
            ITstream& data,
            scalar& lastRes
        ) const;

public:


    // Static Data Members

        //- Run-time type information
        TypeName("steadyStateControl");


    // Constructors

        //- Construct from mesh
        steadyStateControl(fvMesh& mesh);


    //- Destructor
    virtual ~steadyStateControl() = default;


    // Member Functions

        // Solution control

            //- Loop loop
            virtual bool loop();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
