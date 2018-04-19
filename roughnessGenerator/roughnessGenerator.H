/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::RoughnessGenerator

Description
    General class to generate roughness on 2D surface

SourceFiles
    RoughnessGenerator.C

\*---------------------------------------------------------------------------*/

#ifndef roughnessGenerator_H
#define roughnessGenerator_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class RoughnessGenerator Declaration
\*---------------------------------------------------------------------------*/

class RoughnessGenerator
{
    // Protected data
    
        int seed;
        int majNum;
        int minNum;
        double majLen;
        double minLen;
        double rgh;
        double dHurst;        // Fractal dimension D = 3 - dHurst
        double cutLen;
        double maxDisp;
    
    // Private functions

        //- Return the linear index from the i-j indices
        label index
        (
            const label i,
            const label j
        );

        scalar power(double ksq);

public:

    // Constructors

        //- Construct from fvMesh
        RoughnessGenerator
        (
            int seed_,
            int majNum_,
            int minNum_,
            double majLen_,
            double minLen_,
            double rgh_,
            double dHurst_,
            double cutLen_,
            double maxDisp_
        );

    // Member Functions

        // Access

        void getFFTdisp(scalarField& disp);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
