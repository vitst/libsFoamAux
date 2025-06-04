/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::functionObjects::liParticle

Description
    Computes several post processing values for Li3PO4 particle growth case.

    Operands:
    \table
      Operand        | Type            | Location
      input          | -               | -
      output file    | -               | -
      output field   | volScalarField  | $FOAM_CASE/\<time\>/liParticle
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    liParticle1
    {
        // Mandatory entries (unmodifiable)
        type            liParticle;
        libs            (fieldFunctionObjects);

        // Mandatory entries (runtime modifiable)
        geometry
        {
            motorBike.obj
            {
                type triSurfaceMesh;
                name motorBike;
            }
        }

        // Optional entries (runtime modifiable)
        calculateCells  true;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                           | Type | Req'd | Dflt
      type         | Type name: liParticle             | word |  yes  | -
      libs         | Library name: fieldFunctionObjects    | word |  yes  | -
      patchName    | Name of the patch for surface         | word |  yes  | -
    \endtable

SourceFiles
    liParticle.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_liParticle_H
#define functionObjects_liParticle_H

//#include "fieldExpression.H"

#include "fvMeshFunctionObject.H"
#include "volFields.H"
#include "volFieldsFwd.H"
//#include "point.H"
//#include "scalar.H"
//#include "scalarField.H"
//#include "interpolation.H"
#include "volFieldSelection.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
//class objectRegistry;
//class dictionary;
//class mapPolyMesh;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                  Class liParticle Declaration
\*---------------------------------------------------------------------------*/

class liParticle
:
    public fvMeshFunctionObject
    //public fieldExpression
{
protected:

    // Protected Data
        wordList fieldNames_;

        //- Patch name
        word patchName_;

        // kinetics coefs
        scalar nCA_, nCB_;

        // hardcoded parameters TODO move to dict
        const scalar u0 = 0.074; // cm/s
        const scalar h0 = 1e-4; // cm
        const scalar nu_m = 48.24583333333334; // cm^3/mol
        const scalar c0 = 0.0001; // mol/cm^3

        //- result values to write
        scalar curTimeDim; // in s
        scalar minFieldX; // in um/s
        scalar maxFieldX; // in um/s
        scalar minFieldY; // in um/s
        scalar maxFieldY; // in um/s
        scalar minX; // in um
        scalar maxX; // in um
        scalar minY; // in um
        scalar maxY; // in um

        //- surface velocity field
        //scalarField surfVel_;
        //- Face centers
        pointField faceCenters_;

        // to read C fields
        volFieldSelection fieldSet_;
public:

    //- Runtime type information
    TypeName("liParticle");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        liParticle
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        liParticle(const liParticle&) = delete;

        //- No copy assignment
        void operator=(const liParticle&) = delete;


    //- Destructor
    virtual ~liParticle() = default;


    // Member Functions

        //- Read the controls
        virtual bool read(const dictionary&);

        //- Calculate the surfVel
        virtual bool execute();

        //- Write the surfVel and position
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
