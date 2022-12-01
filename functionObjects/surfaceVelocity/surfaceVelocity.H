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
    Foam::functionObjects::surfaceVelocity

Group
    grpFieldFunctionObjects

Description
    Computes the distance to the nearest surface from a given geometry.

    Operands:
    \table
      Operand        | Type            | Location
      input          | -               | -
      output file    | -               | -
      output field   | volScalarField  | $FOAM_CASE/\<time\>/surfaceVelocity
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    surfaceVelocity1
    {
        // Mandatory entries (unmodifiable)
        type            surfaceVelocity;
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
      type         | Type name: surfaceVelocity             | word |  yes  | -
      libs         | Library name: fieldFunctionObjects    | word |  yes  | -
      patchName    | Name of the patch for surface         | word |  yes  | -
    \endtable

SourceFiles
    surfaceVelocity.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_surfaceVelocity_H
#define functionObjects_surfaceVelocity_H

//#include "fieldExpression.H"

#include "fvMeshFunctionObject.H"
#include "volFields.H"
#include "volFieldsFwd.H"
//#include "point.H"
//#include "scalar.H"
//#include "scalarField.H"
//#include "interpolation.H"
//#include "volFieldSelection.H"


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
                  Class surfaceVelocity Declaration
\*---------------------------------------------------------------------------*/

class surfaceVelocity
:
    public fvMeshFunctionObject
    //public fieldExpression
{
protected:

    // Protected Data

        //- Patch name
        word patchName_;

        //- Curvature field
        scalarField surfVel_;
        //- Face centers
        pointField faceCenters_;

        // to read C field
        //volFieldSelection fieldSet_;
public:

    //- Runtime type information
    TypeName("surfaceVelocity");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        surfaceVelocity
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        surfaceVelocity(const surfaceVelocity&) = delete;

        //- No copy assignment
        void operator=(const surfaceVelocity&) = delete;


    //- Destructor
    virtual ~surfaceVelocity() = default;


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
