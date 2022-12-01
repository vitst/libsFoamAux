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
    Foam::functionObjects::patchSurfaceArea

Group
    grpFieldFunctionObjects

Description
    Computes the distance to the nearest surface from a given geometry.

    Operands:
    \table
      Operand        | Type            | Location
      input          | -               | -
      output file    | -               | -
      output field   | volScalarField  | $FOAM_CASE/\<time\>/patchSurfaceArea
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    patchSurfaceArea1
    {
        // Mandatory entries (unmodifiable)
        type            patchSurfaceArea;
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
      type         | Type name: patchSurfaceArea             | word |  yes  | -
      libs         | Library name: fieldFunctionObjects    | word |  yes  | -
      patchName    | Name of the patch for surface         | word |  yes  | -
    \endtable

SourceFiles
    patchSurfaceArea.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_patchSurfaceArea_H
#define functionObjects_patchSurfaceArea_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"

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
                  Class patchSurfaceArea Declaration
\*---------------------------------------------------------------------------*/

class patchSurfaceArea
:
    public fvMeshFunctionObject
{
protected:

    // Protected Data

        //- Patch name
        word patchName_;

        //- SurfaceArea 
        scalar surfArea_;

        //- TotalVolume
        scalar totalVolume_;

        //- Result file path
        word current_postPr_dir;
        fileName current_file_path;


public:

    //- Runtime type information
    TypeName("patchSurfaceArea");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        patchSurfaceArea
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        patchSurfaceArea(const patchSurfaceArea&) = delete;

        //- No copy assignment
        void operator=(const patchSurfaceArea&) = delete;


    //- Destructor
    virtual ~patchSurfaceArea() = default;


    // Member Functions

        //- Read the controls
        virtual bool read(const dictionary&);

        //- Calculate the curvature
        virtual bool execute();

        //- Write the curvature and position
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
