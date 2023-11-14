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
    Foam::functionObjects::xpatchMotionData

Group
    grpFieldFunctionObjects

Description
    Writes average X pos of a patch and field value time dependence

Usage

SourceFiles
    xpatchMotionData.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_xpatchMotionData_H
#define functionObjects_xpatchMotionData_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"
#include "volFieldsFwd.H"
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
                  Class xpatchMotionData Declaration
\*---------------------------------------------------------------------------*/

class xpatchMotionData
:
    public fvMeshFunctionObject
{
protected:

    // Protected Data
        wordList fieldNames_;

        //- Patch name
        word patchName_;

        //fileName postProc_dir;
        //- ValuesX field
        scalarField values_;
        //- Face centers
        pointField faceCenters_;

        volFieldSelection fieldSet_;

        template<class Type>
        bool foundObject(const word& name, const bool verbose = true) const;
    
        //- get boundary of the field
        template<class Type>
        void getBoundaryValues
        (
            const word& fieldName
        );
public:

    //- Runtime type information
    TypeName("xpatchMotionData");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        xpatchMotionData
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        xpatchMotionData(const xpatchMotionData&) = delete;

        //- No copy assignment
        void operator=(const xpatchMotionData&) = delete;


    //- Destructor
    virtual ~xpatchMotionData() = default;


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
#ifdef NoRepository
    #include "xpatchMotionDataTemplates.C"
#endif

#endif

// ************************************************************************* //
