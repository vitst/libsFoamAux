/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2016 OpenFOAM Foundation
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

Class
    Foam::functionObjects::field2tiff

Group
    grpFieldFunctionObjects

Description
    Calculates the field2tiff, the curl of the velocity.

    The field is stored on the mesh database so that it can be retrieved
    and used for other applications.

Usage
    \verbatim
    field2tiff1
    {
        type        field2tiff;
        libs        ("libfieldFunctionObjects.so");
        ...
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: field2tiff    | yes         |
        U            | Name of velocity field  | no          | U
    \endtable

See also
    Foam::functionObjects::fieldExpression
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    field2tiff.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_field2tiff_H
#define functionObjects_field2tiff_H

#include "fvMeshFunctionObject.H"
#include "volFieldsFwd.H"
#include "point.H"
#include "scalar.H"
#include "scalarField.H"
#include "interpolation.H"
#include "volFieldSelection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class field2tiff Declaration
\*---------------------------------------------------------------------------*/

class field2tiff
:
    public fvMeshFunctionObject
{
protected:
    // Protected member data
        wordList fieldNames_;
        
        point minPoint_, maxPoint_;

        scalar imageUnit_;
        
        //- Fields to integrate
        volFieldSelection fieldSet_;

    // Protected Member Functions

        void printInfo() const;
        
        template<class Type>
        bool foundObject(const word& name, const bool verbose = true) const;
    
        //- Calculate the field integral
        template<class Type>
        void convertFields
        (
            const word& fieldName
        );
        
private:
        //- Disallow default bitwise copy construct
        field2tiff(const field2tiff&);
        
        //- Disallow default bitwise assignment
        void operator=(const field2tiff&);

public:

    //- Runtime type information
    TypeName("field2tiff");

    // Constructors

        //- Construct from Time and dictionary
        field2tiff
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

    //- Destructor
    virtual ~field2tiff();
    
    // Member Functions

        //- Read the fieldExpression data
        virtual bool read(const dictionary&);

        //- Calculate the result field
        virtual bool execute();

        //- Write the result field
        virtual bool write();

        //- Clear the result field from the objectRegistry
        //virtual bool clear();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "field2tiffTemplates.C"
#endif

#endif

// ************************************************************************* //
