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
    Foam::functionObjects::apertureAver

Group
    grpFieldFunctionObjects

Description
    Calculates the apertureAver, the curl of the velocity.

    The field is stored on the mesh database so that it can be retrieved
    and used for other applications.

Usage
    \verbatim
    apertureAver1
    {
        type        apertureAver;
        libs        ("libfieldFunctionObjects.so");
        ...
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: apertureAver    | yes         |
        U            | Name of velocity field  | no          | U
        result       | Name of Courant number field | no    | \<function name\>
        log          | Log to standard output  | no          | yes
    \endtable

See also
    Foam::functionObjects::fieldExpression
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    apertureAver.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_apertureAver_H
#define functionObjects_apertureAver_H

#include "fvMeshFunctionObject.H"
#include "volFieldsFwd.H"
#include "point.H"
#include "scalar.H"
#include "scalarField.H"
#include "interpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class apertureAver Declaration
\*---------------------------------------------------------------------------*/

class apertureAver
:
    public fvMeshFunctionObject
{
protected:
    // Protected member data
        wordList patchNames_;
        wordList fieldNames_;
        
        point minPoint_, maxPoint_;
        //int flowDirection_;
        //int integrationDirection_;
        
        int expectedNumberOfIntersections_;
        List<int> integrationPoints_;
        
        //pointList intersections;
        
        label N_, M_, K_;   // N - major, M - lateral
        label N1_, M1_, K1_;
        
        int N1M1;
        
        // directions
        label majDir, latDir, intDir;
        
        scalar maxPosMaj, minPosMaj;
        scalar maxPosLat, minPosLat;
        // in case of one surface minPosIntegr is a surface confining 
        // a reservoir instead of second wall
        scalar maxPosInt, minPosInt;
        
        scalar dx, dy; // dx - major, dy - lateral
        
        List<point> pointsXYonsurface;
        
        ///////////////////////
        const static int maxNumProcPoints = 1000000; // max size of the pattern
        int thisTimeSize;
        int totNumLoop;
        int sizeAA;
        int curNum;
        int curBlock;
        int curEnd;        
        ///////////////////////
    
        autoPtr<interpolation<vector> > interpolatorU;
        autoPtr<interpolation<scalar> > interpolatorC;
        
    // Protected Member Functions

        void printInfo() const;

        template<class Type>
        bool foundObject(const word& name, const bool verbose = true) const;
    
        scalar primitive_simpson_integration
        (
            scalarField& x,
            scalarField& y  
        );
        
        void build_surface_points
        (
          const fvMesh& mesh
        );
        
        void integrate_and_write
        (
            const fvMesh& mesh,
            const Time& runTime
        );
        
private:
        //- Disallow default bitwise copy construct
        apertureAver(const apertureAver&);
        
        //- Disallow default bitwise assignment
        void operator=(const apertureAver&);

public:

    //- Runtime type information
    TypeName("apertureAver");

    // Constructors

        //- Construct from Time and dictionary
        apertureAver
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

    //- Destructor
    virtual ~apertureAver();
    
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
    #include "apertureAverTemplates.C"
#endif

#endif

// ************************************************************************* //
