/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
    Foam::functionObjects::curv

Group
    grpFieldFunctionObjects

Description
    Computes the second invariant of
    the velocity gradient tensor \f$[s^{-2}]\f$.

    \f[
        \vec curv = \frac{1}{2}
        [ (tr (\nabla \vec U) )^2
        - tr (\nabla \vec U \cdot \nabla \vec U) ]
    \f]

    where
    \vartable
        \vec U       | velocity [m/s]
    \endvartable

    Operands:
    \table
      Operand        | Type           | Location
      input          | volVectorField | $FOAM_CASE/\<time\>/\<inpField\>
      output file    | -              | -
      output field   | volScalarField | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    curv1
    {
        // Mandatory entries (unmodifiable)
        type        curv;
        libs        (fieldFunctionObjects);

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: curv                       | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fieldExpression.H \endlink

    Minimal example by using the \c postProcess utility:
    \verbatim
        postProcess -func curv
    \endverbatim

See also
    - Foam::functionObject
    - Foam::functionObjects::fieldExpression
    - Foam::functionObjects::fvMeshFunctionObject
    - ExtendedCodeGuide::functionObjects::field::curv

SourceFiles
    curv.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_curv_H
#define functionObjects_curv_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class curv Declaration
\*---------------------------------------------------------------------------*/

class curv
:
    public fieldExpression
{
    // Private Member Functions

        //- Calculate the curv field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("curv");


    // Constructors

        //- Construct from Time and dictionary
        curv
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        curv(const curv&) = delete;

        //- No copy assignment
        void operator=(const curv&) = delete;


    //- Destructor
    virtual ~curv() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
