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
    Foam::CoupledPatchInterpolation

Description
    Interpolation class within a coupled patch. Allows interpolation from
    points to faces and vice versa
    Modification of PrimitivePatchInterpolation class to synchronize across
    a processor or other couple patch
    ***Only faceToPointInterpolate() has been implemented at present***
    Vitaliy Starchenko (Feb 2015)

SourceFiles
    CoupledPatchInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef CoupledPatchInterpolation_H
#define CoupledPatchInterpolation_H

#include "scalarList.H"
#include "Field.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class CoupledPatchInterpolation Declaration
\*---------------------------------------------------------------------------*/

template<class Patch>
class CoupledPatchInterpolation
{
    // Private data

        //- Reference to patch
        const Patch& patch_;
        
        //- Reference to mesh
        const fvMesh& mesh_;


    // Demand-driven data

        //- Face-to-point weights
        mutable scalarListList* faceToPointWeightsPtr_;
        
        //- Face-to-point sum weights
        mutable scalarList* faceToPointSumWeightsPtr_;

        //- Face-to-edge weights
        mutable scalarList* faceToEdgeWeightsPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        CoupledPatchInterpolation(const CoupledPatchInterpolation&);

        //- Disallow default bitwise assignment
        void operator=(const CoupledPatchInterpolation&);


        //- Face-to-point weights
        const scalarListList& faceToPointWeights() const;

        //- Face-to-point sum weights
        const scalarList& faceToPointSumWeights() const;

        //- Make face-to-point weights
        void makeFaceToPointWeights() const;

        //- Face-to-edge weights
        const scalarList& faceToEdgeWeights() const;

        //- Make face-to-edge weights
        void makeFaceToEdgeWeights() const;

        //- Clear weights
        void clearWeights();


public:

    // Constructors
      //- Construct from PrimitivePatch, and fvMesh
      CoupledPatchInterpolation(const Patch& p, const fvMesh& mesh);

    //- Destructor
      ~CoupledPatchInterpolation();


    // Member Functions

        //- Interpolate from faces to points
        template<class Type>
        tmp<Field<Type> > faceToPointInterpolate
        (
            const Field<Type>& ff
        ) const;

        template<class Type>
        tmp<Field<Type> > faceToPointInterpolate
        (
            const tmp<Field<Type> >& tff
        ) const;

        //- Interpolate from points to faces
        template<class Type>
        tmp<Field<Type> > pointToFaceInterpolate
        (
            const Field<Type>& pf
        ) const;

        template<class Type>
        tmp<Field<Type> > pointToFaceInterpolate
        (
            const tmp<Field<Type> >& tpf
        ) const;

        //- Interpolate from faces to edges
        template<class Type>
        tmp<Field<Type> > faceToEdgeInterpolate
        (
            const Field<Type>& ff
        ) const;

        template<class Type>
        tmp<Field<Type> > faceToEdgeInterpolate
        (
            const tmp<Field<Type> >& tff
        ) const;

        //- Do what is neccessary if the mesh has moved
        bool movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "CoupledPatchInterpolation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
