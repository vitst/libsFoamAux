/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::dissolMotionPointPatchVectorField

Description
  A boundary for the mesh motion according to the dissolution model.
  Currently a neighbor patches should be either: cyclic, symmetry
  or fixedNormalSlip. Processor patch is OK.

See also
    Foam::dissolMotionFvField

SourceFiles
    dissolMotionPointPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef dissolMotionPointPatchVectorField_H
#define dissolMotionPointPatchVectorField_H

#include "fixedValuePointPatchField.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  
typedef List<vectorField> vectorFieldList;
typedef List<scalarField> scalarFieldList;


/*---------------------------------------------------------------------------*\
        Class dissolMotionPointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class dissolMotionPointPatchVectorField
:
    public fixedValuePointPatchField<vector>
{
    // Private data

  vectorFieldList surfWeights;
  //scalarFieldList edgeWeights;
  
  vectorField pinnedPointsNorm;
  labelList pinnedPoints;
  labelList fixedPoints;

  vectorField pinnedPointsNormE;
  labelList pinnedPointsE;
  
public:

    //- Runtime type information
    TypeName("dissolMotion");

    // Constructors

        //- Construct from patch and internal field
        dissolMotionPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        dissolMotionPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patch field onto a new patch
        dissolMotionPointPatchVectorField
        (
            const dissolMotionPointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct as copy
        dissolMotionPointPatchVectorField
        (
            const dissolMotionPointPatchVectorField&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<vector>> clone() const
        {
            return autoPtr<pointPatchField<vector>>
            (
                new dissolMotionPointPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        dissolMotionPointPatchVectorField
        (
            const dissolMotionPointPatchVectorField&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<vector>> clone
        (
            const DimensionedField<vector, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<vector>>
            (
                new dissolMotionPointPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

        void relaxEdges
        (
          vectorField&
        );
        
        void relaxPatchMesh
        (
          vectorField&
        );
        vectorField faceNormals
        (
          const pointField& points, 
          const List<face>& flist
        ) const;
        pointField faceCentres
        (
          const pointField& points,
          const List<face>& flist
        ) const;
        vectorFieldList calc_weights_surface();
        
        void neighborListEdge
        (
          const labelList&,
          const edgeList&,
          const labelListList&,
          labelListList&
        );
        
        
        void getPointMotion
        (
          vectorField&
        );

        void fixCommonNeighborPatchPoints
        (
          vectorField&
        );
        
        void commonPoints
        (
          const labelList&,
          const labelList&, 
          labelList&,
          labelList&,
          labelList&
        );
        

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//    #include "dissolMotionPointPatchVectorField.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
