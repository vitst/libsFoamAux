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
  scalar scaleMotion;
  
  // read parameters
  scalar rlxTol;
  
  bool surfaceRlx;
  int q_norm_recalc;
  scalar k_1, k_2;
  int q_2;
  
  bool edgeRlx;
  int q_norm_recalc_edge;
  scalar k_1edge, k_2edge;
  int q_2edge;
  
  bool pinnedPoint;


  // auxiliary parameters
  bool listsUpdated, weightsUpdated;
  
  // fields
  vectorFieldList surfWeights;
  
  vectorFieldList surfPointWeights;

  labelList pinnedPoints;
  vectorField pinnedPointsNorm;
  
  //labelList rlxEdgePoints;
  labelList fixedPoints;
  labelList globalFixedPoints;
  vectorField fixedPointNorms;
  scalarFieldList rlxEdgeWeights;
  
  labelListList nepe;
  
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
        void fixed_p_edges(labelList&, vectorField&);
        
        void calc_weights_surface();
        
        void calc_point_weights_surface();
        
        void make_lists_and_normals();

        void relaxEdges
        (
          vectorField&
        );
        
        void relaxPatchMesh
        (
          vectorField&
        );
        void relaxPatchMeshPoints
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
        void fixPinnedPoints
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
        
        // set / get
        bool getListsUpdated()
        {
          return listsUpdated;
        }
        void setListsUpdated(bool v)
        {
          listsUpdated = v;
        }
        
        bool getWeightsUpdated()
        {
          return weightsUpdated;
        }
        void setWeightsUpdated(bool v)
        {
          weightsUpdated = v;
        }

        bool getSurfaceRlx()
        {
          return surfaceRlx;
        }
        void setSurfaceRlx(bool v)
        {
          Info<<"surfaceRlx in dissolMotion boundary is set to: "<< v << nl;
          surfaceRlx = v;
        }
  
        bool getEdgeRlx()
        {
          return edgeRlx;
        }
        void setEdgeRlx(bool v)
        {
          Info<<"edgeRlx in dissolMotion boundary is set to: "<< v << nl;
          edgeRlx = v;
        }

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
