/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "dissolMotionPointPatchVectorField.H"
#include "Time.H"
#include "IFstream.H"

#include "coupledPatchInterpolation.H"
#include "volFields.H"
#include "pointPatchFields.H"

#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "polyMesh.H"

//#include "wedgePolyPatch.H"
#include "cyclicPolyPatch.H"
#include "symmetryPolyPatch.H"
#include "processorPolyPatch.H"
//#include "emptyPolyPatch.H"
//#include "fixedNormalSlipPointPatchField.H"

#include "plane.H"
#include "face.H"
#include "syncTools.H"

//#include "triSurface.H"
//#include "triSurfaceTools.H"
//#include "triSurfaceSearch.H"
//#include "meshSearch.H"
//#include "pyramidPointFaceRef.H"
//#include "triPointRef.H"

#include "timeSelector.H"
#include "argList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
dissolMotionPointPatchVectorField::
dissolMotionPointPatchVectorField
(
  const pointPatch& p,
  const DimensionedField<vector, pointMesh>& iF
)
:
  fixedValuePointPatchField<vector>(p, iF),
  scaleMotion(1.0),
  rlxTol(1e-6),
  surfaceRlx(true),
  q_norm_recalc(1),
  k_1(1.0),
  k_2(1.0),
  q_2(1),
  edgeRlx(true),
  q_norm_recalc_edge(1),
  k_1edge(1.0),
  k_2edge(1.0),
  q_2edge(1),
  pinnedPoint(false)
{
  if(debug) 
  {
    Info << "dissolMotionPointPatchVectorField constructor 0"<<endl;
  }
  this->setListsUpdated(false);
  this->setWeightsUpdated(false);
}


Foam::
dissolMotionPointPatchVectorField::
dissolMotionPointPatchVectorField
(
    const dissolMotionPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
  fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
  scaleMotion(ptf.scaleMotion),
  rlxTol(ptf.rlxTol),
  surfaceRlx(ptf.surfaceRlx),
  q_norm_recalc(ptf.q_norm_recalc),
  k_1(ptf.k_1),
  k_2(ptf.k_2),
  q_2(ptf.q_2),
  edgeRlx(ptf.edgeRlx),
  q_norm_recalc_edge(ptf.q_norm_recalc_edge),
  k_1edge(ptf.k_1edge),
  k_2edge(ptf.k_2edge),
  q_2edge(ptf.q_2edge),
  pinnedPoint(ptf.pinnedPoint)
{
  if(debug) 
  {
    Info << "dissolMotionPointPatchVectorField constructor 1"<<endl;
  }
  
  this->setListsUpdated(false);
  this->setWeightsUpdated(false);
  //this->operator==(timeSeries_(this->db().time().timeOutputValue()));
}


Foam::
dissolMotionPointPatchVectorField::
dissolMotionPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
  fixedValuePointPatchField<vector>(p, iF)
{
  if(debug) {
    Info << "dissolMotionPointPatchVectorField constructor 2"<<endl;
  }
  
  if (dict.found("value"))
  {
    this->operator==
    (
      vectorField("value", dict, p.size())
    );
  }
  else
  {
    WarningIn(
        "dissolMotionPointPatchVectorField constructor 2"
        "("
        "const pointPatch& p,"
        "const DimensionedField<vector, pointMesh>& iF,"
        "const dictionary& dict"
        ")"
    ) << "No value defined for " << this->internalField().name()
        << " on " << this->patch().name()
        << endl;
  }
  
  
  if (!dict.readIfPresent<scalar>("scaleMotion", scaleMotion))
  {
    scaleMotion = 1.0;
    WarningIn
    (
      "dissolMotionPointPatchVectorField constructor 2"
      "("
      "const fvPatch& p,"
      "const DimensionedField<Type, volMesh>& iF,"
      "const dictionary& dict"
      ")"
    )
      << "No value defined for scaleMotion"
      << " on " << this->patch().name() << " therefore using "
      << scaleMotion
      << endl;
  }

  if (!dict.readIfPresent<scalar>("rlxTol", rlxTol))
  {
    rlxTol = 1.0e-6;
    WarningIn
    (
      "dissolMotionPointPatchVectorField constructor 2"
      "("
      "const fvPatch& p,"
      "const DimensionedField<Type, volMesh>& iF,"
      "const dictionary& dict"
      ")"
    )
      << "No value defined for rlxTol"
      << " on " << this->patch().name() << " therefore using "
      << rlxTol
      << endl;
  }

  if (!dict.readIfPresent<bool>("surfaceRlx", surfaceRlx))
  {
    surfaceRlx = true;
    WarningIn
    (
      "dissolMotionPointPatchVectorField constructor 2"
      "("
      "const fvPatch& p,"
      "const DimensionedField<scalar, volMesh>& iF,"
      "const dictionary& dict"
      ")"
    )
      << "No value defined for surfaceRlx"
      << " on " << this->patch().name() << " therefore using true"
      << endl;
  }

  if (!dict.readIfPresent<int>("q_norm_recalc", q_norm_recalc))
  {
    q_norm_recalc = 1;
    WarningIn
    (
      "dissolMotionPointPatchVectorField constructor 2"
      "("
      "const fvPatch& p,"
      "const DimensionedField<Type, volMesh>& iF,"
      "const dictionary& dict"
      ")"
    )
      << "No value defined for q_norm_recalc"
      << " on " << this->patch().name() << " therefore using "
      << q_norm_recalc
      << endl;
  }
  if (!dict.readIfPresent<scalar>("k_1", k_1))
  {
    k_1 = 1.0;
    WarningIn
    (
      "dissolMotionPointPatchVectorField constructor 2"
      "("
      "const fvPatch& p,"
      "const DimensionedField<Type, volMesh>& iF,"
      "const dictionary& dict"
      ")"
    )
      << "No value defined for k_1"
      << " on " << this->patch().name() << " therefore using "
      << k_1
      << endl;
  }
  if (!dict.readIfPresent<scalar>("k_2", k_2))
  {
    k_2 = 1.0;
    WarningIn
    (
      "dissolMotionPointPatchVectorField constructor 2"
      "("
      "const fvPatch& p,"
      "const DimensionedField<Type, volMesh>& iF,"
      "const dictionary& dict"
      ")"
    )
      << "No value defined for k_2"
      << " on " << this->patch().name() << " therefore using "
      << k_2
      << endl;
  }
  if (!dict.readIfPresent<int>("q_2", q_2))
  {
    q_2 = 1;
    WarningIn
    (
      "dissolMotionPointPatchVectorField constructor 2"
      "("
      "const fvPatch& p,"
      "const DimensionedField<Type, volMesh>& iF,"
      "const dictionary& dict"
      ")"
    )
      << "No value defined for q_2"
      << " on " << this->patch().name() << " therefore using "
      << q_2
      << endl;
  }


  if (!dict.readIfPresent<bool>("edgeRlx", edgeRlx))
  {
    edgeRlx = true;
    WarningIn
    (
      "dissolMotionPointPatchVectorField constructor 2"
      "("
      "const fvPatch& p,"
      "const DimensionedField<scalar, volMesh>& iF,"
      "const dictionary& dict"
      ")"
    )
      << "No value defined for edgeRlx"
      << " on " << this->patch().name() << " therefore using true"
      << endl;
  }
  if (!dict.readIfPresent<int>("q_norm_recalc_edge", q_norm_recalc_edge))
  {
    q_norm_recalc_edge = 1;
    WarningIn
    (
      "dissolMotionPointPatchVectorField constructor 2"
      "("
      "const fvPatch& p,"
      "const DimensionedField<Type, volMesh>& iF,"
      "const dictionary& dict"
      ")"
    )
      << "No value defined for q_norm_recalc_edge"
      << " on " << this->patch().name() << " therefore using "
      << q_norm_recalc_edge
      << endl;
  }
  if (!dict.readIfPresent<scalar>("k_1edge", k_1edge))
  {
    k_1edge = 1.0;
    WarningIn
    (
      "dissolMotionPointPatchVectorField constructor 2"
      "("
      "const fvPatch& p,"
      "const DimensionedField<Type, volMesh>& iF,"
      "const dictionary& dict"
      ")"
    )
      << "No value defined for k_1edge"
      << " on " << this->patch().name() << " therefore using "
      << k_1edge
      << endl;
  }
  if (!dict.readIfPresent<scalar>("k_2edge", k_2edge))
  {
    k_2edge = 1.0;
    WarningIn
    (
      "dissolMotionPointPatchVectorField constructor 2"
      "("
      "const fvPatch& p,"
      "const DimensionedField<Type, volMesh>& iF,"
      "const dictionary& dict"
      ")"
    )
      << "No value defined for k_2edge"
      << " on " << this->patch().name() << " therefore using "
      << k_2edge
      << endl;
  }
  if (!dict.readIfPresent<int>("q_2edge", q_2edge))
  {
    q_2edge = 1;
    WarningIn
    (
      "dissolMotionPointPatchVectorField constructor 2"
      "("
      "const fvPatch& p,"
      "const DimensionedField<Type, volMesh>& iF,"
      "const dictionary& dict"
      ")"
    )
      << "No value defined for q_2edge"
      << " on " << this->patch().name() << " therefore using "
      << q_2edge
      << endl;
  }
  
  if (!dict.readIfPresent<bool>("pinnedPoint", pinnedPoint))
  {
    pinnedPoint = false;
    WarningIn
    (
      "dissolMotionPointPatchVectorField constructor 2"
      "("
      "const fvPatch& p,"
      "const DimensionedField<scalar, volMesh>& iF,"
      "const dictionary& dict"
      ")"
    ) 
      << "No value defined for pinnedPoint"
      << " on " << this->patch().name() << " therefore using false"
      << endl;
  }

  this->setListsUpdated(false);
  this->setWeightsUpdated(false);
}

Foam::
dissolMotionPointPatchVectorField::
dissolMotionPointPatchVectorField
(
  const dissolMotionPointPatchVectorField& ptf
)
:
  fixedValuePointPatchField<vector>(ptf),
  scaleMotion(ptf.scaleMotion),
  rlxTol(ptf.rlxTol),
  surfaceRlx(ptf.surfaceRlx),
  q_norm_recalc(ptf.q_norm_recalc),
  k_1(ptf.k_1),
  k_2(ptf.k_2),
  q_2(ptf.q_2),
  edgeRlx(ptf.edgeRlx),
  q_norm_recalc_edge(ptf.q_norm_recalc_edge),
  k_1edge(ptf.k_1edge),
  k_2edge(ptf.k_2edge),
  q_2edge(ptf.q_2edge),
  pinnedPoint(ptf.pinnedPoint)
{
  if(debug)
  {
    Info << "dissolMotionPointPatchVectorField constructor 3"<<endl;
  }
  this->setListsUpdated(false);
  this->setWeightsUpdated(false);
}


Foam::
dissolMotionPointPatchVectorField::
dissolMotionPointPatchVectorField
(
  const dissolMotionPointPatchVectorField& ptf,
  const DimensionedField<vector, pointMesh>& iF
)
:
  fixedValuePointPatchField<vector>(ptf, iF),
  scaleMotion(ptf.scaleMotion),
  rlxTol(ptf.rlxTol),
  surfaceRlx(ptf.surfaceRlx),
  q_norm_recalc(ptf.q_norm_recalc),
  k_1(ptf.k_1),
  k_2(ptf.k_2),
  q_2(ptf.q_2),
  edgeRlx(ptf.edgeRlx),
  q_norm_recalc_edge(ptf.q_norm_recalc_edge),
  k_1edge(ptf.k_1edge),
  k_2edge(ptf.k_2edge),
  q_2edge(ptf.q_2edge),
  pinnedPoint(ptf.pinnedPoint)
{
  if(debug) {
    Info << "dissolMotionPointPatchVectorField constructor 4"<<endl;
  }
  this->setListsUpdated(false);
  this->setWeightsUpdated(false);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dissolMotionPointPatchVectorField::updateCoeffs()
{
  if (this->updated())
  {
    return;
  }
  
  if(debug) 
  {
    Info << "dissolMotionPointPatchVectorField::updateCoeffs()"<<endl;
  }
  
  if( this->db().foundObject<IOdictionary>("transportProperties") )
  {
    if(debug) 
    {
      Info << "dissolMotionPointPatchVectorField::updateCoeffs() "
              "calc_weights_surface"<<nl;
    }

    if( !this->getWeightsUpdated() )
      calc_weights_surface();

    if(debug) 
    {
      Info << "dissolMotionPointPatchVectorField::updateCoeffs() "
              "make_lists_and_normals"<<nl;
    }

    if( !this->getListsUpdated() )
      make_lists_and_normals();
    
    
    //label patchID = this->patch().index();
    const scalar dt = this->db().time().deltaTValue();

    //const polyMesh& mesh = this->internalField().mesh()();

    vectorField pointMotion;
    getPointMotion(pointMotion);

    pointMotion *= (dt * scaleMotion);
    
    /*
    vectorField faceMotion;
    getPointMotion(faceMotion);
    faceMotion *= dt;
     */
    
    /*
    Info<<"PPPoint motion: "<<nl;
    for(int i = 0; i<10; i++)
    {
      Info<< pointMotion[i] << nl;
    }
    if(Pstream::myProcNo()==3)
      Pout<< "QQQQQQQQ: fixedPoints[74]: "<< fixedPoints[74]
              <<" pointMotion[fp74]: "<<pointMotion[fixedPoints[74]]
              <<" pointMotion[2507]: "<<pointMotion[2507]
              <<nl<<nl;
     */
    

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // TODO the rest of relaxation should be implemented
    Info<<"fixCommonNeighborPatchPoints"<<nl;
    fixCommonNeighborPatchPoints(pointMotion);

    /*
    if(Pstream::myProcNo()==3)
      Pout<< "FFFFFFFFF: fixedPoints[74]: "<< fixedPoints[74]
              <<" pointMotion[fp74]: "<<pointMotion[fixedPoints[74]]
              <<" pointMotion[2507]: "<<pointMotion[2507]
              <<nl<<nl;
    */
    
    if(edgeRlx)
    {
      Info<<nl<<"relaxEdges"<<nl;
      relaxEdges(pointMotion);
    }
    if(surfaceRlx)
    {
      Info<<nl<<"relaxPatchMesh"<<nl;
      calc_weights_surface();
      relaxPatchMesh(pointMotion);
    }
    
    //calc_point_weights_surface();
    //relaxPatchMeshPoints(pointMotion);
    
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    /*
    Info<<"Aftern: "<<nl;
    for(int i = 0; i<10; i++)
    {
      Info<< pointMotion[i] << nl;
    }
     */
    

    // set the velocity for the point motion
    this->operator==( pointMotion/dt );
    
    /*
    faceMotion /= dt;
    
    label patchID = this->patch().index();
    const polyMesh& mesh = this->internalField().mesh()();
    //- set-up interpolator
    const fvMesh& fvmesh_ = refCast<const fvMesh>(mesh);
    coupledPatchInterpolation patchInterpolator
    (
      mesh.boundaryMesh()[patchID], fvmesh_
    );
    vectorField pointMotion = 
            patchInterpolator.faceToPointInterpolate(faceMotion);
    
    this->operator==( pointMotion );
     */

  }
  fixedValuePointPatchField<vector>::updateCoeffs();
}

void Foam::
dissolMotionPointPatchVectorField::
relaxEdges(vectorField& pointMotion)
{
  label patchID = this->patch().index();
  const polyMesh& mesh = this->internalField().mesh()();

  const polyPatch& curPatch = mesh.boundaryMesh()[patchID];
  
  // current patch geometry
  const pointField& curPP  = curPatch.localPoints();
  //const labelList& curMeshPoints = curPatch.meshPoints();
  //const labelListList& curPointEdges = curPatch.pointEdges();
  //const edgeList& ee = curPatch.edges();
  
  const List<face>& llf = curPatch.localFaces();
  const labelListList& plistFaces = curPatch.pointFaces();

  //const polyBoundaryMesh& bMesh = mesh.boundaryMesh();
  
  labelList fixedEdgePoint;
  
  fixed_p_edges(fixedEdgePoint, pointMotion);
  
  int NN = fixedPoints.size();
  
  /*
  if(fixedEdgePoint.size()>0)
    Pout  << " fixEdgePoint: " << fixedEdgePoint
          << " fixPointsSize: " << fixedPoints.size()
          << " pntPos: " << curPP[fixedPoints[fixedEdgePoint[0]]]
          << " pntMotion: " << pointMotion[fixedPoints[fixedEdgePoint[0]]]
          <<endl;

  if(Pstream::myProcNo()==3)
    Pout<< "AAAAAAAA: fixedPoints[74]: "<< fixedPoints[74]
            <<" pointMotion[fp74]: "<<pointMotion[fixedPoints[74]]
            <<" curPP[2507]: "<<curPP[2507]
            <<" pointMotion[2507]: "<<pointMotion[2507]
            <<nl<<nl;
  */



  labelList loc_corners;
  //labelList loc_corners_fP;
  labelListList frecBin(curPP.size());
  forAll(fixedPoints, i)
  {
    //labelList& fpl = frecBin[fixedPoints[i]];
    frecBin[fixedPoints[i]].append(i);
  }
  forAll(frecBin, i)
  {
    labelList& fpl = frecBin[i];
    if(fpl.size()>1)
    {
      loc_corners.append(i);
    }
  }
  vectorField loc_corners_norms(loc_corners.size());
  forAll(loc_corners,i)
  {
    labelList& fpl = frecBin[loc_corners[i]];
    
    forAll(fpl, j)
    {
      label crnp = fpl[j];
      vector crnv = fixedPointNorms[crnp];
      if(j==0)
      {
        loc_corners_norms[i] = crnv;
      }
      else
      {
        loc_corners_norms[i] = loc_corners_norms[i] ^ crnv;
      }
    }
  }
  
  /*
  if(Pstream::myProcNo()==3)
    Pout<< "CCCCCCCC: "<<loc_corners.size()
            <<" pointMotion[fp74]: "<<pointMotion[fixedPoints[74]]
            <<" curPP[2507]: "<<curPP[2507]
            <<" pointMotion[2507]: "<<pointMotion[2507]
            <<nl<<nl;
   */
  // fix corners
  forAll(loc_corners, i)
  {
    pointMotion[loc_corners[i]] = 
            pointMotion[loc_corners[i]] 
            & 
            (loc_corners_norms[i] * loc_corners_norms[i]);
    
    /*
    if(Pstream::myProcNo()==3)
      Pout<< "Corners: "<<loc_corners[i]
          <<"  loc_corners_norms "<<loc_corners_norms[i]
            << "  pM : " << pointMotion[loc_corners[i]]
            <<nl<<nl;
    / *
    Info<< "Corners: "<<loc_corners[i]
          <<"  loc_corners_norms "<<loc_corners_norms[i]
            << "  pM: " << pointMotion[loc_corners[i]]
            <<nl;
     */
  }
  //std::exit(0);
  
  vectorField pointNorm( NN, vector::zero );
  scalarList faceToPointSumWeights( NN, 0.0 );

  /*
  if(Pstream::myProcNo()==3)
    Pout<< "DDDDDDDD: "<<loc_corners.size()
            <<" pointMotion[2507]: "<<pointMotion[2507]
            <<" curPP[2507]: "<<curPP[2507]
            <<nl<<nl;
   */
  
  if(debug) 
  {
    Info << "  Number of corner points "<< loc_corners.size() << nl;
  }
  
  pointField movedPoints = curPP + pointMotion;
  
  if(debug) 
  {
    Info << "  Relaxation loop:"<< nl;
  }
  double displ_tol = 1.0;
  int itt = 0;
  while(displ_tol>rlxTol)
  {
    if(itt%q_norm_recalc_edge==0)
    {
      // set fields to zero
      pointNorm = vector::zero;
      faceToPointSumWeights = 0.0;
    }

    pointField faceCs = faceCentres(movedPoints, llf);
    vectorField faceNs = faceNormals(movedPoints, llf);

    vectorField displacement(NN, vector::zero);
    scalarField tol(NN, 0.0);
    scalarList sumWeights( NN, 0.0 );

    //forAll(nepe, i)
    forAll(fixedPoints, i)
    {
      label  curI = fixedPoints[i];
      point& curP = movedPoints[curI];
      const labelList& pNeib = nepe[i];

      const scalarField& curwv = rlxEdgeWeights[i]; //weights[i];

      forAll(pNeib, ii)
      {
        label ind = pNeib[ii];
        point& neibP = movedPoints[ind];
        vector d2 = neibP - curP;

        scalar mag_d = mag(d2);

        displacement[i] += curwv[ii] * d2;
        tol[i] += mag_d;
        sumWeights[i] += 1.0;
      }
      
      /*
      if(itt%1000==0 && i>72 && i<76)
      {
        if(Pstream::myProcNo()==3)
          Pout<<i<<"  beforeSync: "<<displacement[i]
                  <<"  curP " << curP
                  <<"  pointMotion: "<<pointMotion[curI]
                  <<"  neP1 " << movedPoints[pNeib[0]]
                  <<"  neP2 " << movedPoints[pNeib[1]]
                  <<"  pNei " <<  pNeib
                  //<<"  curw " <<  curwv 
                  //<< " tolr "<< tol[i] 
                  << "            itt "<< itt
                  <<nl;
      }
      */
      

      // stick to cyclic boundary

      if(itt%q_norm_recalc_edge==0)
      {

        //label curIpp = local_pp_EdgePoints[i];
        const vector& curNormPP = fixedPointNorms [ i ];

        const labelList& pFaces = plistFaces[curI];
        forAll(pFaces, j)
        {
          label ind = pNeib[j];
          point& neibP = movedPoints[ind];
          vector d2 = neibP - curP;     // vector from the current point to its neighbor

          label faceI = pFaces[j]; 
          vector fnn = faceNs[ faceI ]; // face normal
          fnn = transform(I - curNormPP*curNormPP, fnn);

          // correction of the normal, otherwise realN = fnn
          scalar middi = mag(d2)/2.0;   // half distance to neighbor
          point midpo = curP + d2/2.0;  // midpoint
          plane pll(midpo, d2);           // plane perpendicular to the edge via midpoint
          point endNorm = midpo + fnn;  // end of the face normal projected on the inlet
          point projp = pll.nearestPoint(endNorm); // projection of endNorm onto pll plane

          vector realN = projp - midpo; // normal to the edge between current point
                                        // and its neighbor in the inlet plane facing
                                        // outside the fracture (mag(realN)!=1 !!!)

          scalar nw = 1.0 / middi;
          pointNorm[i] += nw * realN;
          pointNorm[i] = transform(I - curNormPP*curNormPP, pointNorm[i]);
          faceToPointSumWeights[i] += nw;
        }
      }
    }
    
    syncTools::syncPointList(mesh, globalFixedPoints, displacement, plusEqOp<vector>(), vector::zero);
    syncTools::syncPointList(mesh, globalFixedPoints, tol, plusEqOp<scalar>(), 0.0);
    syncTools::syncPointList(mesh, globalFixedPoints, sumWeights, plusEqOp<scalar>(), 0.0);

    /*
    if(itt%1000==0)
    {
      int iii = findMax(displacement);
      if(iii>=0)
        Pout<<" afterSyncMax: "<<max(displacement)<<"  "
                << displacement[0]<< "  "
                << iii << "  "
                << sumWeights[iii] << "  "
                << movedPoints[iii]
                <<nl;
    }
    */
    
    forAll(pointNorm, i)
    {
      displacement[i] /= sumWeights[i];
    }

    if(itt%q_norm_recalc_edge==0)
    {
      syncTools::syncPointList(mesh, globalFixedPoints, pointNorm, plusEqOp<vector>(), vector::zero);
      syncTools::syncPointList(mesh, globalFixedPoints, faceToPointSumWeights, plusEqOp<scalar>(), 0.0);
      // normalization
      forAll(pointNorm, i)
      {
        pointNorm[i] /= mag( pointNorm[i] );
      }
    }


    forAll(fixedEdgePoint, i1)
    {
      //displacement[fixedEdgePoint[i1]] = vector::zero;
    }

    vectorField projectedDisplacement = 
            transform(I - pointNorm*pointNorm, displacement);

    forAll(loc_corners, i)
    {
      labelList& fpl = frecBin[loc_corners[i]];
      forAll(fpl, j)
      {
        label lid = fpl[j];
        projectedDisplacement[lid] = vector::zero;
      }
    }
    
    scalar factor = (itt%q_2edge==0) ? k_2edge : k_1edge;
    projectedDisplacement *= factor;

    /*
    for(int i=0; i<3; i++){
      if(projectedDisplacement.size()>0)
        Info<<"i="<<i<<"  pr: "<<projectedDisplacement[i]
                <<"  pm: "<< pointMotion[i]
                <<endl;
    }
    */
    // stick to cyclic boundary
    //fixPinnedPoints(projectedDisplacement);
    
    forAll(fixedPoints, i)
    {
      label ind = fixedPoints[i];
      movedPoints[ind] += projectedDisplacement[i];
    }
    
    // To avoid WARNING message
    // @TODO need to be checked for correctness
    scalar gAmpD = 0.0;
    scalar gAmt  = 0.0;
    //if(projectedDisplacement.size()>0){
      gAmpD = gAverage( mag(projectedDisplacement) );
      gAmt  = gAverage( mag(tol) );
    //}
    

    if(gAmt>SMALL)
      displ_tol = gAmpD /factor/ gAmt;
    else
      displ_tol = 0.0;
      
    if(itt%100==0)
    {
      //Pout<<" maxDDD: "<<max(projectedDisplacement)<<nl;
      Info << "  edge rlx iter " << itt
           << " tolerance: " << displ_tol << endl;
    }

    itt++;
  }
  Info << nl << "  Edge relaxation done in " << itt
       << " iterations. Tolerance: " << displ_tol << endl;

  pointMotion = (movedPoints-curPP);

  fixPinnedPoints(pointMotion);
  
  /*
  forAll(fixedPoints, i)
  {
    Info<< " pntPos: " << curPP[fixedPoints[i]]
        << "  pntMotion: " << pointMotion[fixedPoints[i]]
        << nl
        << "     norms: "<< fixedPointNorms[i]
        << "  weight: " << rlxEdgeWeights[i]
        <<endl;
  }
  */
  
}

void Foam::
dissolMotionPointPatchVectorField::
fixed_p_edges(labelList& fixedPEdges, vectorField& pointMotion)
{
  label patchID = this->patch().index();
  const polyMesh& mesh = this->internalField().mesh()();

  const polyPatch& curPatch = mesh.boundaryMesh()[patchID];

  // current patch geometry
  const pointField& curPP  = curPatch.localPoints();
  const labelList& curMeshPoints = curPatch.meshPoints();
  const labelListList& curPointEdges = curPatch.pointEdges();
  const edgeList& ee = curPatch.edges();

  const List<face>& llf = curPatch.localFaces();
  const labelListList& plistFaces = curPatch.pointFaces();

  const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

  forAll(bMesh, patchi)
  {
    if ( (patchi != patchID) ) //bMesh[patchi].size() && 
    {
      if
      (
        isA<cyclicPolyPatch>(bMesh[patchi]) ||
        isA<symmetryPolyPatch>(bMesh[patchi]) ||
        isA<processorPolyPatch>(bMesh[patchi])
      )
      {
        // skip
      }
      else
      {
        //Info<<nl<< "!!!!!!!!!!!!!!!!!!!!!!!!!!!! patch type:  "<< bMesh[patchi].type()<<nl<<endl;
        const polyPatch& pp = refCast<const polyPatch>(bMesh[patchi]);
        //const pointField& pfPP = pp.localPoints();
        const labelList& ppMeshPoints = pp.meshPoints();
        const vectorField& ppPointNormals = pp.pointNormals();  // 3
        
        labelList local_EdgePoints, global_EdgePoints;
        labelList local_pp_EdgePoints;

        commonPoints
        (
          curMeshPoints,
          ppMeshPoints,
          local_EdgePoints,
          global_EdgePoints,
          local_pp_EdgePoints
        );
        
        /*
        if(local_EdgePoints.size()>0)
          Pout<<"Patch type:  "<< bMesh[patchi].type()
            <<" "<<bMesh[patchi].name()
            <<" lep "<<local_EdgePoints.size()
            <<" gep "<<global_EdgePoints.size()
            <<endl;
         */

        //if(local_EdgePoints.size()>0)
        //{
          labelListList nepe11;
          neighborListEdge(local_EdgePoints, ee, curPointEdges, nepe11);

          label NN = local_EdgePoints.size();


          vectorField pointNorm( NN, vector::zero );
          scalarList faceToPointSumWeights( NN, 0.0 );

          pointField movedPoints = curPP + pointMotion;

          // *************************************************************************
          // calculate old points normals once
          vectorField currentNorms( NN, vector::zero );
          {
            pointField faceCs = faceCentres(curPP, llf);
            vectorField faceNs = faceNormals(curPP, llf);

            scalarList sumWeights( NN, 0.0 );

            forAll(nepe11, i)
            {
              label  curI = local_EdgePoints[i];
              const point& curP = curPP[curI];
              const labelList& pNeib = nepe11[i];
              label curIpp = local_pp_EdgePoints[i];
              const vector& curNormPP = ppPointNormals[ curIpp ];

              const labelList& pFaces = plistFaces[curI];
              forAll(pFaces, j)
              {
                label ind = pNeib[j];
                const point& neibP = curPP[ind];
                vector d2 = neibP - curP;     // vector from the current point 
                                              // to its neighbor

                label faceI = pFaces[j];
                vector fnn = faceNs[ faceI ]; // face normal

                fnn = transform(I - curNormPP*curNormPP, fnn);

                // correction of the normal, otherwise realN = fnn
                scalar middi = mag(d2)/2.0;   // half distance to neighbor
                point midpo = curP + d2/2.0;  // midpoint
                plane pll(midpo, d2);           // plane perpendicular to the edge via midpoint
                point endNorm = midpo + fnn;  // end of the face normal projected on the inlet
                point projp = pll.nearestPoint(endNorm); // projection of endNorm onto pll plane

                vector realN = projp - midpo; // normal to the edge between current point
                                              // and its neighbor in the inlet plane facing
                                              // outside the fracture (mag(realN)!=1 !!!)

                scalar nw = 1.0 / middi;
                currentNorms[i] += nw * realN;

                currentNorms[i] = transform(I - curNormPP*curNormPP, currentNorms[i]);
                faceToPointSumWeights[i] += nw;
              }
            }

            syncTools::syncPointList(mesh, global_EdgePoints, sumWeights, plusEqOp<scalar>(), 0.0);
            syncTools::syncPointList(mesh, global_EdgePoints, currentNorms, plusEqOp<vector>(), vector::zero);
            syncTools::syncPointList(mesh, global_EdgePoints, faceToPointSumWeights, plusEqOp<scalar>(), 0.0);
                                                                                                      
            // normalization
            forAll(currentNorms, i)
              currentNorms[i] /= mag( currentNorms[i] );
          }

          vectorField movedNorms( NN, vector::zero );
          {
            pointField faceCs = faceCentres(movedPoints, llf);
            vectorField faceNs = faceNormals(movedPoints, llf);

            scalarList sumWeights( NN, 0.0 );

            forAll(nepe11, i)
            {
              label  curI = local_EdgePoints[i];
              const point& curP = movedPoints[curI];
              const labelList& pNeib = nepe11[i];

              label curIpp = local_pp_EdgePoints[i];
              const vector& curNormPP = ppPointNormals[ curIpp ];

              const labelList& pFaces = plistFaces[curI];
              forAll(pFaces, j)
              {
                label ind = pNeib[j];
                const point& neibP = movedPoints[ind];
                vector d2 = neibP - curP;     // vector from the current point to its neighbor

                label faceI = pFaces[j];
                vector fnn = faceNs[ faceI ]; // face normal
                fnn = transform(I - curNormPP*curNormPP, fnn);

                // correction of the normal, otherwise realN = fnn
                scalar middi = mag(d2)/2.0;   // half distance to neighbor
                point midpo = curP + d2/2.0;  // midpoint
                plane pll(midpo, d2);           // plane perpendicular to the edge via midpoint
                point endNorm = midpo + fnn;  // end of the face normal projected on the inlet
                point projp = pll.nearestPoint(endNorm); // projection of endNorm onto pll plane

                vector realN = projp - midpo; // normal to the edge between current point
                                              // and its neighbor in the inlet plane facing
                                              // outside the fracture (mag(realN)!=1 !!!)

                scalar nw = 1.0 / middi;
                movedNorms[i] += nw * realN;
                movedNorms[i] = transform(I - curNormPP*curNormPP, movedNorms[i]);
                faceToPointSumWeights[i] += nw;
              }
            }

            syncTools::syncPointList(mesh, global_EdgePoints, sumWeights, plusEqOp<scalar>(), 0.0);
            syncTools::syncPointList(mesh, global_EdgePoints, movedNorms, plusEqOp<vector>(), vector::zero);
            syncTools::syncPointList(mesh, global_EdgePoints, faceToPointSumWeights, plusEqOp<scalar>(), 0.0);
            // normalization
            forAll(movedNorms, i)
            {
              movedNorms[i] /= mag( movedNorms[i] );
            }
          }
          scalarField aa = mag(currentNorms ^ movedNorms);

          // TODO synchronize for processor
          label fixedEdgePoint = findMin(aa);
          //Pout<<nl<< " 111Patch name:  "<< pp.name() << nl
          //        << "  number of common points " << NN
          //        << "     fixedEdgePoint: " << fixedEdgePoint
          //        <<endl;
          // *************************************************************************
          //fixedPEdges.append( local_EdgePoints[fixedEdgePoint] );
          if(fixedEdgePoint>=0)
            fixedPEdges.append( fixedEdgePoint );
        //}
      }
    }
  }
}                


void Foam::
dissolMotionPointPatchVectorField::
neighborListEdge
(
  const labelList& list,
  const edgeList& eList,
  const labelListList& pointEdgesLi,
  labelListList& neLi
)
{
  forAll(list, i)
  {
    label ind = list[i];

    const labelList& lll = pointEdgesLi[ind];
    labelList nel;
    forAll(lll, ii)
    {
      const edge& ed = eList[ lll[ii] ];
      label edb = ed.start();
      
      // @TODO use otherVertex function instead
      if( edb!=ind && findIndex(list, edb) != -1 )
        nel.append(edb);

      label ede = ed.end();
      if( ede!=ind && findIndex(list, ede) != -1)
        nel.append(ede);
    }
    neLi.append(nel);
  }
}


void Foam::
dissolMotionPointPatchVectorField::
relaxPatchMesh(vectorField& pointMotion)
{
  pointField newPointsPos = this->patch().localPoints() + pointMotion;
  
  label patchID = this->patch().index();
  const polyMesh& mesh = this->internalField().mesh()();
  const polyBoundaryMesh& bMesh = mesh.boundaryMesh();
  const polyPatch& pPatch = refCast<const polyPatch>(bMesh[patchID]);
  const List<face>& llf = pPatch.localFaces();
  const labelListList& plistFaces = pPatch.pointFaces();
  const labelList& meshPoints = pPatch.meshPoints();
  
  int N = newPointsPos.size();
  
  /*
  // + + + + + + +  + + + + + + +  + + + + + + +  + + + + + + +  + + + + + + + 
  vectorFieldList updW(N);
  vectorField sumUpdW( N );

  // update weights
  vectorField faceNs = faceNormals(newPointsPos, llf);
  forAll(updW, i)
  {
    const labelList& pFaces = plistFaces[i];
    const vectorField& pw = surfWeights[i];
    vectorField& updPW = updW[i];

    updPW.setSize(pFaces.size());
    vector sumw = vector::zero;
    forAll(pFaces, j)
    {
      label faceI = pFaces[j];

      updPW[j] = 
              transform
              (
                I-faceNs[ faceI ]*faceNs[ faceI ],
                pw[j]
              );
      sumw += updPW[j];
      //Info<<" FF: "<< pw[j] <<"  NEWFF "<<updPW[j]<<"  faceN "<<faceNs[ faceI ]<<nl;
    }
    sumUpdW[i] = sumw;
    //Info<<nl;
  }
  //std::exit(0);

  syncTools::syncPointList(mesh, meshPoints, sumUpdW, plusEqOp<vector>(), vector::zero);

  forAll(updW, i)
  {
    vectorField& uw = updW[i];
    vector &suw = sumUpdW[i];
    forAll(uw, j)
    {
      vector &uuw = uw[j];
      for(int ii=0; ii<3; ii++)
      {
        if( mag(suw[ii])>SMALL )
        {
          uuw[ii] /= suw[ii];
        }
      }
    }
  }
  // + + + + + + +  + + + + + + +  + + + + + + +  + + + + + + +  + + + + + + + 
  */

  if(debug) 
  {
    Info << "    Relaxation loop: "<< nl;
  }
  vectorField pointNorm( N, vector::zero );
  scalarList faceToPointSumWeights( N, 0.0 );
  double displ_tol = 1.0;
  int itt = 0;
  while(displ_tol>rlxTol)
  {

    if(itt%q_norm_recalc==0){
      // set fields to zero
      pointNorm = vector::zero;
      faceToPointSumWeights = 0.0;
    }
      
    // calculate current face centers
    pointField faceCs = faceCentres(newPointsPos, llf);
    vectorField faceNs = faceNormals(newPointsPos, llf);

    // create a displacement field for points
    vectorField displacement(N, vector::zero);
    scalarField tol(N, 0.0);
    
    
    forAll(newPointsPos, i)
    {
      point& curP = newPointsPos[i];
      const labelList& pFaces = plistFaces[i];
      const vectorField& pw = surfWeights[i];
      //const vectorField& puw = updW[i];

      forAll(pFaces, j)
      {
        label faceI = pFaces[j];
        point& faceC = faceCs[faceI];

        vector d = faceC - curP;

        scalar mag_d = mag(d);
        
        /*
        vector newW = 
                transform
                (
                  I-faceNs[ faceI ]*faceNs[ faceI ],
                  pw[j]
                );
        */

        vector disp;
        //disp.x() = d.x() * puw[j].x(); //* newW.x();
        //disp.y() = d.y() * puw[j].y() ; //* newW.y();
        //disp.z() = d.z() * puw[j].z() ; //* newW.z();
        disp.x() = d.x() * pw[j].x();
        disp.y() = d.y() * pw[j].y();
        disp.z() = d.z() * pw[j].z();
        displacement[i] += disp;

        //tol[i] += mag_d * mag(pw[j]);
        tol[i] += mag( disp );

        // this is for normal
        if(itt%q_norm_recalc==0)
        {
          scalar nw = 1.0 / mag_d;
          pointNorm[i] += nw * faceNs[ faceI ];
          faceToPointSumWeights[i] += nw;
        }
      }

      // TODO stick to cyclic boundary
    }

    // synchronizing over the cyclic and processor boundaries
    syncTools::syncPointList(mesh, meshPoints, displacement, plusEqOp<vector>(), vector::zero);
    syncTools::syncPointList(mesh, meshPoints, tol, plusEqOp<scalar>(), 0.0);

    if(itt%q_norm_recalc==0){
      syncTools::syncPointList(mesh, meshPoints, pointNorm, plusEqOp<vector>(), vector::zero);
      syncTools::syncPointList(mesh, meshPoints, faceToPointSumWeights, plusEqOp<scalar>(), 0.0);

      // normalization
      forAll(pointNorm, pointi)
      {
        pointNorm[pointi] /= faceToPointSumWeights[pointi];
        pointNorm[pointi] /= mag(pointNorm[pointi]);
      }
    }

    // apply edge boundary for relxation (edges between patches are fixed)
    forAll(fixedPoints, i)
    {
      label pointI = fixedPoints[i];
      displacement[ pointI ] = vector::zero;
    }
    
    fixPinnedPoints(displacement);
    
    vectorField projectedDisplacement = transform(I - pointNorm*pointNorm, displacement);

    // TODO stick to cyclic boundary
    /*
    forAll(pinnedPoints, ii)
    {
      label ind = pinnedPoints[ii];
      projectedDisplacement[ind] = 
              transform
              (
                I-pinnedPointsNorm[ii]*pinnedPointsNorm[ii],
                projectedDisplacement[ind]
              );
    }
    */
    
    scalar factor = (itt%q_2==0) ? k_2 : k_1;

    vectorField finalDisplacement = factor * projectedDisplacement;

    newPointsPos += finalDisplacement;

    displ_tol = gAverage( mag(finalDisplacement/factor)/tol );

    if(debug) 
    {
      //Info <<"    "<< itt << "  Displ_tol: "<< displ_tol << nl;
      //if(itt==0)
      //  break;
    }
    
    if(itt%1000==0)
    {
      Info << "  " << this->patch().name() << "  rlx iter " << itt
           << "  tolerance: " << displ_tol << endl;
    }
    
    itt++;
  }
  Info << nl << "  " << this->patch().name() << "  rlx converged in " << itt 
       << " iterations. Tolerance: " << displ_tol<< nl << endl;

  pointMotion = (newPointsPos - this->patch().localPoints());
  
  fixPinnedPoints(pointMotion);
  
  scalarField cc(pointMotion.size(), 1.0);
  syncTools::syncPointList(mesh, meshPoints, pointMotion, plusEqOp<vector>(), vector::zero);
  syncTools::syncPointList(mesh, meshPoints, cc, plusEqOp<scalar>(), 0.0);

  forAll(pointMotion, i){
    pointMotion[i] /= cc[i];
  }

  /*
  Pout<<nl<< "BBB1:  " << this->patch().localPoints()[0] << "   " << pointMotion[0] <<endl;
  forAll(pinnedPoints, ii)
  {
    label ind = pinnedPoints[ii];
    
    //if(  this->patch().localPoints()[ind].x() < 1
    //  && this->patch().localPoints()[ind].x() > 0
    //)
    if(ind==4884 || ind==5139)
    {
      Pout<<ind<<"  "
              << this->patch().localPoints()[ind] << "  "
              << pointMotion[ind]<<endl;
    }
  }
   */
  //std::exit(0);
}



void Foam::
dissolMotionPointPatchVectorField::
relaxPatchMeshPoints(vectorField& pointMotion)
{
  pointField newPointsPos = this->patch().localPoints() + pointMotion;
  
  label patchID = this->patch().index();
  const polyMesh& mesh = this->internalField().mesh()();
  const polyBoundaryMesh& bMesh = mesh.boundaryMesh();
  const polyPatch& pPatch = refCast<const polyPatch>(bMesh[patchID]);
  const List<face>& llf = pPatch.localFaces();
  const labelListList& plistFaces = pPatch.pointFaces();
  //const labelListList& flistFaces = pPatch.faceFaces();
  const labelList& meshPoints = pPatch.meshPoints();
  
  int N = newPointsPos.size();
  
  //const pointField& boundaryPoints = pPatch.localPoints();
  
  //const pointField&  faceCs = pPatch.faceCentres();
  //const vectorField& faceNs = pPatch.faceNormals();
  
  if(debug) 
  {
    Info << "    Relaxation loop (only points): "<< nl;
  }
  vectorField pointNorm( N, vector::zero );
  scalarList faceToPointSumWeights( N, 0.0 );
  double displ_tol = 1.0;
  int itt = 0;
  while(displ_tol>rlxTol)
  {

    if(itt%q_norm_recalc==0){
      // set fields to zero
      pointNorm = vector::zero;
      faceToPointSumWeights = 0.0;
    }
      
    // calculate current face centers
    pointField faceCs = faceCentres(newPointsPos, llf);
    vectorField faceNs = faceNormals(newPointsPos, llf);

    // create a displacement field for points
    vectorField displacement(N, vector::zero);
    scalarField tol(N, 0.0);
    
    // !!!!!!!!!!!!!!!!!!!!!!!! @HERE
    forAll(newPointsPos, i)
    {
      point& curP = newPointsPos[i];
      const labelList& pFaces = plistFaces[i];
      const vectorField& pw = surfWeights[i];
      //const vectorField& puw = updW[i];

      // recalculate normals
      if(itt%q_norm_recalc==0)
      {
        forAll(pFaces, j)
        {
          label faceI = pFaces[j];
          point& faceC = faceCs[faceI];

          vector d = faceC - curP;

          scalar mag_d = mag(d);

          scalar nw = 1.0 / mag_d;
          pointNorm[i] += nw * faceNs[ faceI ];
          faceToPointSumWeights[i] += nw;
        }
      }

      
      forAll(pFaces, j)
      {
        label faceI = pFaces[j];
        point& faceC = faceCs[faceI];

        vector d = faceC - curP;

        scalar mag_d = mag(d);
        
        /*
        vector newW = 
                transform
                (
                  I-faceNs[ faceI ]*faceNs[ faceI ],
                  pw[j]
                );
        */

        vector disp;
        //disp.x() = d.x() * puw[j].x(); //* newW.x();
        //disp.y() = d.y() * puw[j].y() ; //* newW.y();
        //disp.z() = d.z() * puw[j].z() ; //* newW.z();
        disp.x() = d.x() * pw[j].x();
        disp.y() = d.y() * pw[j].y();
        disp.z() = d.z() * pw[j].z();
        displacement[i] += disp;

        //tol[i] += mag_d * mag(pw[j]);
        tol[i] += mag( disp );

        // this is for normal
        if(itt%q_norm_recalc==0)
        {
          scalar nw = 1.0 / mag_d;
          pointNorm[i] += nw * faceNs[ faceI ];
          faceToPointSumWeights[i] += nw;
        }
      }
      
      
      // TODO stick to cyclic boundary
    }

    // synchronizing over the cyclic and processor boundaries
    syncTools::syncPointList(mesh, meshPoints, displacement, plusEqOp<vector>(), vector::zero);
    syncTools::syncPointList(mesh, meshPoints, tol, plusEqOp<scalar>(), 0.0);

    if(itt%q_norm_recalc==0){
      syncTools::syncPointList(mesh, meshPoints, pointNorm, plusEqOp<vector>(), vector::zero);
      syncTools::syncPointList(mesh, meshPoints, faceToPointSumWeights, plusEqOp<scalar>(), 0.0);

      // normalization
      forAll(pointNorm, pointi)
      {
        pointNorm[pointi] /= faceToPointSumWeights[pointi];
        pointNorm[pointi] /= mag(pointNorm[pointi]);
      }
    }

    // apply edge boundary for relxation (edges between patches are fixed)
    forAll(fixedPoints, i)
    {
      label pointI = fixedPoints[i];
      displacement[ pointI ] = vector::zero;
    }
    
    fixPinnedPoints(displacement);
    
    vectorField projectedDisplacement = transform(I - pointNorm*pointNorm, displacement);

    // TODO stick to cyclic boundary
    /*
    forAll(pinnedPoints, ii)
    {
      label ind = pinnedPoints[ii];
      projectedDisplacement[ind] = 
              transform
              (
                I-pinnedPointsNorm[ii]*pinnedPointsNorm[ii],
                projectedDisplacement[ind]
              );
    }
    */
    
    scalar factor = (itt%q_2==0) ? k_2 : k_1;

    vectorField finalDisplacement = factor * projectedDisplacement;

    newPointsPos += finalDisplacement;

    displ_tol = gAverage( mag(finalDisplacement/factor)/tol );

    if(debug) 
    {
      //Info <<"    "<< itt << "  Displ_tol: "<< displ_tol << nl;
      //if(itt==0)
      //  break;
    }
    
    if(itt%1000==0)
    {
      Info << "  " << this->patch().name() << "  rlx iter " << itt
           << "  tolerance: " << displ_tol << endl;
    }
    
    itt++;
  }
  Info << nl << "  " << this->patch().name() << "  rlx converged in " << itt 
       << " iterations. Tolerance: " << displ_tol<< nl << endl;

  pointMotion = (newPointsPos - this->patch().localPoints());
  
  fixPinnedPoints(pointMotion);
  
  scalarField cc(pointMotion.size(), 1.0);
  syncTools::syncPointList(mesh, meshPoints, pointMotion, plusEqOp<vector>(), vector::zero);
  syncTools::syncPointList(mesh, meshPoints, cc, plusEqOp<scalar>(), 0.0);

  forAll(pointMotion, i){
    pointMotion[i] /= cc[i];
  }

  /*
  Pout<<nl<< "BBB1:  " << this->patch().localPoints()[0] << "   " << pointMotion[0] <<endl;
  forAll(pinnedPoints, ii)
  {
    label ind = pinnedPoints[ii];
    
    //if(  this->patch().localPoints()[ind].x() < 1
    //  && this->patch().localPoints()[ind].x() > 0
    //)
    if(ind==4884 || ind==5139)
    {
      Pout<<ind<<"  "
              << this->patch().localPoints()[ind] << "  "
              << pointMotion[ind]<<endl;
    }
  }
   */
  //std::exit(0);
}






void Foam::
dissolMotionPointPatchVectorField::
make_lists_and_normals()
{
  const Foam::Time& time = this->db().time();
  Foam::Time timeTmp(Foam::Time::controlDictName, time.rootPath(), time.caseName());
  //Foam::instantList timeDirs = Foam::timeSelector::select(time.times());
  Foam::instantList timeDirs = time.times();
  timeTmp.setTime(timeDirs[0], 0);
  if( timeTmp.timeName()=="constant" )
  {
    timeTmp.setTime(timeDirs[1], 0);
  }
  
/*
  // in case 0 time does not exist
  if( timeTmp.timeName()!="0" ){
    SeriousErrorIn("meshRelax")
            <<"There is no 0 time directory! "
            "If you run in parallel please check a decomposition."
            <<exit(FatalError);
  }
*/
  
  Foam::fvMesh meshTmp
  (
    Foam::IOobject
    (
      Foam::fvMesh::defaultRegion,
      timeTmp.timeName(),
      timeTmp,
      Foam::IOobject::MUST_READ
    )
  );
  
  label patchID = this->patch().index();
  //const polyMesh& mesh = this->internalField().mesh()();
  const polyBoundaryMesh& bMesh = meshTmp.boundaryMesh();
  const polyPatch& pPatch = refCast<const polyPatch>(bMesh[patchID]);
  //const List<face>& llf = pPatch.localFaces();
  //const labelListList& plistFaces = pPatch.pointFaces();
  const labelList& meshPoints = pPatch.meshPoints();
  
  
  const pointField& curPP  = pPatch.localPoints();
  const labelListList& curPointEdges = pPatch.pointEdges();
  const edgeList& ee = pPatch.edges();
 
  
  // list fixed boundary edges
  //labelList pinnedPointsLocal;
  //vectorField pinnedPointsNormLocal;
  //labelList fixedPointsLocal;

  //labelList pinnedPoints;
  //vectorField pinnedPointsNorm;
  forAll(bMesh, patchi)
  {
    if ( (patchi != patchID) ) //bMesh[patchi].size() &&
    {
      //Info<< "patch type:  "<< bMesh[patchi].type()<<endl;
      // TODO redesign this if condition
      if 
      (
        isA<cyclicPolyPatch>(bMesh[patchi]) 
        ||
        isA<symmetryPolyPatch>(bMesh[patchi])
      )
      {
        const polyPatch& cpp = refCast<const polyPatch>(bMesh[patchi]);
        // TODO
        const labelList& ppMeshPoints = cpp.meshPoints();
        const vectorField& ppPointNormals = cpp.pointNormals(); //1

        labelList local_EdgePoints, global_EdgePoints;
        labelList local_pp_EdgePoints;

        commonPoints
        (
          meshPoints,
          ppMeshPoints,
          local_EdgePoints,
          global_EdgePoints,
          local_pp_EdgePoints
        );

        pinnedPoints.append(local_EdgePoints);
        vectorField locNormals(local_EdgePoints.size(), vector::zero);
        forAll(local_EdgePoints, i)
        {
          locNormals[i] = ppPointNormals[ local_pp_EdgePoints[i] ];
        }
        pinnedPointsNorm.append(locNormals);
      }
      else if (isA<processorPolyPatch>(bMesh[patchi]))
      {
        // skip
      }
      else
      {
        const polyPatch& pp = refCast<const polyPatch>(bMesh[patchi]);
        const labelList& ppMeshPoints = pp.meshPoints();
        const vectorField& ppPointNormals = pp.pointNormals(); //2

        labelList local_EdgePoints, global_EdgePoints;
        labelList local_pp_EdgePoints;

        commonPoints
        (
          meshPoints,
          ppMeshPoints,
          local_EdgePoints,
          global_EdgePoints,
          local_pp_EdgePoints
        );

        fixedPoints.append(local_EdgePoints);
        
        vectorField locNormals(local_EdgePoints.size(), vector::zero);
        forAll(local_EdgePoints, i)
        {
          locNormals[i] = ppPointNormals[ local_pp_EdgePoints[i] ];
        }
        
        fixedPointNorms.append(locNormals);
        globalFixedPoints.append(global_EdgePoints);
        
        // * * * * * * * * * * * * 

        labelListList nepeLocal;
        neighborListEdge(local_EdgePoints, ee, curPointEdges, nepeLocal);

        label NN = local_EdgePoints.size();

        scalarFieldList weights( NN );
        scalarField sumWeights( NN, 0.0 );

        forAll(nepeLocal, i)
        {
          label  curI = local_EdgePoints[i];
          const point& curP = curPP[curI];
          const labelList& pNeib = nepeLocal[i];

          scalarField& pw = weights[i];
          pw.setSize(pNeib.size());

          scalar sumw = 0.0;
          forAll(pNeib, ii)
          {
            label ind = pNeib[ii];
            const point& neibP = curPP[ind];
            vector d2 = neibP - curP;
            scalar magd2 = mag(d2);

            scalar w = 0.0;
            if( magd2>SMALL) w = 1.0 / magd2;

            pw[ii] = w;
            sumw += pw[ii];
          }
          sumWeights[i] = sumw;
        }
        
        syncTools::syncPointList
        (
          meshTmp,
          global_EdgePoints,
          sumWeights,
          plusEqOp<scalar>(),
          0.0
        );

        forAll(weights, i)
        {
          scalarField& pw = weights[i];
          scalar &sw = sumWeights[i];
          forAll(pw, j)
            if( mag(sw)>SMALL ) pw[j] /= sw;
        }

        rlxEdgeWeights.append(weights);
        nepe.append(nepeLocal);
        // * * * * * * * * * * * * 
        
      }
    }
  }

  //pinnedPoints = pinnedPointsLocal;
  //pinnedPointsNorm = pinnedPointsNormLocal;
  //fixedPoints = fixedPointsLocal;

/*  
  forAll(pinnedPoints, i)
    Pout<<"pinnedPoints: " << pinnedPoints[i] << "  " 
          << pinnedPointsNorm[i] << "  " 
          << curPP[pinnedPoints[i]]<< "  " << surfWeights[pinnedPoints[i]] << endl;
 */
  
  this->setListsUpdated(true);
}


void Foam::
dissolMotionPointPatchVectorField::
calc_weights_surface()
{
  const Foam::Time& time = this->db().time();
  Foam::Time timeTmp(Foam::Time::controlDictName, time.rootPath(), time.caseName());
  //Foam::instantList timeDirs = Foam::timeSelector::select(time.times());
  Foam::instantList timeDirs = time.times();
  timeTmp.setTime(timeDirs[0], 0);
  if( timeTmp.timeName()=="constant" )
  {
    timeTmp.setTime(timeDirs[1], 0);
  }
  
/*
  // in case 0 time does not exist
  if( timeTmp.timeName()!="0" ){
    SeriousErrorIn("meshRelax")
            <<"There is no 0 time directory! "
            "If you run in parallel please check a decomposition."
            <<exit(FatalError);
  }
*/
  
  Foam::fvMesh meshTmp
  (
    Foam::IOobject
    (
      Foam::fvMesh::defaultRegion,
      timeTmp.timeName(),
      timeTmp,
      Foam::IOobject::MUST_READ
    )
  );
  
  label patchID = this->patch().index();
  //const polyMesh& mesh = this->internalField().mesh()();
  const polyBoundaryMesh& bMesh = meshTmp.boundaryMesh();
  const polyPatch& pp = refCast<const polyPatch>(bMesh[patchID]);
  //const List<face>& flist = pp.localFaces();
  const labelListList& plistFaces = pp.pointFaces();
  const labelList& meshPoints = meshTmp.boundaryMesh()[patchID].meshPoints();
  
  const pointField& boundaryPoints = pp.localPoints();
  const pointField& faceCs = pp.faceCentres();
  
  vectorFieldList weights( boundaryPoints.size() );
  vectorField sumWeights( boundaryPoints.size() );

  forAll(boundaryPoints, i)
  {
    point curP = boundaryPoints[i];

    const labelList& pFaces = plistFaces[i];

    vectorField& pw = weights[i];
    pw.setSize(pFaces.size());

    vector sumw = vector::zero;
    forAll(pFaces, j)
    {
      label faceI = pFaces[j];
      point faceC = faceCs[faceI];

      vector d = faceC - curP;

      vector w = vector::zero;

      for(int ii=0; ii<3; ii++){
        if(mag(d[ii])>SMALL)
        {
          w[ii] = 1.0 / mag( d[ii] );
        }
        else
        {
          w[ii] = GREAT;
        }
      }

      pw[j] = w;
      sumw += pw[j];
    }
    // sum of all distances to face centres
    sumWeights[i] = sumw;
  }

  syncTools::syncPointList( meshTmp, meshPoints, sumWeights, plusEqOp<vector>(), vector::zero);

  scalarField cc(boundaryPoints.size(), 1.0);
  syncTools::syncPointList(meshTmp, meshPoints, cc, plusEqOp<scalar>(), 0.0);

  //const vectorField& faceNs = pp.faceNormals();

  // TODO!!! Re weights
  forAll(weights, i)
  {
    vectorField& pw = weights[i];
    vector &sw = sumWeights[i];
    //const labelList& pFaces = plistFaces[i];
    forAll(pw, j)
    {
      vector &ppw = pw[j];
      for(int ii=0; ii<3; ii++)
      {
        if( mag(sw[ii])>SMALL )
        {
          ppw[ii] /= sw[ii];
        }
        else
        {
          ppw[ii] = GREAT;
        }
        //ppw[ii] = 1.0 / cc[i];
        //if(ii==1) ppw[ii] = 0.0;
      }
      
      //label faceI = pFaces[j];

      /*
      pw[j] = 
              transform
              (
                I-faceNs[ faceI ]*faceNs[ faceI ],
                pw[j]
              );
       */
      //pw[j] /= mag(pw[j]);
    }
    //Info<<" i "<<i<<"  pw "<<pw<<nl;
  }
  
  //std::exit(0);

  surfWeights = weights;
  this->setWeightsUpdated(true);
}


void Foam::
dissolMotionPointPatchVectorField::
calc_point_weights_surface()
{
  label patchID = this->patch().index();
  const polyMesh& mesh = this->internalField().mesh()();
  
  const polyBoundaryMesh& bMesh = mesh.boundaryMesh();
  const polyPatch& pPatch = refCast<const polyPatch>(bMesh[patchID]);
  //const List<face>& llf = pPatch.localFaces();
  //const labelListList& plistFaces = pPatch.pointFaces();
  //const labelListList& flistFaces = pPatch.faceFaces();
  const labelListList& plistEdges = pPatch.pointEdges();
  const labelList& meshPoints = pPatch.meshPoints();
  
  const edgeList& edges = pPatch.edges();

  const pointField& boundaryPoints = pPatch.localPoints();
  
  //const pointField&  faceCs = pPatch.faceCentres();
  //const vectorField& faceNs = pPatch.faceNormals();
  
  vectorFieldList weights( boundaryPoints.size() );
  vectorField sumWeights( boundaryPoints.size() );

  scalarField cc(boundaryPoints.size(), 1.0);
  syncTools::syncPointList(mesh, meshPoints, cc, plusEqOp<scalar>(), 0.0);
  
  forAll(boundaryPoints, i)
  {
    point curP = boundaryPoints[i];

    const labelList& pEdges = plistEdges[i];

    vectorField& pw = weights[i];
    pw.setSize(pEdges.size());

    vector sumw = vector::zero;
    forAll(pEdges, j)
    {
      label edgeI = pEdges[j];
      
      label nPoint = edges[edgeI].otherVertex(i);
      
      point neiP = boundaryPoints[nPoint];

      vector d = neiP - curP;

      vector w = vector::zero;

      for(int ii=0; ii<3; ii++){
        if(mag(d[ii])>SMALL){
          w[ii] = 1.0 / mag( d[ii] );
        }
        else{
          //w[ii] = 1.0;
          w[ii] = GREAT;
        }
      }

      pw[j] = w / cc[nPoint];
      sumw += pw[j];
    }
    // sum of all distances to face centres
    sumWeights[i] = sumw;
  }

  syncTools::syncPointList( mesh, meshPoints, sumWeights, plusEqOp<vector>(), vector::zero);
  

  // @TODO!!! Review weights
  forAll(weights, i)
  {
    vectorField& pw = weights[i];
    vector &sw = sumWeights[i];
    //const labelList& pEdges = plistEdges[i];
    forAll(pw, j)
    {
      vector &ppw = pw[j];
      for(int ii=0; ii<3; ii++)
      {
        if( mag(sw[ii])>SMALL )
        {
          ppw[ii] /= sw[ii];
        }
        else
        {
          ppw[ii] *= GREAT;
        }
        //ppw[ii] *= 1.0 / cc[i];
      }
      
      /*
      pw[j] = 
              transform
              (
                I-faceNs[ faceI ]*faceNs[ faceI ],
                pw[j]
              );
       */
      //pw[j] /= mag(pw[j]);
    }
    //Info<<" i "<<i<<"  pw "<<pw<<nl;
  }
  
  //std::exit(0);

  surfPointWeights = weights;
}



Foam::vectorField Foam::
dissolMotionPointPatchVectorField::
faceNormals(const pointField& points, const List<face>& flist) const
{
  vectorField fn( flist.size() );
  forAll(fn, facei)
  {
    fn[facei]  = flist[facei].normal(points);
    fn[facei] /= mag(fn[facei]) + VSMALL;
  }
  return fn;
}
Foam::pointField Foam::
dissolMotionPointPatchVectorField::
faceCentres(const pointField& points, const List<face>& flist) const
{
  pointField fc( flist.size() );
  forAll(fc, facei)
  {
    fc[facei] = flist[facei].centre(points);
  }
  return fc;
}

void Foam::
dissolMotionPointPatchVectorField::
fixPinnedPoints( vectorField& pointMotion )
{
  forAll(pinnedPoints, ppi)
  {
    label pointI = pinnedPoints[ppi];
    pointMotion[pointI] = 
            transform
            (
              I-pinnedPointsNorm[ppi]*pinnedPointsNorm[ppi],
              pointMotion[pointI]
            );
  }
}

void Foam::
dissolMotionPointPatchVectorField::
fixCommonNeighborPatchPoints( vectorField& pointMotion )
{
  label patchID = this->patch().index();
  const polyMesh& mesh = this->internalField().mesh()();

  // current patch geometry
  const pointField& curPP  = mesh.boundaryMesh()[patchID].localPoints();
  const labelList& curMeshPoints = mesh.boundaryMesh()[patchID].meshPoints();
  
  const polyBoundaryMesh& bMesh = mesh.boundaryMesh();
  
  fixPinnedPoints(pointMotion);
  
  forAll(bMesh, patchi)
  {
    /*
    Info<< "WWWWW patch type:  "<< bMesh[patchi].type() 
            << "  " << bMesh[patchi].name()
            <<endl;
     */
    if (bMesh[patchi].size() && (patchi != patchID) )
    {
      if(
        isA<cyclicPolyPatch>(bMesh[patchi]) 
        ||
        isA<symmetryPolyPatch>(bMesh[patchi])
        )
      {
        if(debug) 
        {
          Info << "  This is cyclic or symmetry patch, so it is skipped"<<nl;
        }

        /*
        const polyPatch& cpp = 
            refCast<const polyPatch>(bMesh[patchi]);
        // TODO
        
        const labelList& ppMeshPoints = cpp.meshPoints();

        labelList local_EdgePoints, global_EdgePoints;
        labelList local_pp_EdgePoints;
    
        commonPoints
        (
          curMeshPoints,
          ppMeshPoints,
          local_EdgePoints,
          global_EdgePoints,
          local_pp_EdgePoints
        );
        
        //vectorField
        forAll(local_EdgePoints, i)
        {
          label pointI = local_EdgePoints[i];
          
          vector nrm(vector::zero);
          label iii = findIndex(pinnedPoints, pointI);
          if( iii != -1 )
          {
            nrm = pinnedPointsNorm[iii];
          }

          pointMotion[pointI] = 
                  transform
                  (
                    I-nrm*nrm,
                    pointMotion[pointI]
                  );
        }
        */
        // skip
        
      }
      else if (isA<processorPolyPatch>(bMesh[patchi]))
      {
        if(debug) 
        {
          Info << "  This is processor patch, so it is skipped"<<nl;
        }
        // skip
      }
      else
      {
        const polyPatch& pp = refCast<const polyPatch>(bMesh[patchi]);
        const pointField& pfPP = pp.localPoints();
        const labelList& ppMeshPoints = pp.meshPoints();
        const vectorField& ppPointNormals = pp.pointNormals();  // 2

        labelList local_EdgePoints, global_EdgePoints;
        labelList local_pp_EdgePoints;
    
        commonPoints
        (
          curMeshPoints,
          ppMeshPoints,
          local_EdgePoints,
          global_EdgePoints,
          local_pp_EdgePoints
        );
        
        /*
        if(local_EdgePoints.size()>0)
        {
          Pout<< "BBBBBpatch type: "<< bMesh[patchi].type()
                  << " name: "<< bMesh[patchi].name()
                  <<" KKK: "<< pointMotion[local_EdgePoints[0]]
                  <<"  "<< ppPointNormals[ local_pp_EdgePoints[0] ]
                  << endl;
        }
        */
        
        if(debug) 
        {
          Info << "  This is "<<bMesh[patchi].name()
                  <<" patch, of type "<<bMesh[patchi].type()<<nl;
          Info << "  Number of common points is " << local_EdgePoints.size()<<nl;
        }
        
        
        //vectorField
        forAll(local_EdgePoints, i)
        {
          label pointI = local_EdgePoints[i];
          point Ap = curPP[pointI] + pointMotion[pointI];
          label ppI = local_pp_EdgePoints[i];
          
          // plane perpendicular to the edge via midpoint
          plane pll(pfPP[ ppI ], ppPointNormals[ ppI ]);

          // projection of endNorm onto pll plane
          point projp = pll.nearestPoint(Ap);

          vector AA = projp - curPP[pointI];

          scalar cosa = 0.0;
          
          if( mag(pointMotion[pointI]) > SMALL && mag(AA) > SMALL )
          {
            cosa = 
              (pointMotion[pointI] & AA) / (mag(pointMotion[pointI])*mag(AA));
          }
          
          if( mag(cosa)>SMALL && mag(AA) > SMALL)
          {
            scalar L = mag(pointMotion[pointI]) / cosa;
            pointMotion[pointI] = L * AA / mag(AA);
          }
        }
        
        /*
        if(local_EdgePoints.size()>0)
        {
          Pout<< "AAAAApatch type: "<< bMesh[patchi].type()
                  << " name: "<< bMesh[patchi].name()
                  <<" KKK: "<< pointMotion[local_EdgePoints[0]]
                  <<"  "<< ppPointNormals[ local_pp_EdgePoints[0] ]
                  << endl;
        }
         */
        
      
      }
      
      /*
      Info<< "patch type: "<< bMesh[patchi].type()
              << " name: "<< bMesh[patchi].name()
              <<" KKK: "<< pointMotion[local_EdgePoints[0]] << endl;
       */
      
      
      
    }
  }
  
}

void Foam::
dissolMotionPointPatchVectorField::
commonPoints
(
  const labelList& list1, 
  const labelList& list2, 
  labelList& localList1,
  labelList& globalList,
  labelList& localList2
)
{
  forAll(list1, loclLabel1)
  {
    label globLabel1 = list1[loclLabel1];
    label loclLabel2 = findIndex(list2, globLabel1);
    if( loclLabel2 != -1 )
    {
      localList1.append( loclLabel1 );
      globalList.append( globLabel1 );
      localList2.append( loclLabel2 );
    }
  }
}

void Foam::
dissolMotionPointPatchVectorField::
getPointMotion( vectorField& pointMotion )
{
  label patchID = this->patch().index();
  
  const polyMesh& mesh = this->internalField().mesh()();

  const IOdictionary& IOd
        = this->db().lookupObject<IOdictionary>("transportProperties");
  scalar lR =  (new dimensionedScalar(IOd.lookup("lR")))->value();

  const volScalarField& C = 
    this->db().objectRegistry::lookupObject<volScalarField>("C");
  
  //- set-up interpolator
  const fvMesh& fvmesh_ = refCast<const fvMesh>(mesh);
  coupledPatchInterpolation patchInterpolator
  (
    mesh.boundaryMesh()[patchID], fvmesh_
  );
  
  scalarField pointCface = -C.boundaryField()[patchID].snGrad();
  vectorField pointNface = mesh.boundaryMesh()[patchID].faceNormals();
  
  scalarField motionC    = patchInterpolator.faceToPointInterpolate(pointCface);
  vectorField motionN    = patchInterpolator.faceToPointInterpolate(pointNface);
  
  // normalize point normals to 1
  forAll(motionN, ii) motionN[ii]/=mag(motionN[ii]);
  
  // set the velocity for the point motion
  pointMotion = lR * motionC * motionN;
}

void Foam::dissolMotionPointPatchVectorField::write
(
  Ostream& os
) const
{
  fixedValuePointPatchField<vector>::write(os);

  os.writeKeyword("scaleMotion")<< scaleMotion << token::END_STATEMENT << nl;
  os.writeKeyword("rlxTol")<< rlxTol << token::END_STATEMENT << nl;

  os.writeKeyword("surfaceRlx")<< surfaceRlx << token::END_STATEMENT << nl;
  os.writeKeyword("q_norm_recalc")<<q_norm_recalc<<token::END_STATEMENT<<nl;
  os.writeKeyword("k_1")<< k_1 << token::END_STATEMENT << nl;
  os.writeKeyword("k_2")<< k_2 << token::END_STATEMENT << nl;
  os.writeKeyword("q_2")<< q_2 << token::END_STATEMENT << nl;
  os.writeKeyword("edgeRlx")<< edgeRlx << token::END_STATEMENT << nl;
  os.writeKeyword("q_norm_recalc_edge")<<q_norm_recalc_edge
          <<token::END_STATEMENT<<nl;
  os.writeKeyword("k_1edge")<< k_1edge << token::END_STATEMENT << nl;
  os.writeKeyword("k_2edge")<< k_2edge << token::END_STATEMENT << nl;
  os.writeKeyword("q_2edge")<< q_2edge << token::END_STATEMENT << nl;

  os.writeKeyword("pinnedPoint")<< pinnedPoint << token::END_STATEMENT << nl;
}

namespace Foam
{
  makePointPatchTypeField
  (
    pointPatchVectorField,
    dissolMotionPointPatchVectorField
  );
}


// ************************************************************************* //
