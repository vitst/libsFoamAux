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
    fixedValuePointPatchField<vector>(p, iF)
{
  if(debug) 
  {
    Info << "dissolMotionPointPatchVectorField constructor 0"<<endl;
  }
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
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper)
    //timeSeries_(ptf.timeSeries_)
{
  if(debug) 
  {
    Info << "dissolMotionPointPatchVectorField constructor 1"<<endl;
  }
  
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
    //timeSeries_(dict)
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
  
  if(debug) 
  {
    Info << "dissolMotionPointPatchVectorField constructor 2 "
            "calc_weights_surface"<<nl;
  }
  calc_weights_surface();
  
  if(debug) 
  {
    Info << "dissolMotionPointPatchVectorField constructor 2 "
            "make_lists_and_normals"<<nl;
  }
  make_lists_and_normals();
  
  /*
  for(int i=0; i<5; i++)
  {
    Info<<i<<" "<<surfWeights[i]<<"  "<<pinnedPointsNorm[i]<<nl;
  }
  std::exit(0);
  */
  
}

Foam::
dissolMotionPointPatchVectorField::
dissolMotionPointPatchVectorField
(
    const dissolMotionPointPatchVectorField& ptf
)
:
    fixedValuePointPatchField<vector>(ptf)
    //timeSeries_(ptf.timeSeries_)
{
  if(debug)
  {
    Info << "dissolMotionPointPatchVectorField constructor 3"<<endl;
  }
}


Foam::
dissolMotionPointPatchVectorField::
dissolMotionPointPatchVectorField
(
    const dissolMotionPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF)
    //timeSeries_(ptf.timeSeries_)
{
  if(debug) {
    Info << "dissolMotionPointPatchVectorField constructor 4"<<endl;
  }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dissolMotionPointPatchVectorField::updateCoeffs()
{
  if(debug) 
  {
    Info << "dissolMotionPointPatchVectorField::updateCoeffs()"<<endl;
  }
  
  if (this->updated())
  {
    return;
  }
  
  if( this->db().foundObject<IOdictionary>("transportProperties") )
  {
    //label patchID = this->patch().index();
    const scalar dt = this->db().time().deltaTValue();

    //const polyMesh& mesh = this->internalField().mesh()();

    vectorField pointMotion;
    getPointMotion(pointMotion);
    pointMotion *= dt;
    
    /*
    Info<<"PPPoint motion: "<<nl;
    for(int i = 0; i<10; i++)
    {
      Info<< pointMotion[i] << nl;
    }
    */

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // TODO the rest of relaxation should be implemented
    Info<<"fixCommonNeighborPatchPoints"<<nl;
    fixCommonNeighborPatchPoints(pointMotion);

    Info<<nl<<"relaxEdges"<<nl;
    relaxEdges(pointMotion);

    Info<<nl<<"relaxPatchMesh"<<nl;
    relaxPatchMesh(pointMotion);
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

    fixedValuePointPatchField<vector>::updateCoeffs();
  }
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
  const labelList& curMeshPoints = curPatch.meshPoints();
  const labelListList& curPointEdges = curPatch.pointEdges();
  const edgeList& ee = curPatch.edges();
  
  const List<face>& llf = curPatch.localFaces();
  const labelListList& plistFaces = curPatch.pointFaces();

  const polyBoundaryMesh& bMesh = mesh.boundaryMesh();
  
  // TODO global variables
  int q_2edge = 1;
  double k_1edge = 1.0, k_2edge = 1.2;
  int q_edge_norm_recalc = 5;
  
  forAll(bMesh, patchi)
  {
    if ( (patchi != patchID) ) //bMesh[patchi].size() && 
    {
      if (isA<cyclicPolyPatch>(bMesh[patchi]))
      {
        //const cyclicPolyPatch& cpp = 
        //    refCast<const cyclicPolyPatch>(bMesh[patchi]);
        // TODO
      }
      else if (isA<symmetryPolyPatch>(bMesh[patchi]))
      {
        //const symmetryPolyPatch& spp = 
        //    refCast<const symmetryPolyPatch>(bMesh[patchi]);
        
        // TODO
      }
      else if (isA<processorPolyPatch>(bMesh[patchi]))
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
        
        labelListList nepe;
        neighborListEdge(local_EdgePoints, ee, curPointEdges, nepe);
        
        label NN = local_EdgePoints.size();
        
        //if( NN>0){ //ppMeshPoints.size()==0 ||
        
          scalarFieldList weights( NN );
          scalarField sumWeights( NN, 0.0 );

          forAll(nepe, i)
          {
            label  curI = local_EdgePoints[i];
            const point& curP = curPP[curI];
            const labelList& pNeib = nepe[i];

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
            mesh,
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
          //return weights;





          double displ_tol = 1.0;
          int itt = 0;
          vectorField pointNorm( NN, vector::zero );
          scalarList faceToPointSumWeights( NN, 0.0 );

          pointField movedPoints = curPP + pointMotion;

          // *************************************************************************
          // calculate old points normals once
          //vectorField savePointNorm( NN, vector::zero );
          vectorField currentNorms( NN, vector::zero );

          {
            pointField faceCs = faceCentres(curPP, llf);
            vectorField faceNs = faceNormals(curPP, llf);

            scalarList sumWeights( NN, 0.0 );

            forAll(nepe, i)
            {
              label  curI = local_EdgePoints[i];
              const point& curP = curPP[curI];
              const labelList& pNeib = nepe[i];

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

            forAll(nepe, i)
            {
              label  curI = local_EdgePoints[i];
              const point& curP = movedPoints[curI];
              const labelList& pNeib = nepe[i];

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

          // TODO synchronize for periodic
          label fixedEdgePoint = findMin(aa);
          Pout << "fixedEdgePoint:  " << fixedEdgePoint << endl;

          Info<<nl<< "Patch name:  "<< pp.name() << nl
                  << "  number of common points " << NN
                  << "     fixedEdgePoint: " << fixedEdgePoint
                  <<endl;
          // *************************************************************************



          // TODO global var
          double rlxTol = 0.000001;

          while(displ_tol>rlxTol)
          {
            if(itt%q_edge_norm_recalc==0)
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

            forAll(nepe, i)
            {
              label  curI = local_EdgePoints[i];
              point& curP = movedPoints[curI];
              const labelList& pNeib = nepe[i];

              const scalarField& curwv = weights[i];

              forAll(pNeib, ii)
              {
                label ind = pNeib[ii];
                point& neibP = movedPoints[ind];
                vector d2 = neibP - curP;

                scalar mag_d = mag(d2);

                displacement[i] += curwv[ii] * d2;
                tol[i] += mag_d;
                sumWeights[i] += 1;
              }

              // stick to cyclic boundary

              if(itt%q_edge_norm_recalc==0)
              {

                label curIpp = local_pp_EdgePoints[i];
                const vector& curNormPP = ppPointNormals[ curIpp ];

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

            syncTools::syncPointList(mesh, global_EdgePoints, displacement, plusEqOp<vector>(), vector::zero);
            syncTools::syncPointList(mesh, global_EdgePoints, tol, plusEqOp<scalar>(), 0.0);
            syncTools::syncPointList(mesh, global_EdgePoints, sumWeights, plusEqOp<scalar>(), 0.0);

            forAll(pointNorm, i)
            {
              displacement[i] /= sumWeights[i];
            }

            if(itt%q_edge_norm_recalc==0)
            {
              syncTools::syncPointList(mesh, global_EdgePoints, pointNorm, plusEqOp<vector>(), vector::zero);
              syncTools::syncPointList(mesh, global_EdgePoints, faceToPointSumWeights, plusEqOp<scalar>(), 0.0);
              // normalization
              forAll(pointNorm, i)
              {
                pointNorm[i] /= mag( pointNorm[i] );
              }
            }

            if(fixedEdgePoint>=0)
            {
              displacement[fixedEdgePoint] = vector::zero;
            }

            vectorField projectedDisplacement = 
                    transform(I - pointNorm*pointNorm, displacement);

            scalar factor = (itt%q_2edge==0) ? k_2edge : k_1edge;
            projectedDisplacement *= factor;

            if(fixedEdgePoint>=0)
            {
              projectedDisplacement[fixedEdgePoint] = vector::zero;
            }
            
            // stick to cyclic boundary
            //fixPinnedPoints(projectedDisplacement);
            /*
            forAll(projectedDisplacement, ii)
            {
              label ind = local_EdgePoints[ii];
              //label pinnedind = findIndex(pinnedPointsE, ind);
              label pinnedind = findIndex(pinnedPoints, ind);

              if( pinnedind != -1)
              {
                //Info<< "AAAAAAAAAAAAAAAAAAAAA:   " << pinnedind << "   " << ii <<endl;
                
                projectedDisplacement[ii] = 
                        transform
                        (
                          I-pinnedPointsNorm[pinnedind]*pinnedPointsNorm[pinnedind],
                          projectedDisplacement[ii]
                        );
              }
            }
            */
            

            forAll(local_EdgePoints, i)
            {
              label ind = local_EdgePoints[i];
              movedPoints[ind] += projectedDisplacement[i];
            }

            //if(projectedDisplacement.size()>0)
              displ_tol = gAverage( mag(projectedDisplacement/factor)/tol );
            //else
            //  displ_tol = gAverage( 0.0 );

            if(itt%1000==0)
            {
              Info << pp.name() << "  edge rlx iter " << itt
                   << " tolerance: " << displ_tol << endl;
            }

            itt++;
          }

          pointMotion = (movedPoints-curPP);
          
          fixPinnedPoints(pointMotion);
          
          forAll(local_EdgePoints, i)
          {
            label curIpp = local_pp_EdgePoints[i];
            label ind = local_EdgePoints[i];
            Pout<<"ind1111: "<< ind << "  "
                    << curPP[ind]
                    << "  "
                    << pointMotion[ind] 
                    << "  " << ppPointNormals[ curIpp ]
                    << nl;
          }

        //}
          //Pout<< "EEE111:   " << curPP[0] << "   " << pointMotion[0] <<endl;
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

  double displ_tol = 1.0;
  int itt = 0;
  
  // TODO global parameters
  double rlxTol = 0.000001;
  int q_norm_recalc = 1;
  double k_1 = 1.0, k_2 = 2.0;
  int q_2 = 5;
  
  vectorField pointNorm( N, vector::zero );
  scalarList faceToPointSumWeights( N, 0.0 );
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

      forAll(pFaces, j)
      {
        label faceI = pFaces[j];
        point& faceC = faceCs[faceI];

        vector d = faceC - curP;

        scalar mag_d = mag(d);

        vector disp;
        disp.x() = d.x() * pw[j].x();
        disp.y() = d.y() * pw[j].y();
        disp.z() = d.z() * pw[j].z();
        displacement[i] += disp;

        tol[i] += mag_d * mag(pw[j]);

        // this is for normal
        if(itt%q_norm_recalc==0)
        {
          scalar nw = 1.0 / mag_d;
          pointNorm[i] += nw * faceNs[ faceI ];
          faceToPointSumWeights[i] += nw;
        }
      }

      // TODO stick to cyclic boundary
      /*
      label iii = findIndex(pinnedPoints, i);
      if( iii != -1 )
      {
        displacement[i] = 
                transform
                (
                  I-pinnedPointsNorm[iii]*pinnedPointsNorm[iii],
                  displacement[i]
                );
        
      }
       */
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

    if(itt%1000==0)
    {
      Info << this->patch().name() << "  rlx iter " << itt
           << "  tolerance: " << displ_tol << endl;
    }

    itt++;
  }
  Info << this->patch().name() << "  rlx converged in 11 " << itt 
       << " iterations. Tolerance: " << displ_tol<< endl;

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
  
  // in case 0 time does not exist
  if( timeTmp.timeName()!="0" ){
    SeriousErrorIn("meshRelax")
            <<"There is no 0 time directory! "
            "If you run in parallel please check a decomposition."
            <<exit(FatalError);
  }
  
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

        //pinnedPointsLocal.append(local_EdgePoints);
        pinnedPoints.append(local_EdgePoints);
        vectorField locNormals(local_EdgePoints.size(), vector::zero);
        forAll(local_EdgePoints, i)
        {
          locNormals[i] = ppPointNormals[ local_pp_EdgePoints[i] ];
        }
        //pinnedPointsNormLocal.append(locNormals);
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
        //fixedPointsLocal.append(local_EdgePoints);
      }
    }
  }

  //pinnedPoints = pinnedPointsLocal;
  //pinnedPointsNorm = pinnedPointsNormLocal;
  //fixedPoints = fixedPointsLocal;
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
  
  // in case 0 time does not exist
  if( timeTmp.timeName()!="0" ){
    SeriousErrorIn("meshRelax")
            <<"There is no 0 time directory! "
            "If you run in parallel please check a decomposition."
            <<exit(FatalError);
  }
  
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
        if(mag(d[ii])>SMALL){
          w[ii] = 1.0 / mag( d[ii] );
        }
        else{
          w[ii] = 0.0;
        }
      }

      pw[j] = w;
      sumw += pw[j];
    }
    // sum of all distances to face centres
    sumWeights[i] = sumw;
  }

  syncTools::syncPointList( meshTmp, meshPoints, sumWeights, plusEqOp<vector>(), vector::zero);

  forAll(weights, i)
  {
    vectorField& pw = weights[i];
    vector &sw = sumWeights[i];
    forAll(pw, j)
    {
      vector &ppw = pw[j];
      for(int ii=0; ii<3; ii++)
      {
        if( mag(sw[ii])>SMALL )
        {
          ppw[ii] /= sw[ii];
        }
        if( ii==2 )
        {
          ppw[ii] = 0.0;
        }
      }

    }
  }

  surfWeights = weights;
}

/*
Foam::scalarFieldList Foam::
dissolMotionPointPatchVectorField::
calc_edge_weights(const fvMesh& meshTmp, const Patch& patchIO_)
{
  label patchID = this->patch().index();
  const polyMesh& mesh = this->internalField().mesh()();
  const polyBoundaryMesh& bMesh = mesh.boundaryMesh();
  const polyPatch& patch_ = refCast<const polyPatch>(bMesh[patchID]);
  const List<face>& flist = patch_.localFaces();
  const labelListList& plistFaces = patch_.pointFaces();
  const labelList& meshPoints = mesh.boundaryMesh()[patchID].meshPoints();
  

  const pointField& boundaryPoints = patch_.localPoints();
  const pointField& faceCs = patch_.faceCentres();
  
  
  
  
  const labelList& wallsTo  = patch_.meshPoints();
  const labelList& meshPoints = patchIO_.meshPoints();
  labelList local_wall_WallEdges, global_WallEdges;
  commonPoints(wallsTo, meshPoints, local_wall_WallEdges, global_WallEdges);

  const labelListList& pointEdgesLi = patch_.pointEdges();
  const edgeList& ee = patch_.edges();
  labelListList nepe;
  neighborListEdge(local_wall_WallEdges, ee, pointEdgesLi, nepe);

  label NN = local_wall_WallEdges.size();
  scalarFieldList weights( NN );
  scalarField sumWeights( NN, 0.0 );

  forAll(nepe, i){
    label  curI = local_wall_WallEdges[i];
    point& curP = boundaryPoints[curI];
    const labelList& pNeib = nepe[i];

    scalarField& pw = weights[i];
    pw.setSize(pNeib.size());

    scalar sumw = 0.0;
    forAll(pNeib, ii){
      label ind = pNeib[ii];
      point& neibP = boundaryPoints[ind];
      vector d2 = neibP - curP;
      scalar magd2 = mag(d2);

      scalar w = 0.0;
      if( magd2>SMALL) w = 1.0 / magd2;

      pw[ii] = w;
      sumw += pw[ii];
    }
    sumWeights[i] = sumw;
  }

  syncTools::syncPointList(meshTmp, global_WallEdges, sumWeights, plusEqOp<scalar>(), 0.0);

  forAll(weights, i){
    scalarField& pw = weights[i];
    scalar &sw = sumWeights[i];
    forAll(pw, j){
      if( mag(sw)>SMALL ) pw[j] /= sw;
    }
  }

  return weights;
}
*/




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
    if (bMesh[patchi].size() && (patchi != patchID) )
    {
      //Info<< "patch type:  "<< bMesh[patchi].type()<<endl;
      if(
              isA<cyclicPolyPatch>(bMesh[patchi]) 
              ||
              isA<symmetryPolyPatch>(bMesh[patchi])
        )
      {
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
          
          if( mag(pointMotion[pointI])*mag(AA) > SMALL )
          {
            cosa = 
              (pointMotion[pointI] & AA) / (mag(pointMotion[pointI])*mag(AA));
          }
          
          if( mag(cosa)>SMALL )
          {
            scalar L = mag(pointMotion[pointI]) / cosa;
            pointMotion[pointI] = L * AA / mag(AA);
          }
        }
      }
      
      /*
      Info<< "patch type: "<< bMesh[patchi].type()
              << " name: "<< bMesh[patchi].name()
              <<" KKK: "<< pointMotion[0] << endl;
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

  scalar l_T = 0.0;
  const IOdictionary& iod = this->db().objectRegistry::template 
                            lookupObject<IOdictionary>("transportProperties");
  if( !iod.readIfPresent<scalar>("l_T", l_T) )
  {
    SeriousErrorIn("dissolMotionPointPatchVectorField::updateCoeffs()")
            <<"There is no l_T parameter in transportProperties dictionary"
            <<exit(FatalError);
  }

  const volScalarField& C = 
    this->db().objectRegistry::lookupObject<volScalarField>("C");
  
  //const vectorField& motionN = this->patch().pointNormals();

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
  pointMotion = l_T * motionC * motionN;
}

void Foam::dissolMotionPointPatchVectorField::write
(
  Ostream& os
) const
{
    fixedValuePointPatchField<vector>::write(os);
    //timeSeries_.write(os);
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
