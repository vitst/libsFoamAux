/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2016 OpenFOAM Foundation
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

#include "apertureAver.H"
#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"
#include "volFields.H"

#include "memInfo.H"
#include "OFstreamMod.H"
#include "interpolation.H"
#include "interpolationCellPoint.H"
#include "triSurfaceSearch.H"
#include "meshSearch.H"

#define NUMBER_OF_COLUMNS 10

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(apertureAver, 0);
    addToRunTimeSelectionTable(functionObject, apertureAver, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //
void Foam::functionObjects::apertureAver::printInfo() const
{
    Info<< "patch names: " << patchNames_ << nl;
    Info<< "field names: " << fieldNames_ << nl;
        
    Info<< "min point: " << minPoint_ <<nl;
    Info<< "max point: " << maxPoint_ <<nl;
    Info<< "flowDirection: " << majDir <<nl;
    Info<< "integrationDirection: " << intDir <<nl;
        
    Info<< "expected number of intersections: " 
            << expectedNumberOfIntersections_ <<nl;
    Info<< "integration points: " << integrationPoints_ <<nl;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::apertureAver::apertureAver
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict)
{
    if (debug)
    {
        Info<< "functionObjects::apertureAver Constructor"<<nl;
    }
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::apertureAver::~apertureAver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::apertureAver::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    /*
    if (fieldName_.empty() || dict.found("field"))
    {
        dict.lookup("field") >> fieldName_;
    }

    if (dict.found("result"))
    {
        dict.lookup("result") >> resultName_;
    }
     */
    
    if (debug)
    {
        Info<< "Read dictionary"<<nl;
    }
    
    if( !dict.readIfPresent<wordList>("fields", fieldNames_) ){
        SeriousErrorIn("apertureAver::read")
              << "There is no fields parameter in apertureAver dictionary"
              << exit(FatalError);
    }
    
    if( !dict.readIfPresent<wordList>("patchNames", patchNames_) ){
        SeriousErrorIn("apertureAver::read")
              << "There is no patchNames parameter in apertureAver dictionary"
              << exit(FatalError);
    }
    
    if( !dict.readIfPresent<point>("minPoint", minPoint_) ){
        SeriousErrorIn("apertureAver::read")
              << "There is no minPoint parameter in apertureAver dictionary"
              << exit(FatalError);
    }
    
    if( !dict.readIfPresent<point>("maxPoint", maxPoint_) ){
        SeriousErrorIn("apertureAver::read")
              << "There is no maxPoint parameter in apertureAver dictionary"
              << exit(FatalError);
    }
    
    if( !dict.readIfPresent<label>("flowDirection", majDir) ){
        SeriousErrorIn("apertureAver::read")
              << "There is no flowDirection parameter "
                "in apertureAver dictionary"
              << exit(FatalError);
    }
    
    if( !dict.readIfPresent<label>("integrationDirection", intDir) ){
        SeriousErrorIn("apertureAver::read")
              << "There is no integrationDirection parameter "
                "in apertureAver dictionary"
              << exit(FatalError);
    }
    
    latDir = 3 - majDir - intDir;
    
    if( !dict.readIfPresent<int>("expectedNumberOfIntersections", 
            expectedNumberOfIntersections_) ){
        SeriousErrorIn("apertureAver::read")
              << "There is no expectedNumberOfIntersections parameter "
                "in apertureAver dictionary"
              << exit(FatalError);
    }
    
    if( !dict.readIfPresent< List<int> >("integrationPoints", 
            integrationPoints_) ){
        SeriousErrorIn("apertureAver::read")
              << "There is no integrationPoints parameter "
                "in apertureAver dictionary"
              << exit(FatalError);
    }
    
    // prepare variables
    
    N_ = integrationPoints_[majDir];
    M_ = integrationPoints_[latDir];
    K_ = integrationPoints_[intDir];
    
    N1_ = integrationPoints_[majDir]+1;
    M1_ = integrationPoints_[latDir]+1;
    K1_ = integrationPoints_[intDir]+1;

    N1M1 = N1_*M1_;
    
    minPosMaj = minPoint_.component(majDir);
    minPosLat = minPoint_.component(latDir);
    minPosInt = minPoint_.component(intDir);

    maxPosMaj = maxPoint_.component(majDir);
    maxPosLat = maxPoint_.component(latDir);
    maxPosInt = maxPoint_.component(intDir);

    // z becomes x; x becomes y
    dx = (maxPosMaj - minPosMaj) / static_cast<scalar>(N_);
    dy = (maxPosLat - minPosLat) / static_cast<scalar>(M_);
    
    thisTimeSize = min(N1M1, maxNumProcPoints);
    curNum = 0;
    curEnd = thisTimeSize;

    totNumLoop = (N1M1-1) / maxNumProcPoints + 1;
    
    Info<<"END READ"<<nl<<nl;
    
    return true;
}


bool Foam::functionObjects::apertureAver::execute()
{
    if(debug)
    {
        Info<<"Calculate averaged fields"<<nl;
    }
    
    Info<<"START EXECUTE"<<nl<<nl;
    
    printInfo();
    
    for(int cI=0; cI<totNumLoop; cI++)
    {
        curNum = cI;
        curBlock = thisTimeSize * curNum;

        sizeAA = thisTimeSize;
        if(cI==totNumLoop-1){
          sizeAA = N1M1 - (totNumLoop-1)*thisTimeSize;
        }

        Info << "Find the points on the surface"<<nl;
        pointsXYonsurface.clear();
        pointsXYonsurface.setSize( expectedNumberOfIntersections_ * sizeAA );

        build_surface_points( mesh_ );
        
        Info << "build_surface_points done"<<nl;

        fileName current_postPr_dir;
        current_postPr_dir = "postProcessing/apertureAver" / time_.timeName();
        if ( !isDir(current_postPr_dir) ) mkDir(current_postPr_dir);

        integrate_and_write(mesh_, time_);
    }
    
    return true;
    /*
    if (foundObject<volVectorField>("U"))
    {
        const volVectorField& uF = lookupObject<volVectorField>("U");
    }
    const volVectorField& uF = lookupObject<volVectorField>("U");
    if (foundObject<volScalarField>("C"))
    {
        const volScalarField& uF = lookupObject<volScalarField>("C");
    }
    
    return true;
    */
}


Foam::scalar Foam::functionObjects::apertureAver::primitive_simpson_integration
(
    scalarField& x,
    scalarField& y  
)
{
    scalar result = 0.0;

    // TODO check len y == len x
    int step = 2;

    scalarField dlt_x( x.size()-1 );

    forAll(dlt_x, i)
    {
      dlt_x = x[i+1]-x[i];
    }

    int len_lab_list = static_cast<int>( std::floor( dlt_x.size() / 2.0 ) );

    labelList slice0( len_lab_list );
    labelList slice1( len_lab_list );
    labelList slice2( len_lab_list );

    int count = 0;
    forAll(slice0, i)
    {
      slice0[i] = count;
      slice1[i] = count+1;
      slice2[i] = count+2;
      count += step;
    }

    scalarField dx0(len_lab_list);
    scalarField dx1(len_lab_list);
    forAll(dx0, j)
    {
      dx0[j] = dlt_x[slice0[j]];
      dx1[j] = dlt_x[slice1[j]];
    }

    scalarField hsum(len_lab_list);     //= h0 + h1
    scalarField hprod(len_lab_list);    //= h0 * h1
    scalarField h0divh1(len_lab_list);  //= h0 / h1


    forAll(dx0, i)
    {
      hsum[i] = dx0[i] + dx1[i];
      hprod[i] = dx0[i] * dx1[i];
      h0divh1[i] = dx0[i] / dx1[i];
    }

    forAll(hsum, i){
      result += hsum[i]/6.0 *
                            (
                              y[slice0[i]] * (2-1.0/h0divh1[i]) +
                              y[slice1[i]] * hsum[i]*hsum[i]/hprod[i] +
                              y[slice2[i]] * (2-h0divh1[i])
                            );
    }

    return result;
}

void Foam::functionObjects::apertureAver::build_surface_points
(
  const fvMesh& mesh
)
{
    // triangulation
    // @TODO in case of parallel calculations see surfaceMeshTriangulate.C
    // @TODO add the error handling for "walls"
    
    bool setMinMax = false;
    scalar maxMaj = 0.0;
    scalar minMaj = 0.0;
    scalar maxLat = 0.0;
    scalar minLat = 0.0;

    labelHashSet includePatches(patchNames_.size());
    forAll(patchNames_, pi)
    {
        label patchID = mesh.boundaryMesh().findPatchID(patchNames_[pi]);

        if(patchID==-1)
        {
          SeriousErrorIn("build_surface_points")
                  <<"There is no "
                  << patchNames_[pi]
                  << exit(FatalError);
        }
        includePatches.insert(patchID);
        
        const pointField &lP = mesh.boundaryMesh()[patchID].localPoints();
        
        if(!setMinMax)
        {
            maxMaj = max( lP.component(majDir) );
            minMaj = min( lP.component(majDir) );
            maxLat = max( lP.component(latDir) );
            minLat = min( lP.component(latDir) );
            setMinMax = true;
        }
        else
        {
            maxMaj = max( maxMaj, max(lP.component(majDir)) );
            minMaj = min( minMaj, min(lP.component(majDir)) );
            maxLat = max( maxLat, max(lP.component(latDir)) );
            minLat = min( minLat, min(lP.component(latDir)) );
        }
    }

    triSurface wallTriSurface
    (
      triSurfaceTools::triangulate( mesh.boundaryMesh(), includePatches )
    );

    // intersection
    const triSurfaceSearch querySurf(wallTriSurface);
    const indexedOctree<treeDataTriSurface>& tree = querySurf.tree();

    scalar curMaj = 0.0;
    scalar curLat = 0.0;


    // N+1 and M+1 because we need points at x=minPosX and x=maxPosX
    for(int ijk=0; ijk<sizeAA; ijk++)
    {
        int totIJK = ijk + curNum * thisTimeSize;

        label i = totIJK / M1_;
        label j = totIJK % M1_;
        curMaj = minPosMaj + i*dx;
        curLat = minPosLat + j*dy;

        label ind = ijk;

        scalar eps=1e-2;
        point searchStart;
        searchStart.component(majDir) = curMaj;
        searchStart.component(latDir) = curLat;
        searchStart.component(intDir) = minPosInt;

        while(searchStart.component(majDir) >= maxMaj ) 
            searchStart.component(majDir) -= eps;
        while(searchStart.component(majDir) <= minMaj ) 
            searchStart.component(majDir) += eps;
        while(searchStart.component(latDir) >= maxLat ) 
            searchStart.component(latDir) -= eps;
        while(searchStart.component(latDir) <= minLat ) 
            searchStart.component(latDir) += eps;

        point searchEnd = searchStart;
        searchEnd.component(intDir) = maxPosInt +
                mag(maxPosInt-minPosInt)/(maxPosInt-minPosInt);

        point hitPoint(0.0, 0.0, 0.0);

        pointIndexHit pHit = tree.findLine(searchStart, searchEnd);

        if ( pHit.hit() )
        {
            hitPoint = pHit.hitPoint();
        }

        label hitWallLabel = expectedNumberOfIntersections_*ind;

        pointsXYonsurface[hitWallLabel] = hitPoint;

        if(expectedNumberOfIntersections_>1)
        {
            // search for second intersection
            searchStart.component(intDir) = hitPoint.component(intDir) + eps;
            pHit = tree.findLine(searchStart, searchEnd);
            label secondHitWallLabel = expectedNumberOfIntersections_ * ind + 1;
            if ( pHit.hit() )
            {
                pointsXYonsurface[secondHitWallLabel] = pHit.hitPoint();
            }
            else
            {
                pointsXYonsurface[hitWallLabel] = point::zero;
                pointsXYonsurface[secondHitWallLabel] = point::zero;
            }
        }

    }
}


void Foam::functionObjects::apertureAver::integrate_and_write
(
    const fvMesh& mesh,
    const Time& runTime
)
{
    /*
    typedef GeometricField<vector, fvPatchField, volMesh> fieldU;
    typedef GeometricField<scalar, fvPatchField, volMesh> fieldC;

    IOobject headerU
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );
    IOobject headerC
    (
        "C",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );
    */


    //if (headerU.headerOk() && headerC.headerOk())
    {
        //fieldU field_u(headerU, mesh);
        //fieldC field_c(headerC, mesh);
        
        //if (foundObject<volScalarField>(fieldName))
        //{
        const volScalarField& cF = lookupObject<volScalarField>("C");
        const volVectorField& uF = lookupObject<volVectorField>("U");
        //}
        

        fileName current_file_path_cup =
                "postProcessing/apertureAver" / runTime.timeName() / "c";
        fileName current_file_path_qx, current_file_path_qy;
        current_file_path_qx =
                "postProcessing/apertureAver" / runTime.timeName() / "qx";
        current_file_path_qy =
                "postProcessing/apertureAver" / runTime.timeName() / "qy";
        fileName current_file_path_h =
                "postProcessing/apertureAver" / runTime.timeName() / "h";
        fileName current_file_path_csurf =
                "postProcessing/apertureAver" / runTime.timeName() / "csurf";

        /*
        static autoPtr<interpolation<vector> > interpolatorU;
        if (!interpolatorU.ptr())
        {
           interpolatorU = 
            interpolation<vector>::New("cellPoint",field_u);
        }

        static autoPtr<interpolation<scalar> > interpolatorC;
        if (!interpolatorC.ptr())
        {
           interpolatorC = 
            interpolation<scalar>::New("cellPoint",field_c);
        }
        */
        
        //const autoPtr<interpolation<scalar> > interpolatorC = 
        //interpolation<scalar> interpolatorC = 
        //        interpolation<scalar>::New("cellPoint", cF);
        
        //const autoPtr<interpolation<vector> > interpolatorU = 
        //interpolation<vector> interpolatorU = 
        //        interpolation<vector>::New("cellPoint", uF);
        
        
        //const interpolationCellPoint<scalar> interpolatorC(cF);
        //const interpolationCellPoint<vector> interpolatorU(uF);

        interpolatorC = interpolation<scalar>::New("cellPoint", cF);
        interpolatorU = interpolation<vector>::New("cellPoint", uF);

        
        
        
        
        

        ios_base::openmode mode =
                (curNum==0) ? 
                  ios_base::out|ios_base::trunc : 
                  ios_base::out|ios_base::app;  
        OFstreamMod mapXYccup( current_file_path_cup, mode );
        OFstreamMod mapXYqx( current_file_path_qx, mode );
        OFstreamMod mapXYqy( current_file_path_qy, mode );
        OFstreamMod mapXYcsurf( current_file_path_csurf, mode);
        OFstreamMod apertureMapFile( current_file_path_h, mode);

        if(curNum==0)
        {
            mapXYccup << N_ << "   " << M_ << "   " << runTime.value() 
                    << "   " << dx << "   " << dy << endl;
            mapXYcsurf << N_ << "   " << M_ << "   " << runTime.value() 
                    << "   " << dx << "   " << dy << endl;
            mapXYqx << N_ << "   " << M_ << "   " << runTime.value()
                    << "   " << dx << "   " << dy << endl;
            mapXYqy << N_ << "   " << M_ << "   " << runTime.value()
                    << "   " << dx << "   " << dy << endl;
            apertureMapFile << N_ << "   " << M_ << "   " << runTime.value() 
                    << "   " << dx << "   " << dy << endl;
        }

        meshSearch searchEngine(mesh);

        int count = 0;
        int iSurf = 0;
        while (iSurf<pointsXYonsurface.size() )
        {
            point first = pointsXYonsurface[iSurf];
            point second = first;

            if(expectedNumberOfIntersections_ == 1)
            {
                second.component(intDir) = minPosInt;
            }
            else
            {
                second = pointsXYonsurface[iSurf+1];
            }

            vector dirc =  first - second;
            scalar dist = mag( dirc );

            scalar csurf = 0.0;
            scalar integratedUx = 0.0;
            scalar integratedUy = 0.0;
            scalar integratedUCx = 0.0;
            scalar integratedUCy = 0.0;
            scalar Ccup = 0.0;
            if( mag(dirc) > SMALL )
            {
                scalarField variable(K1_);
                scalarField Ux(K1_);
                scalarField Uy(K1_);
                scalarField UCx(K1_);
                scalarField UCy(K1_);

                vector dr = 1.0 / static_cast<scalar>(K_) * dirc;

                for(int i=0; i<K1_; i++)
                {
                    point samp_point = second + i*dr;
                    label cellI = searchEngine.findCell( samp_point );
                    if (cellI==-1)
                    {
                        cellI = searchEngine.findNearestCell( samp_point );
                    }
                    label faceI = -1;

                    vector interp_fieldU = 
                            interpolatorU().interpolate(samp_point, cellI, faceI);
                    // velocity = 0 on the surface
                    if(i==0 || i==K_)
                    {
                        interp_fieldU = vector::zero;
                    }

                    variable[i] = samp_point.component(intDir);
                    Ux[i] = interp_fieldU.component(majDir);
                    Uy[i] = interp_fieldU.component(latDir);

                    scalar interp_fieldC = 
                            interpolatorC().interpolate(samp_point, cellI, faceI);
                    UCx[i] = Ux[i] * interp_fieldC;
                    UCy[i] = Uy[i] * interp_fieldC;
                    if(i==K_)
                    {
                        csurf = interp_fieldC;
                    }
                }

                integratedUx = primitive_simpson_integration(variable, Ux);
                integratedUy = primitive_simpson_integration(variable, Uy);
                integratedUCx = primitive_simpson_integration(variable, UCx);
                integratedUCy = primitive_simpson_integration(variable, UCy);
                scalar qSqr = integratedUx*integratedUx+integratedUy*integratedUy;
                if( qSqr > SMALL )
                {
                    Ccup = std::sqrt
                           (
                            (integratedUCx*integratedUCx+integratedUCy*integratedUCy)
                            /
                            qSqr
                           );
                }
            }

            int ind1 = iSurf + curBlock;
            int prcnt  = 100 * (ind1    ) / N1M1 / static_cast<float>(expectedNumberOfIntersections_);
            int prcnt1 = 100 * (ind1 - 1) / N1M1 / static_cast<float>(expectedNumberOfIntersections_);
            if( prcnt%1==0 && prcnt!=prcnt1 )
            {
              Info<<"\r"<< prcnt << "%  "<<flush;
            }

            mapXYccup << Ccup << "  ";
            mapXYcsurf << csurf << "  ";
            mapXYqx << integratedUx << "  ";
            mapXYqy << integratedUy << "  ";
            apertureMapFile << dist << "  ";

            iSurf += expectedNumberOfIntersections_;
            count++;
            if(count >= NUMBER_OF_COLUMNS)
            { 
              mapXYccup <<"\n";
              mapXYcsurf <<"\n";
              mapXYqx <<"\n";
              mapXYqy <<"\n";
              apertureMapFile <<"\n";
              count=0;
            }
        }  
    }
    /*
    else if( !headerU.headerOk() && headerC.headerOk() )
    {
      FatalError<<"There is no U field"<<nl<<nl<<exit(FatalError);
    }
    else if( headerU.headerOk() && !headerC.headerOk() )
    {
      FatalError<<"There is no C field"<<nl<<nl<<exit(FatalError);
    }
    else
    {
      FatalError<<"There is no U and C field"<<nl<<nl<<exit(FatalError);
    }
    */
}



bool Foam::functionObjects::apertureAver::write()
{
    //return writeObject(resultName_);
    Info<<"Write cup averaged field into the file"<<nl;
    return true;
}

/*
bool Foam::functionObjects::fieldExpression::clear()
{
    return clearObject(resultName_);
}
*/


// ************************************************************************* //
