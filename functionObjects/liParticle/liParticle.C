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

\*---------------------------------------------------------------------------*/

#include "liParticle.H"
//#include "dictionary.H"
#include "addToRunTimeSelectionTable.H"
//#include "faCFD.H"
//#include "fvCFD.H"
#include "OFstreamMod.H"
//#include "transportModel.H"

//#include "turbulenceModel.H"
//#include "surfaceInterpolate.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(liParticle, 0);
    addToRunTimeSelectionTable(functionObject, liParticle, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::liParticle::liParticle
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    //fieldExpression(name, runTime, dict, "phi")
    fvMeshFunctionObject(name, runTime, dict),
    fieldSet_(mesh_)
{
    curTimeDim = 0;
    minFieldX = 0;
    maxFieldX = 0;
    minFieldY = 0; 
    maxFieldY = 0;
    minX = 0;
    maxX = 0;
    minY = 0; 
    maxY = 0;
        
    read(dict);

    // initialize write file
    word current_postPr_dir = "postProcessing/liParticle";
    if ( !isDir(current_postPr_dir) ) mkDir(current_postPr_dir);
    
    fileName current_file_path = "postProcessing/liParticle/values.csv";

    ios_base::openmode mode = ios_base::out|ios_base::trunc;
    //ios_base::openmode mode = ios_base::out|ios_base::app;
    OFstreamMod curv_stream(current_file_path, mode);

    curv_stream << "time,time_s,minX,maxX,minY,maxY,dX,dY,vfront,vtail,vlow,vhigh"<<endl;

    /*
    volScalarField* procFieldPtr
    (
        new volScalarField
        (
            IOobject
            (
                "liParticle",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimLength, Zero)
        )
    );

    mesh_.objectRegistry::store(procFieldPtr);
    */
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::liParticle::read
(
    const dictionary& dict
)
{
    fvMeshFunctionObject::read(dict);

    if (debug)
    {
        Info<< "Read dictionary"<<nl;
    }

    if( !dict.readIfPresent<word>("patchName", patchName_) ){
        SeriousErrorIn("liParticle::read")
              << "There is no 'patchName' parameter in liParticle dictionary"
              << exit(FatalError);
    }

    if( !dict.readIfPresent<wordList>("fields", fieldNames_) ){
        SeriousErrorIn("liParticle::read")
              << "There is no fields parameter in liParticle dictionary"
              << exit(FatalError);
    }

    if( !dict.readIfPresent<scalar>("nCA", nCA_) ){
        SeriousErrorIn("liParticle::read")
              << "There is no nCA parameter in liParticle dictionary"
              << exit(FatalError);
    }
    if( !dict.readIfPresent<scalar>("nCB", nCB_) ){
        SeriousErrorIn("liParticle::read")
              << "There is no nCB parameter in liParticle dictionary"
              << exit(FatalError);
    }

    fieldSet_.read(dict);

    return true;
}

bool Foam::functionObjects::liParticle::execute()
{
    label patchID=mesh_.boundaryMesh().findPatchID(patchName_);
    if(patchID==-1)
    {
        SeriousErrorIn("liParticle::execute")
              << "Patch name " << patchName_ << " does not exist. Check patch names."
              << exit(FatalError);
    }

    const polyPatch& curPatch = mesh_.boundaryMesh()[patchID];
    //const labelListList& ff = curPatch.faceFaces();

    faceCenters_ = curPatch.faceCentres();

    const volScalarField& fieldA = mesh_.lookupObject<volScalarField>("CA");
    //scalarField gradField = -field.boundaryField()[patchID].snGrad();
    //const volVectorField& field = mesh_.lookupObject<volVectorField>("grad(C)");
    //scalarField gradField = mag(field.boundaryField()[patchID]);
    const volScalarField& fieldB = mesh_.lookupObject<volScalarField>("CB");

    const scalarField& CAb = fieldA.boundaryField()[patchID];
    const scalarField& CBb = fieldB.boundaryField()[patchID];


    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );


    // scaling of displacement (effectively scaled timestep)
    /*
    scalar aux_dxf =  
        (new dimensionedScalar( "aux_dxf", 
                                 dimensionSet(0, 0, 0, 0, 0, 0, 0), 
                                 transportProperties))->value();
    */
    scalar aux_dxf = readScalar(transportProperties.lookup("aux_dxf"));

    scalar D_ =  
        (new dimensionedScalar( "DA", 
                                 dimensionSet(0, 2, -1, 0, 0, 0, 0), 
                                 transportProperties))->value();
    scalar l_R_ =  
        (new dimensionedScalar( "l_RA", 
                                 dimensionSet(0, 0, 0, 0, 0, 0, 0), 
                                 transportProperties))->value();
    scalar theta =  
        (new dimensionedScalar( "theta", 
                                 dimensionSet(0, 0, 0, 0, 0, 0, 0), 
                                 transportProperties))->value();

    // dimensionless time
    //const scalar dt = mesh_.time().deltaTValue();
    scalar curTime = mesh_.time().value();

    scalar k = D_ * u0 * c0 / l_R_;
    // surface velocity unit in cm/s
    scalar v0 = nu_m * k;

    // time unit in s
    scalar td = h0 / v0;

    // time in s
    curTimeDim = curTime * aux_dxf * td;

    // dimensionless velosity is displacement / dt
    //scalarField surfVel_ = aux_dxf * dt * (theta * Foam::pow(CAb, nCA_) * Foam::pow(CBb, nCB_) - 1);
    scalarField surfVel_ = (theta * Foam::pow(CAb, nCA_) * Foam::pow(CBb, nCB_) - 1);
    
    // find min max in X and Y
    scalarField posX = faceCenters_.component(vector::X);
    scalarField posY = faceCenters_.component(vector::Y);
    labelPair minMaxIdsX = findMinMax(posX);
    labelPair minMaxIdsY = findMinMax(posY);

    label minIdX = minMaxIdsX.first();
    if (minIdX != -1)
    {
        minX = posX[minIdX];
        minFieldX = surfVel_[minIdX];
    }

    label maxIdX = minMaxIdsX.second();
    if (maxIdX != -1)
    {
        maxX = posX[maxIdX];
        maxFieldX = surfVel_[maxIdX];
    }

    label minIdY = minMaxIdsY.first();
    if (minIdY != -1)
    {
        minY = posY[minIdY];
        minFieldY = surfVel_[minIdY];
    }

    label maxIdY = minMaxIdsY.second();
    if (maxIdY != -1)
    {
        maxY = posY[maxIdY];
        maxFieldY = surfVel_[maxIdY];
    }

    return true;
}


bool Foam::functionObjects::liParticle::write()
{
    Log << "    functionObjects::" << type() << " " << name()
        << " writing post-proc values for Li3PO4 particle" << endl;
    //Info<< surfVel_.size() << "  " << faceCenters_.size()<<nl;

    fileName current_file_path = "postProcessing/liParticle/values.csv";

    //ios_base::openmode mode = ios_base::out|ios_base::trunc;
    ios_base::openmode mode = ios_base::out|ios_base::app;

    OFstreamMod stream(current_file_path, mode);

    stream << mesh_.time().timeName() << "," 
        << curTimeDim << ","
        << minX << ","
        << maxX << ","
        << minY << ","
        << maxY << ","
        << (maxX-minX) << ","
        << (maxY-minY) << ","
        << minFieldX << ","
        << maxFieldX << ","
        << minFieldY << ","
        << maxFieldY
        << endl;

    return true;
}


// ************************************************************************* //
