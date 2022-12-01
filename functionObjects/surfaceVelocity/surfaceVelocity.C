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

#include "surfaceVelocity.H"
//#include "dictionary.H"
#include "addToRunTimeSelectionTable.H"
#include "faCFD.H"
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
    defineTypeNameAndDebug(surfaceVelocity, 0);
    addToRunTimeSelectionTable(functionObject, surfaceVelocity, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::surfaceVelocity::surfaceVelocity
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    //fieldExpression(name, runTime, dict, "phi")
    fvMeshFunctionObject(name, runTime, dict)
    //fieldSet_(mesh_)
{
    read(dict);

    /*
    volScalarField* procFieldPtr
    (
        new volScalarField
        (
            IOobject
            (
                "surfaceVelocity",
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

bool Foam::functionObjects::surfaceVelocity::read
(
    const dictionary& dict
)
{
    fvMeshFunctionObject::read(dict);

    if( !dict.readIfPresent<word>("patchName", patchName_) ){
        SeriousErrorIn("volumeIntegrate::read")
              << "There is no 'patchName' parameter in surfaceVelocity dictionary"
              << exit(FatalError);
    }

    //fieldSet_.read(dict);

    return true;
}

bool Foam::functionObjects::surfaceVelocity::execute()
{
    label patchID=mesh_.boundaryMesh().findPatchID(patchName_);
    faMesh fam(mesh_, patchID);

    const polyPatch& curPatch = mesh_.boundaryMesh()[patchID];
    const labelListList& ff = curPatch.faceFaces();

    //surfVel_ = fam.faceCurvatures();

    faceCenters_ = fam.areaCentres();

    //const volScalarField& field = mesh_.lookupObject<volScalarField>("C");
    //scalarField gradField = -field.boundaryField()[patchID].snGrad();
    const volVectorField& field = mesh_.lookupObject<volVectorField>("grad(C)");
    scalarField gradField = mag(field.boundaryField()[patchID]);

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

    //const IOdictionary& IOd =
    //  mesh_.lookupObject<IOdictionary>("transportProperties");
    //const dictionary& IOd =
    //              mesh_.lookupObject<dictionary>("transportProperties");


    scalar factor =  
        (new dimensionedScalar(
                               "factor_drdt", 
                               dimensionSet(0, 0, 0, 0, 0, 0, 0), 
                               transportProperties
                              ))->value();

    surfVel_ = factor * gradField;


    // averaging
    /*
    scalarField auxCurv(surfVel_.size(), 0);
    scalarField count(surfVel_.size(), 0);
    forAll(surfVel_, i)
    {
        const labelList& ffi = ff[i];
        auxCurv[i] = surfVel_[i];
        count[i] = 1.0 + ffi.size();
        forAll(ffi, j)
        {
            auxCurv[i] += surfVel_[ffi[j]];
        }
    }

    forAll(surfVel_, i)
    {
        surfVel_[i] = auxCurv[i] / count[i];
    }
    */

    return true;
}


bool Foam::functionObjects::surfaceVelocity::write()
{
    Log << "    functionObjects::" << type() << " " << name()
        << " writing surface velocity field and corresponding face centers" << endl;

    Info<< surfVel_.size() << "  " << faceCenters_.size()<<nl;

    word current_postPr_dir = "postProcessing/surfaceVelocity" / mesh_.time().timeName();
    if ( !isDir(current_postPr_dir) ) mkDir(current_postPr_dir);

    fileName current_file_path =
              "postProcessing/surfaceVelocity" / mesh_.time().timeName() / "surfaceVelocity";

    ios_base::openmode mode = ios_base::out|ios_base::trunc;
    //ios_base::openmode mode = ios_base::out|ios_base::app;
    OFstreamMod stream(current_file_path, mode);

    stream << "Velocity on " << patchName_ << " with "
        << faceCenters_.size() << " faces" << endl;
    stream << "face_center        velocity"<<endl;

    forAll(faceCenters_, i)
    {
        stream << faceCenters_[i] << "   " << surfVel_[i] << endl;
    }

    return true;
}


// ************************************************************************* //
