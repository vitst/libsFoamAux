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

#include "patchCurvature.H"
#include "addToRunTimeSelectionTable.H"
#include "faCFD.H"
#include "OFstreamMod.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(patchCurvature, 0);
    addToRunTimeSelectionTable(functionObject, patchCurvature, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::patchCurvature::patchCurvature
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict)
{
    read(dict);

    /*
    volScalarField* procFieldPtr
    (
        new volScalarField
        (
            IOobject
            (
                "patchCurvature",
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

bool Foam::functionObjects::patchCurvature::read
(
    const dictionary& dict
)
{
    fvMeshFunctionObject::read(dict);

    if( !dict.readIfPresent<word>("patchName", patchName_) ){
        SeriousErrorIn("volumeIntegrate::read")
              << "There is no 'patchName' parameter in patchCurvature dictionary"
              << exit(FatalError);
    }

    return true;
}


bool Foam::functionObjects::patchCurvature::execute()
{
    label patchID=mesh_.boundaryMesh().findPatchID(patchName_);
    faMesh fam(mesh_, patchID);

    const polyPatch& curPatch = mesh_.boundaryMesh()[patchID];

    const labelListList& ff = curPatch.faceFaces();

    curvature_ = fam.faceCurvatures();

    faceCenters_ = fam.areaCentres();

    scalarField fa = fam.S();

    Info<<"Min face area: "<<Foam::min(fa)<<nl;
    Info<<"Max face curv: "<<Foam::max(curvature_)<<nl;

    // averaging
    /*
    scalarField auxCurv(curvature_.size(), 0);
    scalarField count(curvature_.size(), 0);
    forAll(curvature_, i)
    {
        const labelList& ffi = ff[i];
        auxCurv[i] = curvature_[i];
        count[i] = 1.0 + ffi.size();
        forAll(ffi, j)
        {
            auxCurv[i] += curvature_[ffi[j]];
        }
    }

    forAll(curvature_, i)
    {
        curvature_[i] = auxCurv[i] / count[i];
    }
    */

    return true;
}


bool Foam::functionObjects::patchCurvature::write()
{
    Log << "    functionObjects::" << type() << " " << name()
        << " writing curvature field and corresponding cell centers" << endl;

    Info<< curvature_.size() << "  " << faceCenters_.size()<<nl;

    word current_postPr_dir = "postProcessing/patchCurvature" / mesh_.time().timeName();
    if ( !isDir(current_postPr_dir) ) mkDir(current_postPr_dir);

    fileName current_file_path =
              "postProcessing/patchCurvature" / mesh_.time().timeName() / "patchCurvature";

    ios_base::openmode mode = ios_base::out|ios_base::trunc;
    //ios_base::openmode mode = ios_base::out|ios_base::app;
    OFstreamMod curv_stream(current_file_path, mode);

    curv_stream << "Curvature on " << patchName_ << " with "
        << faceCenters_.size() << " faces" << endl;
    curv_stream << "face_center        curvature"<<endl;

    forAll(faceCenters_, i)
    {
        curv_stream << faceCenters_[i] << "   " << curvature_[i] << endl;
    }

    //const volScalarField& distance =
    //    mesh_.lookupObject<volScalarField>("patchCurvature");
    //distance.write();

    return true;
}


// ************************************************************************* //
