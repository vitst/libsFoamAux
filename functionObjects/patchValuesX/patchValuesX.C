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

#include "patchValuesX.H"
#include "addToRunTimeSelectionTable.H"
//#include "faCFD.H"
#include "OFstreamMod.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(patchValuesX, 0);
    addToRunTimeSelectionTable(functionObject, patchValuesX, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::patchValuesX::patchValuesX
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fieldSet_(mesh_)
{
    read(dict);
    //postProc_dir = "postProcessing/volumeIntegrateZ";
    //if ( !isDir(postProc_dir) ) mkDir(postProc_dir);

    /*
    volScalarField* procFieldPtr
    (
        new volScalarField
        (
            IOobject
            (
                "patchValuesX",
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

bool Foam::functionObjects::patchValuesX::read
(
    const dictionary& dict
)
{
    fvMeshFunctionObject::read(dict);

    if( !dict.readIfPresent<wordList>("fields", fieldNames_) ){
        SeriousErrorIn("patchValuesX::read")
              << "There is no fields parameter in patchValuesX dictionary"
              << exit(FatalError);
    }

    if( !dict.readIfPresent<word>("patchName", patchName_) ){
        SeriousErrorIn("patchValuesX::read")
              << "There is no 'patchName' parameter in patchValuesX dictionary"
              << exit(FatalError);
    }
    fieldSet_.read(dict);
    Info<<"END READ"<<nl<<nl;
    Info<<nl<<nl;

    return true;
}


bool Foam::functionObjects::patchValuesX::execute()
{
    if(debug)
    {
        Info<<"Get boundary of a field"<<nl;
    }
    
    //for (const word& fieldName : fieldSet_.selection())
    for (const word& fieldName : fieldNames_)
    {
        getBoundaryValues<scalar>(fieldName);
        getBoundaryValues<vector>(fieldName);
        getBoundaryValues<sphericalTensor>(fieldName);
        getBoundaryValues<symmTensor>(fieldName);
        getBoundaryValues<tensor>(fieldName);
    }
    
    return true;
}


bool Foam::functionObjects::patchValuesX::write()
{
    Log << "    functionObjects::" << type() << " " << name()
        << " writing curvature field and corresponding cell centers" << endl;

    Info<< values_.size() << "  " << faceCenters_.size()<<nl;

    word current_postPr_dir = "postProcessing/patchValuesX" / mesh_.time().timeName();
    if ( !isDir(current_postPr_dir) ) mkDir(current_postPr_dir);

    fileName current_file_path =
              "postProcessing/patchValuesX" / mesh_.time().timeName() / "patchValuesX";

    ios_base::openmode mode = ios_base::out|ios_base::trunc;
    //ios_base::openmode mode = ios_base::out|ios_base::app;
    OFstreamMod curv_stream(current_file_path, mode);

    curv_stream << "ValuesX on " << patchName_ << " with "
        << faceCenters_.size() << " faces" << endl;
    curv_stream << "x_coord        value"<<endl;

    forAll(faceCenters_, i)
    {
        curv_stream << faceCenters_[i].x() << "   " << values_[i] << endl;
    }

    //const volScalarField& distance =
    //    mesh_.lookupObject<volScalarField>("patchValuesX");
    //distance.write();

    return true;
}


// ************************************************************************* //
