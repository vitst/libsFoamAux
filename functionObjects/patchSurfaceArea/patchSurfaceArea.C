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

#include "patchSurfaceArea.H"
#include "addToRunTimeSelectionTable.H"
#include "faCFD.H"
#include "OFstreamMod.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(patchSurfaceArea, 0);
    addToRunTimeSelectionTable(functionObject, patchSurfaceArea, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::patchSurfaceArea::patchSurfaceArea
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict)
{
    read(dict);

    current_postPr_dir = "postProcessing/patchSurfaceArea";
    if ( !isDir(current_postPr_dir) ) mkDir(current_postPr_dir);
    current_file_path = current_postPr_dir / "patchSurfaceArea";

    ios_base::openmode mode = ios_base::out|ios_base::trunc;
    OFstreamMod stream(current_file_path, mode);

    stream << "SurfaceArea on " << patchName_ << endl;
    stream << "time   surfaceArea   totalVolume"<<endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::patchSurfaceArea::read
(
    const dictionary& dict
)
{
    fvMeshFunctionObject::read(dict);

    if( !dict.readIfPresent<word>("patchName", patchName_) ){
        SeriousErrorIn("volumeIntegrate::read")
              << "There is no 'patchName' parameter in patchSurfaceArea dictionary"
              << exit(FatalError);
    }

    return true;
}


bool Foam::functionObjects::patchSurfaceArea::execute()
{
    label patchID=mesh_.boundaryMesh().findPatchID(patchName_);
    faMesh fam(mesh_, patchID);

    const polyPatch& curPatch = mesh_.boundaryMesh()[patchID];

    const labelListList& ff = curPatch.faceFaces();

    scalarField fa = fam.S();

    surfArea_ = 0;
    forAll(fa, i)
    {
        surfArea_ += fa[i];
    }
    totalVolume_ = 0;
    forAll(mesh_.V(), cellI)
    {
        totalVolume_ += mesh_.V()[cellI];
    }

    return true;
}


bool Foam::functionObjects::patchSurfaceArea::write()
{
    Log << "    functionObjects::" << type() << " " << name()
        << " writing patch surface area" << endl;

    Log << "S = "<< surfArea_ << "  Total volume: " << totalVolume_ << endl;

    //ios_base::openmode mode = ios_base::out|ios_base::trunc;
    ios_base::openmode mode = ios_base::out|ios_base::app;
    OFstreamMod stream(current_file_path, mode);

    stream << mesh_.time().value() << "   " << surfArea_ 
                                   << "   " << totalVolume_ << endl;

    return true;
}


// ************************************************************************* //
