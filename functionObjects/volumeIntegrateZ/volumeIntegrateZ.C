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

#include "volumeIntegrateZ.H"
#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"
#include "volFields.H"

#include "memInfo.H"
#include "OFstreamMod.H"
#include "interpolation.H"
#include "interpolationCellPoint.H"
#include "triSurfaceSearch.H"
#include "meshSearch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(volumeIntegrateZ, 0);
    addToRunTimeSelectionTable(functionObject, volumeIntegrateZ, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //
void Foam::functionObjects::volumeIntegrateZ::printInfo() const
{
    Info<<"***************************************************************"<<nl;
    Info<<"Dictionary parameters"<<nl;
    Info<< "field names: " << fieldNames_ << nl;
        
    Info<< "min point: " << minPoint_ <<nl;
    Info<< "max point: " << maxPoint_ <<nl;

    Info<< "number of bins along z direction: " << zBins_ <<nl;
    Info<<"***************************************************************"<<nl;
}

bool Foam::functionObjects::volumeIntegrateZ::cellInsideTheBox(point& pos) const
{
    bool res = true;
    
    for(int i=0; i<3; i++)
    {
        res = res  &&  ( pos[i]>=minPoint_[i] && pos[i]<=maxPoint_[i] );
    }
    
    return res;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::volumeIntegrateZ::volumeIntegrateZ
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fieldSet_(mesh_)
{
    if (debug)
    {
        Info<< "functionObjects::volumeIntegrateZ Constructor"<<nl;
    }

    read(dict);

    postProc_dir = "postProcessing/volumeIntegrateZ";
    if ( !isDir(postProc_dir) ) mkDir(postProc_dir);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::volumeIntegrateZ::~volumeIntegrateZ()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::volumeIntegrateZ::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    if (debug)
    {
        Info<< "Read dictionary"<<nl;
    }
    
    if( !dict.readIfPresent<wordList>("fields", fieldNames_) ){
        SeriousErrorIn("volumeIntegrateZ::read")
              << "There is no fields parameter in volumeIntegrateZ dictionary"
              << exit(FatalError);
    }
    
    if( !dict.readIfPresent<point>("minPoint", minPoint_) ){
        SeriousErrorIn("volumeIntegrateZ::read")
              << "There is no minPoint parameter in volumeIntegrateZ dictionary"
              << exit(FatalError);
    }
    
    if( !dict.readIfPresent<point>("maxPoint", maxPoint_) ){
        SeriousErrorIn("volumeIntegrateZ::read")
              << "There is no maxPoint parameter in volumeIntegrateZ dictionary"
              << exit(FatalError);
    }

    if( !dict.readIfPresent<int>("zBins", zBins_) ){
        SeriousErrorIn("volumeIntegrateZ::read")
              << "There is no zBins parameter in volumeIntegrateZ dictionary"
              << exit(FatalError);
    }

    //const word format(dict.get<word>("setFormat"));
    //formatterPtr_ = writer<scalar>::New(format);
    //formatterPtr_ = writer<scalar>::New("raw");
    
    fieldSet_.read(dict);
    
    Info<<"END READ"<<nl<<nl;
    
    printInfo();
    
    Info<<nl<<nl;
    
    return true;
}


bool Foam::functionObjects::volumeIntegrateZ::execute()
{
    if(debug)
    {
        Info<<"Integrate fields within the box boundaries"<<nl;
        printInfo();
    }
    
    //for (const word& fieldName : fieldSet_.selection())
    for (const word& fieldName : fieldNames_)
    {
        integrateFields<scalar>(fieldName);
        integrateFields<vector>(fieldName);
        integrateFields<sphericalTensor>(fieldName);
        integrateFields<symmTensor>(fieldName);
        integrateFields<tensor>(fieldName);
    }
    
    return true;
}

bool Foam::functionObjects::volumeIntegrateZ::write()
{
    Info<<"Currently there is no write into a file"<<nl;
    return true;
}


// ************************************************************************* //
