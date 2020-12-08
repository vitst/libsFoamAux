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

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <inttypes.h>
#include <tiffio.h>
#include <string>
#include <cstring>

#include <sstream>      // std::stringstream
//using namespace std;

#include "field2tiff.H"
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
    defineTypeNameAndDebug(field2tiff, 0);
    addToRunTimeSelectionTable(functionObject, field2tiff, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //
void Foam::functionObjects::field2tiff::printInfo() const
{
    Info<<"***************************************************************"<<nl;
    Info<<"Dictionary parameters"<<nl;
    Info<< "field names: " << fieldNames_ << nl;

    Info<< "min point: " << minPoint_ <<nl;
    Info<< "max point: " << maxPoint_ <<nl;
        
    Info<< "unit length for tiff: " << imageUnit_ <<nl;
    Info<<"***************************************************************"<<nl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::field2tiff::field2tiff
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
        Info<< "functionObjects::field2tiff Constructor"<<nl;
    }
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::field2tiff::~field2tiff()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::field2tiff::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    if (debug)
    {
        Info<< "Read dictionary"<<nl;
    }
    
    if( !dict.readIfPresent<wordList>("fields", fieldNames_) ){
        SeriousErrorIn("field2tiff::read")
              << "There is no fields parameter in field2tiff dictionary"
              << exit(FatalError);
    }
    
    if( !dict.readIfPresent<point>("minPoint", minPoint_) ){
        SeriousErrorIn("volumeIntegrate::read")
              << "There is no minPoint parameter in volumeIntegrate dictionary"
              << exit(FatalError);
    }

    if( !dict.readIfPresent<point>("maxPoint", maxPoint_) ){
        SeriousErrorIn("volumeIntegrate::read")
              << "There is no maxPoint parameter in volumeIntegrate dictionary"
              << exit(FatalError);
    }

    if( !dict.readIfPresent<scalar>("imageUnit", imageUnit_) ){
        SeriousErrorIn("field2tiff::read")
              << "There is no imageUnit parameter in field2tiff dictionary"
              << exit(FatalError);
    }
    
    fieldSet_.read(dict);
    
    Info<<"End of dictionary read"<<nl<<nl;
    
    printInfo();
    
    Info<<nl<<nl;
    
    return true;
}


bool Foam::functionObjects::field2tiff::execute()
{
    if(debug)
    {
        Info<<"Integrate fields within the box boundaries"<<nl;
        printInfo();
    }
    
    //for (const word& fieldName : fieldSet_.selection())
    for (const word& fieldName : fieldNames_)
    {
        convertFields<scalar>(fieldName);
        convertFields<vector>(fieldName);
        convertFields<sphericalTensor>(fieldName);
        convertFields<symmTensor>(fieldName);
        convertFields<tensor>(fieldName);
    }
    
    return true;
}

bool Foam::functionObjects::field2tiff::write()
{
    return true;
}


// ************************************************************************* //
