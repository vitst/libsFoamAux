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

#include "volumeIntegrate.H"
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
    defineTypeNameAndDebug(volumeIntegrate, 0);
    addToRunTimeSelectionTable(functionObject, volumeIntegrate, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //
void Foam::functionObjects::volumeIntegrate::printInfo() const
{
    Info<< "field names: " << fieldNames_ << nl;
        
    Info<< "min point: " << minPoint_ <<nl;
    Info<< "max point: " << maxPoint_ <<nl;
}

bool Foam::functionObjects::volumeIntegrate::cellInsideTheBox(point& pos) const
{
    bool res = true;
    
    for(int i=0; i<3; i++)
    {
        res = res  &&  ( pos[i]>=minPoint_[i] && pos[i]<=maxPoint_[i] );
    }
    
    return res;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::volumeIntegrate::volumeIntegrate
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
        Info<< "functionObjects::volumeIntegrate Constructor"<<nl;
    }
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::volumeIntegrate::~volumeIntegrate()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::volumeIntegrate::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    if (debug)
    {
        Info<< "Read dictionary"<<nl;
    }
    
    if( !dict.readIfPresent<wordList>("fields", fieldNames_) ){
        SeriousErrorIn("volumeIntegrate::read")
              << "There is no fields parameter in volumeIntegrate dictionary"
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
    
    fieldSet_.read(dict);
    
    Info<<"END READ"<<nl<<nl;
    
    printInfo();
    
    return true;
}


bool Foam::functionObjects::volumeIntegrate::execute()
{
    if(debug)
    {
        Info<<"Integrate fields within the box boundaries"<<nl;
        printInfo();
    }
    
    //for (const word& fieldName : fieldSet_.selection())
    for (const word& fieldName : fieldNames_)
    {
        Info<<"working filed: "<<fieldName<<nl;
        integrateFields<scalar>(fieldName);
        integrateFields<vector>(fieldName);
        integrateFields<sphericalTensor>(fieldName);
        integrateFields<symmTensor>(fieldName);
        integrateFields<tensor>(fieldName);
    }
    
    return true;
}

bool Foam::functionObjects::volumeIntegrate::write()
{
    Info<<"At the current moment there is no write into the file"<<nl;
    return true;
}


// ************************************************************************* //
