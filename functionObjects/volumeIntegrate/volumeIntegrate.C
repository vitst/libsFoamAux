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
    Info<<"***************************************************************"<<nl;
    Info<<"Dictionary parameters"<<nl;
    Info<< "field names: " << fieldNames_ << nl;
        
    Info<< "min point: " << minPoint_ <<nl;
    Info<< "max point: " << maxPoint_ <<nl;
    Info<<"***************************************************************"<<nl;
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

    //volume_ = 0;
    //volumeIntegral_ = pTraits<Type>::zero;

    read(dict);

    // initialize write file
    word current_postPr_dir = "postProcessing/volumeIntegral";
    if ( !isDir(current_postPr_dir) ) mkDir(current_postPr_dir);

    for (const word& fieldName : fieldNames_)
    {
        fileName current_file_path =
                  "postProcessing/volumeIntegral/volumeIntegral_"+fieldName+".csv";

        ios_base::openmode mode = ios_base::out|ios_base::trunc;
        //ios_base::openmode mode = ios_base::out|ios_base::app;
        OFstreamMod curv_stream(current_file_path, mode);

        curv_stream << "Volume integral of " << fieldName << endl;
        curv_stream << "time,    volume,    filedInt"<<endl;
    }
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
    
    Info<<nl<<nl;
    
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
    //Info<<"Currently there is no write into a file"<<nl;
    Info<<"Currently write is in execute"<<nl;
    
    /*
    //word current_postPr_dir = "postProcessing/volumeIntegral" / mesh_.time().timeName();
    word current_postPr_dir = "postProcessing/volumeIntegral";
    if ( !isDir(current_postPr_dir) ) mkDir(current_postPr_dir);

    fileName current_file_path =
              "postProcessing/volumeIntegral/volumeIntegral";

    //ios_base::openmode mode = ios_base::out|ios_base::trunc;
    ios_base::openmode mode = ios_base::out|ios_base::app;
    OFstreamMod curv_stream(current_file_path, mode);

    curv_stream << mesh_.time().timeName()<< "   " << volume_ << "   " << volumeIntegral_ << endl;
    */

    return true;
}


// ************************************************************************* //
