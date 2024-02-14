/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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
#include "OFstreamMod.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::volumeIntegrate::foundObject
(
    const word& name,
    const bool verbose
) const
{
    if (fvMeshFunctionObject::foundObject<Type>(name))
    {
        return true;
    }
    else
    {
        if (debug || verbose)
        {
            Warning
                << "    functionObjects::" << type() << " " << this->name()
                << " cannot find required object " << name << " of type "
                << Type::typeName << endl;
        }

        return false;
    }
}

template<class Type>
void Foam::functionObjects::volumeIntegrate::integrateFields
(
    const word& fieldName
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (obr_.foundObject<fieldType>(fieldName))
    {
        const fieldType& field = lookupObject<fieldType>(fieldName);
        
        Type volIntegral = pTraits<Type>::zero;
        scalar volume = 0;

        // loop over the cells and calculate integral of
        forAll (field, cellI)
        {
            vector pos=mesh_.C()[cellI];
            if( cellInsideTheBox(pos) )
            {
                //volIntegral += field[cellI] / mesh_.V()[cellI];
                volIntegral += field[cellI] * mesh_.V()[cellI];
                //volIntegral += field[cellI];
                volume += mesh_.V()[cellI];
            }
        }
        volIntegral = volIntegral / volume;
        
        Info<<"Volume:  "<<volume<<endl;
        Info<<"Volume integral of "<<fieldName<<":  "<<volIntegral<<endl;

        //volume_ = volume;
        //volumeIntegral_ = volIntegral;
        
        for (const word& fieldName : fieldNames_)
        {
            fileName current_file_path =
                      "postProcessing/volumeIntegral/volumeIntegral_"+fieldName+".csv";
            ios_base::openmode mode = ios_base::out|ios_base::app;
            OFstreamMod curv_stream(current_file_path, mode);
            curv_stream << mesh_.time().timeName()<< ",  " << volume << ",  " << volIntegral << endl;
        }
    }
}



// ************************************************************************* //
