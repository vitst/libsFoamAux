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

#include "volumeIntegrateZ.H"
#include "OFstreamMod.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::volumeIntegrateZ::foundObject
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
void Foam::functionObjects::volumeIntegrateZ::integrateFields
(
    const word& fieldName
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (obr_.foundObject<fieldType>(fieldName))
    {
        const fieldType& field = lookupObject<fieldType>(fieldName);
        
        List<Type> volIntegral(zBins_, pTraits<Type>::zero);
        scalarList volume(zBins_, 0.0);

        scalar hight = maxPoint_.z() - minPoint_.z();
        scalar dz = hight / static_cast<double>(zBins_);

        Info<<"Calculating volume integral"<<endl;
        // loop over the cells and calculate integral of
        forAll (field, cellI)
        {
            vector pos = mesh_.C()[cellI];
            int bin = static_cast<int>(pos.z() / dz);
            volIntegral[bin] += field[cellI] * mesh_.V()[cellI];
            volume[bin] += mesh_.V()[cellI];
        }


        fileName current_time_path = postProc_dir / mesh_.time().timeName();
        if ( !isDir(current_time_path) ) mkDir(current_time_path);
        fileName file_path = current_time_path / fieldName+"_z";
        Info<<"Writing the results into a file:\n  " << file_path << endl;
        //ios_base::openmode mode =
        //        (zBins_==0) ?
        //          ios_base::out|ios_base::trunc :
        //          ios_base::out|ios_base::app;
        ios_base::openmode mode = ios_base::out|ios_base::app;
        OFstreamMod volIntZ( file_path, mode );
        volIntZ << "Zcoord" << "    " << "volume" << "    "
                    << "volumeIntegral" << "    " << "volumeAverage" << nl;

        forAll (volume, binI)
        {
            Type volAver = volIntegral[binI] / volume[binI];

            volIntZ << (binI+0.5)*dz << "   " 
                    << volume[binI] << "   " 
                    << volIntegral[binI] << "   " 
                    << volAver << "   " 
                    << nl;
            
        }
        
    }
}



// ************************************************************************* //
