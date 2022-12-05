/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
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

#include "curv.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"

#include "faCFD.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(curv, 0);
    addToRunTimeSelectionTable(functionObject, curv, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::curv::calc()
{
    faMesh fam(mesh_);

    /*
    areaScalarField curv
    (
        IOobject
        (
            "curvature",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fam,
        dimensionedScalar(
            "curvature", dimensionSet(0, -1, 0, 0, 0, 0, 0), 1.0
        )
    );
    */

    //curv = fam.faceCurvatures();
    const scalarField& K = fam.faceCurvatures().internalField();
    //curv = fam.faceCurvatures();

    //Info<<"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"<<endl;
    //Info<<K.size()<<endl;

    if( K.size() > 0)
    {
        int numPatches = mesh_.boundaryMesh().size();
        wordList boundaryTypes(numPatches, "fixedValue"); // default bc
      
        volScalarField curv
        (
            IOobject
            (
                "curvature",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar(
                "curvature", dimensionSet(0, -1, 0, 0, 0, 0, 0), 0.0
            ),
            boundaryTypes
        );
    
        const label& patchID = mesh_.boundaryMesh().findPatchID("reactive_surface");
        curv.boundaryFieldRef()[patchID] = K;
        //curv.correctBoundaryConditions();
        const labelList& fc = mesh_.boundaryMesh()[patchID].faceCells();
        forAll(fc, i)
        {
            const label& cI = fc[i];
            curv[cI] = K[i];
        }
    
        //Info<<fc.size()<<nl;
        //Info<<curv.internalField()<<endl;
        //
        const tmp<volScalarField> cc(curv);
    
        return store
        (
            resultName_,
            cc
        );
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::curv::curv
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, "curv")
{
    setResultName(typeName, "curv");
}


// ************************************************************************* //
