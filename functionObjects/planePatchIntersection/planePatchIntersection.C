/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

#include "planePatchIntersection.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "wallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(planePatchIntersection, 0);
    addToRunTimeSelectionTable(functionObject, planePatchIntersection, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::planePatchIntersection::writeFileHeader(Ostream& os) const
{
    // Add headers to output data
    writeHeader(os, "Wall shear stress");
    writeCommented(os, "Time");
    writeTabbed(os, "patch");
    writeTabbed(os, "min");
    writeTabbed(os, "max");
    os << endl;
}


void Foam::functionObjects::planePatchIntersection::calcShearStress
(
    const volSymmTensorField& Reff,
    volVectorField& shearStress
)
{
    shearStress.dimensions().reset(Reff.dimensions());

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchi = iter.key();

        vectorField& ssp = shearStress.boundaryFieldRef()[patchi];
        const vectorField& Sfp = mesh_.Sf().boundaryField()[patchi];
        const scalarField& magSfp = mesh_.magSf().boundaryField()[patchi];
        const symmTensorField& Reffp = Reff.boundaryField()[patchi];

        ssp = (-Sfp/magSfp) & Reffp;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::planePatchIntersection::planePatchIntersection
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    patchSet_()
{
    read(dict);

    writeFileHeader(file());

    volVectorField* planePatchIntersectionPtr
    (
        new volVectorField
        (
            IOobject
            (
                typeName,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector
            (
                "0",
                sqr(dimLength)/sqr(dimTime),
                Zero
            )
        )
    );

    mesh_.objectRegistry::store(planePatchIntersectionPtr);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::planePatchIntersection::~planePatchIntersection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::planePatchIntersection::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    
    if( !dict.readIfPresent<bool>("patchName", patchName_) ){
        SeriousErrorIn("planePatchIntersection::read")
              << "There is no patchName parameter in dissolFoam dictionary"
              << exit(FatalError);
    }
    
    if( !dict.readIfPresent<point>("planePoint", planePoint_) ){
        SeriousErrorIn("planePatchIntersection::read")
              << "There is no planePoint parameter in dissolFoam dictionary"
              << exit(FatalError);
    }
    
    if( !dict.readIfPresent<vector>("planeNormal", planeNormal_) ){
        SeriousErrorIn("planePatchIntersection::read")
              << "There is no planeNormal parameter in dissolFoam dictionary"
              << exit(FatalError);
    }

    if( !dict.readIfPresent<int>("numberOfAngles", numberOfAngles_) ){
        SeriousErrorIn("planePatchIntersection::read")
              << "There is no numberOfAngles parameter in dissolFoam dictionary"
              << exit(FatalError);
    }
    
    return true;
}


bool Foam::functionObjects::planePatchIntersection::execute()
{
    label patchID = mesh_.boundaryMesh().findPatchID(patchName_);

    if(patchID==-1)
    {
        SeriousErrorIn("planePatchIntersection::execute")
                <<"There is no "<< patchName<< exit(FatalError);
    }

    labelHashSet includePatches(1);
    includePatches.insert(patchID);

    triSurface wallTriSurface
    (
      triSurfaceTools::triangulate( mesh_.boundaryMesh(), includePatches )
    );

    // intersection
    const triSurfaceSearch querySurf(wallTriSurface);
    const indexedOctree<treeDataTriSurface>& tree = querySurf.tree();
    
    scalar x0 = planePoint_[0];
    scalar y0 = planePoint_[1];
    scalar a = planeNormal_[0];
    scalar b = planeNormal_[1];
    scalar c = planeNormal_[2];
    
    vector r0(-x0,-y0, (a*x0+b*y0)/c );
    r0 = r0 / mag(r0);
    scalar R = 100;
    
    //scalar tet = twoPi / numberOfAngles;
    
    //vector newr = Foam::cos(tet) + Foam::sin(tet);
    
    //vectot dr =

    for(int ijk=0; ijk<numberOfAngles; ijk++)
    {
      point searchStart(planePoint);

      point searchEnd = searchStart;

      point hitPoint(0.0, 0.0, 0.0);

      pointIndexHit pHit = tree.findLine(searchStart, searchEnd);

      if ( pHit.hit() )
      {
        hitPoint = pHit.hitPoint();
      }

      intersections.append( hitPoint );
    }
    

    return true;
}


bool Foam::functionObjects::planePatchIntersection::write()
{
    const volVectorField& planePatchIntersection =
        obr_.lookupObject<volVectorField>(type());

    Log << type() << " " << name() << " write:" << nl
        << "    writing field " << planePatchIntersection.name() << endl;

    planePatchIntersection.write();

    const fvPatchList& patches = mesh_.boundary();

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchi = iter.key();
        const fvPatch& pp = patches[patchi];

        const vectorField& ssp = planePatchIntersection.boundaryField()[patchi];

        vector minSsp = gMin(ssp);
        vector maxSsp = gMax(ssp);

        if (Pstream::master())
        {
            writeTime(file());

            file()
                << token::TAB << pp.name()
                << token::TAB << minSsp
                << token::TAB << maxSsp
                << endl;
        }

        Log << "    min/max(" << pp.name() << ") = "
            << minSsp << ", " << maxSsp << endl;
    }

    return true;
}


// ************************************************************************* //
