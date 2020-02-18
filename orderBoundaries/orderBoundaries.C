/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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

Application
    createPatch

Group
    grpMeshManipulationUtilities

Description
    Utility to create patches out of selected boundary faces. Faces come either
    from existing patches or from a faceSet.

    More specifically it:
    - creates new patches (from selected boundary faces). Synchronise faces
      on coupled patches.
    - synchronises points on coupled boundaries
    - remove patches with 0 faces in them

\*---------------------------------------------------------------------------*/

#include "cyclicPolyPatch.H"
#include "syncTools.H"
#include "argList.H"
#include "polyMesh.H"
#include "Time.H"
#include "SortableList.H"
#include "OFstream.H"
#include "meshTools.H"
#include "faceSet.H"
#include "IOPtrList.H"
#include "polyTopoChange.H"
#include "polyModifyFace.H"
#include "wordReList.H"
#include "processorMeshes.H"

#include "dynamicFvMesh.H"
#include "normalMotionSlipPointPatchVectorField.H"
#include "codedNormalMotionSlipPointPatchVectorField.H"
#include "pointFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTemplateTypeNameAndDebug(IOPtrList<dictionary>, 0);
}

void changePatchID
(
    const polyMesh& mesh,
    const label faceID,
    const label patchID,
    polyTopoChange& meshMod
)
{
    const label zoneID = mesh.faceZones().whichZone(faceID);

    bool zoneFlip = false;

    if (zoneID >= 0)
    {
        const faceZone& fZone = mesh.faceZones()[zoneID];

        zoneFlip = fZone.flipMap()[fZone.whichFace(faceID)];
    }

    meshMod.setAction
    (
        polyModifyFace
        (
            mesh.faces()[faceID],               // face
            faceID,                             // face ID
            mesh.faceOwner()[faceID],           // owner
            -1,                                 // neighbour
            false,                              // flip flux
            patchID,                            // patch ID
            false,                              // remove from zone
            zoneID,                             // zone ID
            zoneFlip                            // zone flip
        )
    );
}

void changePatchIDforPatch
(
  const polyMesh& mesh,
  const polyPatch& pp,
  label toPatchID,
  polyTopoChange& meshMod
)
{
  forAll(pp, fi)
  {
    changePatchID
    (
      mesh,
      pp.start() + fi,
      toPatchID,
      meshMod
    );
  }
}


int main(int argc, char *argv[])
{
    #include "addOverwriteOption.H"
    #include "addRegionOption.H"
    #include "addDictOption.H"
    Foam::argList::addBoolOption
    (
        "writeObj",
        "write obj files showing the cyclic matching process"
    );
    #include "setRootCase.H"
    #include "createTime.H"
    runTime.functionObjects().off();

    Foam::word meshRegionName = polyMesh::defaultRegion;
    args.optionReadIfPresent("region", meshRegionName);

    const bool overwrite = args.optionFound("overwrite");

    Info<< "Create mesh for time = "
        << runTime.timeName() << nl << endl;

    polyMesh mesh
    (
            IOobject
            (
                dynamicFvMesh::defaultRegion,
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            )
    );

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    const word oldInstance = mesh.pointsInstance();
    
    const word fname("pointMotionU");
    word dname;
    if( exists(runTime.path()/"0"/fname) )
    {
      dname = "0";
    }
    else
    {
      dname = "Zero";
    }
    
    pointVectorField pMU
    (
        IOobject
        (
          "pointMotionU",
          dname,
          mesh,
          IOobject::MUST_READ,
          IOobject::NO_WRITE
        ),
        pointMesh::New(mesh)
    );
    
    Info<<"Patch size: "<<patches.size()<<nl;

    DynamicList<polyPatch*> reorderedPatches(patches.size());
    polyTopoChange meshMod(mesh);
    label startFacei = mesh.nInternalFaces();

    Info<<"  first add moving patches"<<nl;
    Info<<"old patch list:"<<nl;
    forAll(patches, patchi)
    {
      const polyPatch& pp = patches[patchi];
      Info<<" PatchID: "<<patchi<<"  name: "<<pp.name()<< nl;
      
      if (!isA<processorPolyPatch>(pp))
      {
        if
        (
          isA<normalMotionSlipPointPatchVectorField>
                ( pMU.boundaryField()[patchi] ) ||
          isA<codedNormalMotionSlipPointPatchVectorField>
                ( pMU.boundaryField()[patchi] )
        )
        {
          reorderedPatches.append
          (
            pp.clone
            (
              patches,
              patchi,
              pp.size(),
              startFacei
            ).ptr()
          );
          changePatchIDforPatch(mesh, pp, (reorderedPatches.size()-1), meshMod);
          startFacei += pp.size();
        }
      }
    }

    Info<<"  then add the rest"<<nl;
    forAll(patches, patchi)
    {
      const polyPatch& pp = patches[patchi];

      if (!isA<processorPolyPatch>(pp))
      {
        if
        (
          !isA<normalMotionSlipPointPatchVectorField>
                ( pMU.boundaryField()[patchi] ) &&
          !isA<codedNormalMotionSlipPointPatchVectorField>
                ( pMU.boundaryField()[patchi] )
        )
        {
          reorderedPatches.append
          (
            pp.clone
            (
              patches,
              patchi,
              pp.size(),
              startFacei
            ).ptr()
          );
          changePatchIDforPatch(mesh, pp, (reorderedPatches.size()-1), meshMod);
          startFacei += pp.size();
        }
      }
    }

    Info<<"  copy processor patches"<<nl;
    forAll(patches, patchi)
    {
      const polyPatch& pp = patches[patchi];

      if (isA<processorPolyPatch>(pp))
      {
        reorderedPatches.append
        (
          pp.clone
          (
            patches,
            patchi,
            pp.size(),
            startFacei
          ).ptr()
        );
        changePatchIDforPatch(mesh, pp, (reorderedPatches.size()-1), meshMod);
        startFacei += pp.size();
      }
    }
    
    Info<<"  done reordering"<<nl;
    reorderedPatches.shrink();
    mesh.removeBoundary();
    mesh.addPatches(reorderedPatches);
    
    Info<< "Doing topology modification to order faces." << nl << endl;
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, true);
    mesh.movePoints(map().preMotionPoints());
    
    Info<<"New Patch order"<<nl;
    forAll(patches, patchi)
    {
      const polyPatch& pp = patches[patchi];
      Info<<" PatchID: "<<patchi<<"  name: "<<pp.name()<< nl;
    }
    
    // Set the precision of the points data to 10
    IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

    mesh.moving(true);
    
    if (!overwrite)
    {
      runTime++;
    }
    else
    {
      mesh.setInstance(oldInstance);
    }
    
    // Write resulting mesh
    Info<< "Writing repatched mesh to " << mesh.facesInstance() << nl << endl;
    if (!mesh.write())
    {
      FatalErrorInFunction
          << "Failed writing polyMesh."
          << exit(FatalError);
    }
    Info<< "Write done" << endl;
    
    {
      const polyPatchList& patches = mesh.boundaryMesh();

      Info<< "----------------" << nl
          << "Mesh Information" << nl
          << "----------------" << nl
          << "  " << "boundingBox: " << boundBox(mesh.points()) << nl
          << "  " << "nPoints: " << mesh.nPoints() << nl
          << "  " << "nCells: " << mesh.nCells() << nl
          << "  " << "nFaces: " << mesh.nFaces() << nl
          << "  " << "nInternalFaces: " << mesh.nInternalFaces() << nl;

      Info<< "----------------" << nl
          << "Patches" << nl
          << "----------------" << nl;

      forAll(patches, patchi)
      {
        const polyPatch& p = patches[patchi];

        Info<< "  " << "patch " << patchi
            << " (start: " << p.start()
            << " size: " << p.size()
            << ") name: " << p.name()
            << nl;
      }
    }
    
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
