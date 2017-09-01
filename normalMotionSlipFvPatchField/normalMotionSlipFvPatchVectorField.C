/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "normalMotionSlipFvPatchVectorField.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "symmTransformField.H"
#include "addToRunTimeSelectionTable.H"
#include "pointFields.H"

#include "coupledPatchInterpolation.H"


#include "normalMotionSlipPointPatchVectorField.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::normalMotionSlipFvPatchVectorField::
normalMotionSlipFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF)
{}


Foam::normalMotionSlipFvPatchVectorField::
normalMotionSlipFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF)
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


Foam::normalMotionSlipFvPatchVectorField::
normalMotionSlipFvPatchVectorField
(
    const normalMotionSlipFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper)
{}


Foam::normalMotionSlipFvPatchVectorField::
normalMotionSlipFvPatchVectorField
(
    const normalMotionSlipFvPatchVectorField& mwvpvf
)
:
    fixedValueFvPatchVectorField(mwvpvf)
{}


Foam::normalMotionSlipFvPatchVectorField::
normalMotionSlipFvPatchVectorField
(
    const normalMotionSlipFvPatchVectorField& mwvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(mwvpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
  
void Foam::normalMotionSlipFvPatchVectorField::evaluate
(
    const Pstream::commsTypes commsType
)
{
  //Info<<"   normalMotionSlipFvPatchVectorField::evaluate"<<nl;
  
  fixedValueFvPatchField<vector>::evaluate(commsType);
}

void Foam::normalMotionSlipFvPatchVectorField::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
  //Info<<"   normalMotionSlipFvPatchVectorField::initEvaluate"<<nl;
  
  fixedValueFvPatchField<vector>::initEvaluate(commsType);
}


void Foam::normalMotionSlipFvPatchVectorField::updateCoeffs()
{
  //Info<<"   normalMotionSlipFvPatchVectorField::updateCoeffs"<<nl;
  /*
    const fvPatch& p = this->patch();
    const polyPatch& pp = p.patch();
    const fvMesh& mesh = this->internalField().mesh();
    const pointField& points = mesh.points();

    word pfName = this->internalField().name();
    pfName.replace("cell", "point");

    const GeometricField<vector, pointPatchField, pointMesh>& pointMotion =
        this->db().objectRegistry::template
            lookupObject<GeometricField<vector, pointPatchField, pointMesh>>
            (pfName);

    forAll(p, i)
    {
        this->operator[](i) = pp[i].average(points, pointMotion);
    }
    
    
    vectorField pif = this->patchInternalField();
    Info<<"FVFV: pif: "<< max(pif) <<nl;

    fixedValueFvPatchField<vector>::updateCoeffs();
   */
  
    if (updated())
    {
        return;
    }
  
    const fvMesh& mesh = internalField().mesh();

      const fvPatch& p = patch();
        const polyPatch& pp = p.patch();
        const pointField& oldPoints = mesh.oldPoints();

        vectorField oldFc(pp.size());

        forAll(oldFc, i)
        {
            oldFc[i] = pp[i].centre(oldPoints);
        }

        const scalar deltaT = mesh.time().deltaTValue();

        const vectorField Up((pp.faceCentres() - oldFc)/deltaT);  // lg: field of wall velocity

        //const volVectorField& U =
        //    static_cast<const volVectorField&>(internalField()); // lg: field of given vector in U input file
        
        //const volVectorField& U =
        //  this->db().objectRegistry::lookupObject<volVectorField>("U");
        const volScalarField& C =
          this->db().objectRegistry::lookupObject<volScalarField>("C");

        /*
        scalarField phip
        (
            p.patchField<surfaceScalarField, scalar>(fvc::meshPhi(U))
        );
        const vectorField n(p.nf()); // lg: assuming this is unit surface normal field of patch
        const scalarField& magSf = p.magSf();
        tmp<scalarField> Un = phip/(magSf + VSMALL);
        */
        
        const vectorField nHat(this->patch().nf());
        label patchID = this->patch().index();
        
        //Info<<"      iniDisp: "<<iniDisp<<nl;
        //if(iniDisp)
        /*
        {
          const IOdictionary& IOd
                = this->db().lookupObject<IOdictionary>("transportProperties");
          scalar lR =  (new dimensionedScalar(IOd.lookup("lR")))->value();
          scalarField grC = -C.boundaryField()[patchID].snGrad();

          displacement = lR * grC * nHat;
          //iniDisp = false;
        }
         */
        const pointVectorField& pvf = 
          this->db().objectRegistry::lookupObject<pointVectorField>("pointMotionU");
        
        const Foam::normalMotionSlipPointPatchVectorField& pvfBC =
                refCast<const Foam::normalMotionSlipPointPatchVectorField>
                ( (pvf.boundaryField()[patchID]) );
        
        
        vectorField disp2 =
            transform(I - sqr(nHat), this->patchInternalField());
        
        vectorField disp(pp.size());
        /*
        vectorField newNFc(pp.size());
        
        Info<<"      iniDisp: "<<iniDisp<<nl;
        if(iniDisp)
        {
            //const fvMesh& fvmesh_ = refCast<const fvMesh>(mesh);
            coupledPatchInterpolation patchInterpolator
            ( 
              mesh.boundaryMesh()[patchID], mesh
            );
            scalarField motionC=patchInterpolator.faceToPointInterpolate(grC);
            vectorField motionN=patchInterpolator.faceToPointInterpolate(nHat);
            // normalize point normals to 1
            forAll(motionN, ii) motionN[ii]/=mag(motionN[ii]);
            vectorField pointMotion = lR * motionC * motionN;
            
            pointField newP = mesh.boundaryMesh()[patchID].localPoints() + pointMotion;


            forAll(oldFc, i)
            {
                newNFc[i] = pp[i].normal(newP);
            }
          
            disp1 = lR * grC * nHat;
        }
        else
        {
          disp1 = this->patchInternalField();
        }
        vectorField disp2 =
            transform(I - sqr(newNFc), this->patchInternalField());

        iniDisp=false;

        */
        
        disp = pvfBC.getDisp() + disp2;
        
        vectorField::operator=
        ( 
          disp
        );
        
        //vectorField::operator=(nHat*(nHat &  (Up + n*(Un - (n & Up))) )
		//+ transform(I - sqr(nHat), this->patchInternalField()) );
        /*
        vectorField::operator=
        ( 
            this->patchInternalField()
            + 
            transform(I - sqr(nHat), this->patchInternalField()) 
        );
         */
        
        //forAll(p, i)
        //{
        //    this->operator[](i) = disp[i];
        //}
        
        
        /*
        vectorField pif = this->patchInternalField();
        Info<<"FVFV: pif:    "<< max(pif) << nl
            << "      disp:   "<<max(disp) << nl
            << "      disp1:  "<<max(disp1) <<nl
            << "      disp2:  "<<max(disp2) <<nl;
         */
        

    
    fixedValueFvPatchVectorField::updateCoeffs();
}

void Foam::normalMotionSlipFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        normalMotionSlipFvPatchVectorField
    );
}

// ************************************************************************* //
