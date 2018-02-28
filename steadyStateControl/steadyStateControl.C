/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "steadyStateControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(steadyStateControl, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::steadyStateControl::read()
{
    solutionControl::read(true);
}


bool Foam::steadyStateControl::criteriaSatisfied()
{
    if (residualControl_.empty())
    {
        return false;
    }

    Time& time = const_cast<Time&>(mesh_.time());
    if (debug)
    {
        Info<<"Time: "<< time.timeName() <<"    Iteration: "<<iter_counter<<nl;
    }

    bool achieved = true;
    bool checked = false;    // safety that some checks were indeed performed

    const dictionary& solverDict = mesh_.solverPerformanceDict();
    forAllConstIter(dictionary, solverDict, iter)
    {
        const word& variableName = iter().keyword();
        const label fieldi = applyToField(variableName);
        if (fieldi != -1)
        {
            scalar lastResidual = 0;
            const scalar residual =
                this->maxResidual(variableName, iter().stream(), lastResidual);

            checked = true;

            bool absCheck = residual < residualControl_[fieldi].absTol;
            achieved = achieved && absCheck;

            if (debug)
            {
                Info<< algorithmName_ << " solution statistics:";

                Info<< "    " << variableName << ": tolerance = " << residual
                    << " (" << residualControl_[fieldi].absTol << ")"
                    << endl;
            }
        }
    }
    
    // clear dictionary
    dictionary& dict = const_cast<dictionary&>(mesh_.solverPerformanceDict());
    dict.clear();

    initialised_ = false;


    return checked && achieved;
}

template<class Type>
void Foam::steadyStateControl::maxTypeResidual
(
    const word& fieldName,
    ITstream& data,
    scalar& firstRes,
    scalar& lastRes
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;
    if (mesh_.foundObject<fieldType>(fieldName))
    {
	    fieldType& field = const_cast<fieldType&>
        (
		    mesh_.lookupObject<fieldType>( fieldName )
        );

        const List<SolverPerformance<Type> > sp(data);

        int sz = field.size();
        reduce(sz, sumOp<int>());
        Type nF = gSumCmptProd( field, field ) / static_cast<double>(sz);

        scalar norm = mag(nF);

        for (direction cmpt=0; cmpt < nF.size(); cmpt++)
        {
            nF[cmpt] = sqrt( nF[cmpt] );
        }
        nF = nF/norm;

        firstRes = cmptMax( cmptMultiply(sp.first().initialResidual(), nF) );
        lastRes  = cmptMax( cmptMultiply(sp.last().initialResidual(),  nF) );

    }
}


Foam::scalar Foam::steadyStateControl::maxResidual
(
    const word& fieldName,
    ITstream& data,
    scalar& lastRes
) const
{
    scalar firstRes = 0;

    if(mesh_.foundObject<volScalarField>(fieldName))
    {
        const List< SolverPerformance<scalar> > sp(data);
        firstRes = sp.first().initialResidual();
        lastRes = sp.last().initialResidual();
    }
    else
    {
        maxTypeResidual<vector>(fieldName, data, firstRes, lastRes);
        maxTypeResidual<sphericalTensor>(fieldName, data, firstRes, lastRes);
        maxTypeResidual<symmTensor>(fieldName, data, firstRes, lastRes);
        maxTypeResidual<tensor>(fieldName, data, firstRes, lastRes);
    }

    return firstRes;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::steadyStateControl::steadyStateControl(fvMesh& mesh)
:
    solutionControl(mesh, "SIMPLE"),
    initialised_(false),
    iter_counter(0)
{
    read();

    Info<<nl;

    if (residualControl_.empty())
    {
        Info<< algorithmName_ << ": no convergence criteria found. Exiting"<< endl;
        exit(FatalIOError);
    }
    else
    {
        Info<< algorithmName_ << ": convergence criteria" << nl;
        forAll(residualControl_, i)
        {
            Info<< "    field " << residualControl_[i].name << token::TAB
                << " tolerance " << residualControl_[i].absTol
                << nl;
        }
        Info<< endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::steadyStateControl::loop()
{
    bool run = true;
    read();

    if (initialised_)
    {
        iter_counter++;
        if (criteriaSatisfied())
        {
            Info<<nl<<algorithmName_<<" solution converged in "<<iter_counter<<" iterations"<<nl<<nl;
            iter_counter = 0;
            run = false;
        }
        else
        {
            storePrevIterFields();
        }
    }
    else
    {
        initialised_ = true;
        storePrevIterFields();
    }

    return run;
}


// ************************************************************************* //
