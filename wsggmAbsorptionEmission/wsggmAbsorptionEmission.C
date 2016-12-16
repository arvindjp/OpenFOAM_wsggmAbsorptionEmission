/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "wsggmAbsorptionEmission.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"
#include "basicSpecieMixture.H"
#include "thermoPhysicsTypes.H"
#include "reactingMixture.H"
#include "surfaceFields.H"
#include "symmetryFvPatch.H"
#include "cyclicFvPatch.H"
#include "processorFvPatch.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(wsggmAbsorptionEmission, 0);

        addToRunTimeSelectionTable
        (
            absorptionEmissionModel,
            wsggmAbsorptionEmission,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::wsggmAbsorptionEmission::wsggmAbsorptionEmission
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    absorptionEmissionModel(dict, mesh),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    thermo_(mesh.lookupObject<fluidThermo>(basicThermo::dictName)),
    emissivityCoeffs_(coeffsDict_.lookup("emissivityCoeffs")),
    fittingFactors_(coeffsDict_.lookup("fittingFactors")),
    pathLength_(coeffsDict_.lookup("pathLength")),
    meanBeamPathAutoCalcMode_(coeffsDict_.lookupOrDefault<bool>("meanBeamPathAutoCalcMode", false)),
    sector_(coeffsDict_.lookup("sector"))    
{

    if (!isA<basicSpecieMixture>(thermo_))
    {
        FatalErrorInFunction
            << "Model requires a multi-component thermo package"
            << abort(FatalError);
    }

    label nD = mesh.nGeometricD();
     
    if (nD == 3)
    {

        if (meanBeamPathAutoCalcMode_)
        {
            scalar totVolume = gSum(mesh.V());
            
            scalar totArea = 0.0; 
            forAll(mesh.boundary(),patchI)
            {
                if ( (!isA<processorFvPatch>(mesh.boundary()[patchI]))
                    && (!isA<symmetryFvPatch>(mesh.boundary()[patchI])) 
                    && (!isA<cyclicFvPatch>(mesh.boundary()[patchI])) ) 
                  {
	                    totArea += sum(mesh.magSf().boundaryField()[patchI]);
	              }
            }            
            reduce(totArea, sumOp<scalar>());
    
            pathLength_.value() = 3.6*(360.0/sector_.value())*totVolume/totArea;
            Info << "using the computed pathLength: " <<  pathLength_ << endl;            
        }
        else
        {
            Info << "using the user-provided pathLength: "
                 <<  pathLength_ << endl;    
        }        
        
    }
    else if (nD == 2)
    {
        WarningInFunction
            << "Case is 2D, pathLength is not strictly applicable, \n using the User provided  pathLength: "  
            <<  pathLength_
            << endl;    
    }
    else
    {
        FatalErrorInFunction
            << "Case is not 3D or 2D, pathLength is not applicable"
            << exit(FatalError);
    }
     	
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::wsggmAbsorptionEmission::~wsggmAbsorptionEmission()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::wsggmAbsorptionEmission::aCont(const label bandI) const
{
    const basicSpecieMixture& mixture =
        dynamic_cast<const basicSpecieMixture&>(thermo_);

    const volScalarField& T = thermo_.T();
    const volScalarField& p = thermo_.p();

    label indexCO2 = mixture.species()["CO2"];
    label indexH2O = mixture.species()["H2O"];

    tmp<volScalarField> ta
    (
        new volScalarField
        (
            IOobject
            (
                "aCont" + name(bandI),
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("a", dimless/dimLength, 0.0)
        )
    );
 
    volScalarField& a = ta.ref();
    
    forAll(a, celli)
    {
	    a[celli]=0.0;
	    scalar invWt = 0.0; 
            forAll(mixture.Y(), s)
            {
                invWt += mixture.Y(s)[celli]/mixture.W(s);
            }

	    scalar XkC02 = mixture.Y(indexCO2)[celli]/(mixture.W(indexCO2)*invWt); 
	    scalar XkH20 = mixture.Y(indexH2O)[celli]/(mixture.W(indexH2O)*invWt); 

	    scalar pressurePathLength = (XkC02*paToAtm(p[celli]) + XkH20*paToAtm(p[celli]))*pathLength_.value();
	    scalar limitedDimlessTemperature = min(T[celli], 2400.0);

	    scalar emissivity= 0.0;
	    forAll(emissivityCoeffs_,i)
        	{
		        scalar weightingFactor = 0.0;
		        for(label j=1; j<fittingFactors_.size(); j++)
		        {
		            weightingFactor += fittingFactors_[i][j]*pow(limitedDimlessTemperature, (j-1));
		        }
		        emissivity += weightingFactor * (1 - exp( (-1) *emissivityCoeffs_[i] * pressurePathLength));
        	}
	    emissivity = min(emissivity,0.9999);
	    a[celli] = (-1)* log(1-emissivity) / pathLength_.value();
    } 

    return ta;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::wsggmAbsorptionEmission::eCont(const label bandI) const
{
    return aCont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::wsggmAbsorptionEmission::ECont(const label bandI) const
{
    tmp<volScalarField> tE
    (
        new volScalarField
        (
            IOobject
            (
                "ECont" + name(bandI),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("E", dimMass/dimLength/pow3(dimTime), 0.0)
        )
    );

    return tE;

}


// ************************************************************************* //
