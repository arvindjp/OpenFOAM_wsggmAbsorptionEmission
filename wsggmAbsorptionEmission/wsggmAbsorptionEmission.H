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

Class
    Foam::radiation::wsggmAbsorptionEmission

Description
    Weighted-sum-of-gray-gases-model (WSGGM) for gas phase absorption/emission
    
SourceFiles
    wsggmAbsorptionEmission.C

\*---------------------------------------------------------------------------*/

#ifndef radiationwsggmAbsorptionEmission_H
#define radiationwsggmAbsorptionEmission_H

#include "absorptionEmissionModel.H"
#include "fluidThermo.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                  Class wsggmAbsorptionEmission Declaration
\*---------------------------------------------------------------------------*/

class wsggmAbsorptionEmission
:
    public absorptionEmissionModel
{

    // Private data

        //- Absorption model dictionary
        dictionary coeffsDict_;

        //- SLG thermo package
        const fluidThermo& thermo_;

        //- absorption coefficitent for individual grey gases
        scalarList emissivityCoeffs_;

        //- fitting Factors
        scalarListList fittingFactors_;

        //characteristic path length
        dimensionedScalar pathLength_;
        
        //Flag for meanBeamPath calculation
        bool meanBeamPathAutoCalcMode_;
        
        //Domain sector //[degree]
        dimensionedScalar sector_;       
                
      

public:

    //- Runtime type information
    TypeName("wsggmAbsorptionEmission");


    // Constructors

        //- Construct from components
        wsggmAbsorptionEmission(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~wsggmAbsorptionEmission();


    // Member Functions

        // Access

            // Absorption coefficient

                //- Absorption coefficient for continuous phase
                tmp<volScalarField> aCont(const label bandI = 0) const;


            // Emission coefficient

                //- Emission coefficient for continuous phase
                tmp<volScalarField> eCont(const label bandI = 0) const;


            // Emission contribution

                //- Emission contribution for continuous phase
                tmp<volScalarField> ECont(const label bandI = 0) const;


    // Member Functions

        inline bool isGrey() const
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
