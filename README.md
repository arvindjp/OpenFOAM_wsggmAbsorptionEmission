# OpenFOAM_wsggmAbsorptionEmission
Weighted Sum of Gray Gases Model (WSGGM) Absorption Emission sub-model for OpenFOAM 4.1

Implemented based on:

 T. F. Smith, Z. F. Shen, and J. N. Friedman. "Evaluation of Coefficients for the Weighted Sum of Gray Gases Model". J. Heat Transfer. 104. 602608. 1982.

 R. Siegel and J. R. Howell. "Thermal Radiation Heat Transfer. Hemisphere Publishing Corporation, Washington DC. 1992.
 
 [FLuent User manual] (http://www.afs.enea.it/project/neptunius/docs/fluent/html/th/node117.htm)

To select WSGGM model use the follwing options in radiationProperties

      absorptionEmissionModel wsggmAbsorptionEmission;

      wsggmAbsorptionEmissionCoeffs
      {

          meanBeamPathAutoCalcMode		true;
          sector sector    [ 0  0  0  0  0  0  0] 360.0; //Degrees
          pathLength pathLength [ 0 1  0  0  0  0  0]  0.1; //m
          emissivityCoeffs 3(  0.4201 6.516 131.9 );
          fittingFactors  
              3        
              (
               4(6.508e-1 -5.551e-4 3.029e-7  -5.353e-11)
               4(-0.2504e-1  6.112e-4  -3.882e-7  6.528e-11)
               4(2.718e-1 -3.118e-4  1.221e-7  -1.612e-11)
              );

      }
