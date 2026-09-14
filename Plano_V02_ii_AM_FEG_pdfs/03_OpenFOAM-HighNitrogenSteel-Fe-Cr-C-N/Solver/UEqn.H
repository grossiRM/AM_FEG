    volVectorField gradAlpha = fvc::grad(alpha1);
    volVectorField nHat = gradAlpha/(mag(gradAlpha) + mixture.deltaN()); 
    volScalarField factor1 = (scalar(2.0)*rho)/(rho1+rho2);
    volScalarField kappa = -fvc::div(mixture.nHatf());
    Info << "min/max(kappa) = " << min(kappa).value() << ", " << max(kappa).value() << endl;
    Info << "min/max(STForce) = " << min(ST * kappa * nHat  * mag(fvc::grad(alpha1)) * factor1).value() << ", " << max(ST * kappa * nHat  * mag(fvc::grad(alpha1)) * factor1).value() << endl;
    
    MRF.correctBoundaryVelocity(U);

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      //+ MRF.DDt(rho, U)
      //+ turbulence->divDevRhoReff(rho, U) 
      - fvm::laplacian(muf, U)
      - fvm::Sp(A,U) 
      - ST * kappa * nHat  * mag(fvc::grad(alpha1)) * factor1 //surface tension -> interfacial force S_I 1st term
      - dsigmadt * (fvc::grad(T) - nHat * (nHat & fvc::grad(T))) * mag(fvc::grad(alpha1)) * factor1 //Maragoni convection -> S_I 2st term
      - 0.54 * P_0 * exp(( L_v * M * (T - T_v))/( R * T * T_v)) * nHat * mag(fvc::grad(alpha1)) * factor1 * (1-Y) //Recoil pressure -> S_I 3rd term
     ==
        fvOptions(rho, U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
		
		fvc::interpolate(rho* beta * (T - Tref)) * (g & mesh.Sf())
		+
                (
                    //mixture.surfaceTensionForce()*fvc::interpolate(factor1)
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
