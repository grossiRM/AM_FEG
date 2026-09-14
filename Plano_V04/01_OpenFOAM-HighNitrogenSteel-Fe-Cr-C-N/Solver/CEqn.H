
#include "propertiesC.H"
    // construct the flux
    Info << "Executing the flux correction" << endl;

    const volScalarField& alpha = alpha1; 
    
    // Make a crude guess of the phase flux using default interpolation
    surfaceScalarField alphaPhi
    (
        "alphaPhi",
        phi*fvc::interpolate(alpha)
    );
    
    /**********************Get the potential field*********************/
    
    autoPtr<volScalarField> PhiPtr_;

    const volScalarField& p_ = p_rgh;

    wordList PhiPatchFieldTypes(mesh.boundaryMesh().size());
    
    forAll(p_.boundaryField(), patchi)
    {
        PhiPatchFieldTypes[patchi] =
            p_.boundaryField()[patchi].fixesValue()
        ? fixedValueFvPatchField<scalar>::typeName
        : zeroGradientFvPatchField<scalar>::typeName;
    }

    PhiPtr_.set
    (
        new volScalarField
        (
            IOobject
            (
                "PhiY",
                mesh.time().timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(phi.dimensions()/dimLength, Zero),
            PhiPatchFieldTypes
        )
    );
    mesh.setFluxRequired(PhiPtr_->name());
    const volScalarField& Phi = PhiPtr_;
	
    /******************************************************************/
    // Solve for the potential and correct alphaPhi with the resulting flux
    
    if (phi.dimensions() == dimVolume/dimTime)
    {
        while (pimple.correctNonOrthogonal()) 
        {
            fvScalarMatrix PhiEqn
            (
                fvm::laplacian(Phi)
              + fvc::ddt(alpha)
              + fvc::div(alphaPhi)
            );

            PhiEqn.solve(); //this equation solve for Phi!!!

            if (pimple.finalNonOrthogonalIter())
            {
                alphaPhi += PhiEqn.flux();
            }
        }
    }
    else
    {
        SeriousErrorInFunction();
    }
    
    /******************************************************************/
	
	surfaceScalarField alphaDN = fvc::interpolate(alpha)*fvc::interpolate(DN);
	magGradAlpha = mag(fvc::grad(alpha1));
	
    Yloss = magGradAlpha * k_rate * (Ysat - Y) * lf;
    
    fvScalarMatrix alphaYEqn
    (
          fvm::ddt(alpha, Y)
        + fvm::div(alphaPhi, Y)
        - fvm::laplacian(alphaDN, Y)
        == Yloss
    );
    alphaYEqn.relax();
    fvOptions.constrain(alphaYEqn);
    alphaYEqn.solve();
    fvOptions.correct(Y);
  
    Y.max(0.0); 
    
    volScalarField mask = pos(alpha - 0.25);
    YN = mask * Y + (scalar(1) - mask) * scalar(Y0.value());
    
    Info << "min/max(Y) = " << min(Y).value() << ", " << max(Y).value() << endl;
    Info << "min/max(YN) = " << min(YN).value() << ", " << max(YN).value() << endl;






