mu1 = FI*mu2;

forAll(mesh.C(), i)
        {
          cps[i] = (470.973*(T[i]+0.5)+7.5e-2*pow((T[i]+0.5),2)+216050+((1.8e+4*x_C.value()+2e+6*pow(x_C.value(),2))/(44*x_C.value()+1200))) - (470.973*(T[i]-0.5)+7.5e-2*pow((T[i]-0.5),2)+216050+((1.8e+4*x_C.value()+2e+6*pow(x_C.value(),2))/(44*x_C.value()+1200)));
          cpl[i] = 825;
          ks[i] = 20.76 + 0.009*T[i] - 3.2627*x_C.value() + (0.0124-2.204e-4*T[i]+1.078e-7*pow(T[i],2) + 7.822e-4*x_Cr.value()-1.741e-7*T[i]*x_Cr.value())*x_Cr.value(); 
          kl[i] = ((35-0.3574*x_Cr.value())-0.0139*Tl.value()) + 0.0139*T[i];
          ST[i] = 1.48206459 - 0.81623e-3*(T[i]-Tl.value()); 
          if (T[i]>Ts.value())
          {
          mu1[i] = 1e-3 * Foam::exp(-2.396 + 7950/T[i]);
          }  
        }
        
        ST= max(ST,ST_limit);
        mu1= max(mu1,mu_limit);        
        kl= min(kl,kl_limit); 

///////////////////////////////////////////////////////////////////////////////
if (nAlphaSubCycles > 1)
{
    dimensionedScalar totalDeltaT = runTime.deltaT();
    surfaceScalarField rhoPhiSum
    (
        IOobject
        (
            "rhoPhiSum",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(rhoPhi.dimensions(), Zero)
    );

    tmp<volScalarField> trSubDeltaT;

    if (LTS)
    {
        trSubDeltaT =
            fv::localEulerDdt::localRSubDeltaT(mesh, nAlphaSubCycles);
    }

    for
    (
        subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
        !(++alphaSubCycle).end();
    )
    {
        #include "alphaEqn.H"
        rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
    }

    rhoPhi = rhoPhiSum;
}
else
{
    #include "alphaEqn.H"
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Info << "min/max(mu) = " << min(mu1).value() << ", " << max(mu1).value() << endl;
mu = alpha1*mu + alpha2*mu2;
muf = fvc::interpolate(mu);

rho = alpha1*rho1 + alpha2*rho2;
rhoCp = alpha1*rho1*(cps*(scalar(1)-lf)+cpl*lf) +  alpha2*rho2*cp2;
Info << "min/max(rhoCp) = " << min(rhoCp).value() << ", " << max(rhoCp).value() << endl;

k =  alpha1*(ks*(scalar(1)-lf)+kl*lf) +  alpha2*k2;
kappaf =  fvc::interpolate(k);
Info << "min/max(k) = " << min(k).value() << ", " << max(k).value() << endl;

beta = alpha1*beta1 +  alpha2*beta2;

A = -C * alpha1 *sqr(scalar(1)-lf)/(pow(lf,scalar(3))+ b);

