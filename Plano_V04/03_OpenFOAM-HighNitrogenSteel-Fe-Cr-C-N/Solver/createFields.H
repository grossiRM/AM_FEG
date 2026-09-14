#include "createRDeltaT.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dimensionedScalar x_C("x_C", dimensionSet(0, 0, 0, 0, 0, 0, 0),transportProperties);
dimensionedScalar x_Cr("x_Cr", dimensionSet(0, 0, 0, 0, 0, 0, 0),transportProperties);
dimensionedScalar ST_ini("ST_ini", dimensionSet(1, 0, -2, 0, 0, 0, 0),transportProperties);

dimensionedScalar cp2("cp2", dimensionSet(0, 2, -2, -1, 0, 0, 0), transportProperties);

dimensionedScalar k2("k2", dimensionSet(1, 1, -3, -1, 0, 0, 0), transportProperties);

dimensionedScalar mu2("mu2", dimensionSet(1, -1, -1, 0, 0, 0, 0),transportProperties);
//dimensionedScalar mu_limit("mu_limit", dimensionSet(1, -1, -1, 0, 0, 0, 0),transportProperties);

dimensionedScalar beta1("beta1", dimensionSet(0, 0, 0, -1, 0, 0, 0), transportProperties);
dimensionedScalar beta2("beta2", dimensionSet(0, 0, 0, -1, 0, 0, 0), transportProperties);

dimensionedScalar Ts("Ts", dimensionSet(0, 0, 0, 1, 0, 0, 0), transportProperties);
dimensionedScalar Tl("Tl", dimensionSet(0, 0, 0, 1, 0, 0, 0), transportProperties);
dimensionedScalar Tref("Tref", dimensionSet(0, 0, 0, 1, 0, 0, 0), transportProperties);
dimensionedScalar C("C", dimensionSet(1, -3, -1, 0, 0, 0, 0), transportProperties);
dimensionedScalar b("b", dimensionSet(0, 0, 0, 0, 0, 0, 0), transportProperties);
dimensionedScalar dsigmadt("dsigmadt", dimensionSet(1, 0, -2, -1, 0, 0, 0), transportProperties);
dimensionedScalar E("E", dimensionSet(0, 0, 0, 0, 0, 0, 0), transportProperties);
dimensionedScalar h1("h1", dimensionSet(1, 0, -3, -1, 0, 0, 0), transportProperties);
dimensionedScalar sig("sig", dimensionSet(1, 0, -3, -4, 0, 0, 0), transportProperties);
dimensionedScalar P_0("P_0", dimensionSet(1, -1, -2, 0, 0, 0, 0), transportProperties);
dimensionedScalar R("R", dimensionSet(1, 2, -2, -1, -1, 0, 0), transportProperties);

dimensionedScalar Y0("Y0", dimensionSet(0, 0, 0, 0, 0, 0, 0), transportProperties);
dimensionedScalar k_init("k_init", dimensionSet(0, 1, -1, 0, 0, 0, 0), transportProperties);
dimensionedScalar Sh("Sh", dimensionSet(0, 0, 0, 0, 0, 0, 0), transportProperties);
dimensionedScalar Sh_L("Sh_L", dimensionSet(0, 1, 0, 0, 0, 0, 0), transportProperties);
dimensionedScalar t_laser_off ("t_laser_off", dimensionSet(0, 0, 1, 0, 0, 0, 0), transportProperties);
dimensionedScalar Dunit("Dunit", dimensionSet(0, 2, -1, 0, 0, 0, 0),transportProperties); // 1
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
dimensionedScalar L_Fe("L_Fe", dimensionSet(0, 2, -2, 0, 0, 0, 0), 247e+3); 
dimensionedScalar L_Cr("L_Cr", dimensionSet(0, 2, -2, 0, 0, 0, 0), 394e+3);
dimensionedScalar L = L_Fe*(1-(x_Cr/100)-(x_C/100)) + L_Cr*(x_Cr/100);

dimensionedScalar Lv_Fe("Lv_Fe", dimensionSet(0, 2, -2, 0, 0, 0, 0), 6.26e+6); 
dimensionedScalar Lv_Cr("Lv_Cr", dimensionSet(0, 2, -2, 0, 0, 0, 0), 6.62e+6); 
dimensionedScalar L_v = Lv_Fe*(1-(x_C/100)-(x_Cr/100)) + Lv_Cr*(x_Cr/100);

dimensionedScalar Tv_Fe("Tv_Fe", dimensionSet(0, 0, 0, 1, 0, 0, 0), 3134); 
dimensionedScalar Tv_Cr("Tv_Cr", dimensionSet(0, 0, 0, 1, 0, 0, 0), 2944); 
dimensionedScalar T_v = Tv_Fe*(1-(x_C/100)-(x_Cr/100)) + Tv_Cr*(x_Cr/100);

dimensionedScalar ST_limit ("ST_limit", dimensionSet(1, 0, -2, 0, 0, 0, 0), 1.48206459 - 0.81623e-3*(T_v.value()-Tl.value())); 
dimensionedScalar mu_limit("mu_limit", dimensionSet(1, -1, -1, 0, 0, 0, 0), 1e-3 * Foam::exp(-2.396 + 7950/T_v.value())); 
dimensionedScalar kl_limit("kl_limit", dimensionSet(1, 1, -3, -1, 0, 0, 0), ((35-0.3574*x_Cr.value())-0.0139*Tl.value()) + 0.0139*T_v.value()); 

dimensionedScalar M_Fe("M_Fe", dimensionSet(1, 0, 0, 0, -1, 0, 0), 0.055845);
dimensionedScalar M_C("M_C", dimensionSet(1, 0, 0, 0, -1, 0, 0), 0.012011);
dimensionedScalar M_Cr("M_Cr", dimensionSet(1, 0, 0, 0, -1, 0, 0), 0.051996);
dimensionedScalar M = 1/(((1-(x_C/100)-(x_Cr/100))/M_Fe)+((x_C/100)/(M_C))+((x_Cr/100)/M_Cr));
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "createPhi.H"

Info<< "Reading transportProperties\n" << endl;
Info<< "Latent heat of fusion L:" << L.value() << "J/kg"<< endl;
Info<< "Latent heat of vap L_v:" << L_v.value() << "J/kg"<< endl;
Info<< "Boiling point T_v:" << T_v.value() << "K"<< endl;
Info<< "Molar mass M:" << M.value() << "kg/mol"<< endl;

immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();


// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();


// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);


// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture)
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

// MULES compressed flux is registered in case scalarTransport FO needs it.
surfaceScalarField alphaPhiUn
(
    IOobject
    (
        "alphaPhiUn",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(phi.dimensions(), Zero)
);

#include "createMRF.H"
#include "createFvOptions.H"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField lf
(
    IOobject
    (
        "lf",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField FI
(
    IOobject
    (
        "FI",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("FI", dimensionSet(0, 0, 0, 0, 0, 0, 0), 1.0)                        
);

volScalarField FI_krate
(
    IOobject
    (
        "FI_krate",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("FI", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0)                        
);
volScalarField k_rate = FI_krate*k_init;  // just initializing
///////////////CEqn//////////////////////
volScalarField magGradAlpha = mag(fvc::grad(alpha1));

volScalarField DN
(
    IOobject
    (
        "DN",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("DN", dimensionSet(0, 2, -1, 0, 0, 0, 0), 0.0)  
);

dimensionedScalar DNs = Dunit;
dimensionedScalar DNl = Dunit;

volScalarField YN
(
    IOobject
    (
        "YN",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField Y
(
    IOobject
    (
        "Y",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField Ysat
(
    IOobject
    (
        "Ysat",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Ysat", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0)   
);

volScalarField lVol
(
    IOobject
    (
        "lVol",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField Yloss
(
    IOobject
    (
        "Yloss",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    k_rate * magGradAlpha * (Ysat - Y)
);

/////////////////////////////////////////
volScalarField ST = FI*ST_ini;  // just initializing
volScalarField cps = FI*cp2;           // just initializing
volScalarField cpl = FI*cp2;           // just initializing

volScalarField ks = FI*k2;  // just initializing
volScalarField kl = FI*k2;  // just initializing

volScalarField k = alpha1*(ks*(scalar(1)-lf)+kl*lf) +  alpha2*k2;
surfaceScalarField kappaf =  fvc::interpolate(k);

volScalarField rhoCp = alpha1*rho1*(cps*(scalar(1)-lf)+cpl*lf) +  alpha2*rho2*cp2;
surfaceScalarField rhoCpphi = fvc::interpolate(rhoCp)*phi;

volScalarField mu1 = FI*mu2;  // just initializing
volScalarField mu = alpha1*mu1 + alpha2*mu2;
surfaceScalarField muf = fvc::interpolate(mu);

volScalarField beta = alpha1*beta1 +  alpha2*beta2;
volScalarField A = -C * alpha1 * sqr(scalar(1)-lf)/ (pow(lf,scalar(3))+b) ;

// latent heat source terms (implicit)
volScalarField lfPRELOOP = lf;
volScalarField lfold = lf;
dimensionedScalar dFdT = 1 / (Tl - Ts);
volScalarField invF = lfold * (Tl - Ts) + Ts;
volScalarField dH = alpha1*rho * L * pos(Tl - T) * pos(T - Ts); //alpha1 multiplication
volScalarField Sp = - dH/runTime.deltaT() * dFdT;
volScalarField Sc = dH/runTime.deltaT() * (lfPRELOOP - lf + dFdT * invF);

///////////////////////////////////////////////////////////////////////////////// Radiation Model /////////////////////////////////////////////////////////////////////////////////////////////


volScalarField alphaT = alpha1;

autoPtr<radiation::radiationModel> radiation
(
        radiation::radiationModel::New(alphaT)
);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


