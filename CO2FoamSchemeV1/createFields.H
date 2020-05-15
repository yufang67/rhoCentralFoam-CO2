Info<< "Reading thermophysical properties\n" << endl;

autoPtr<CO2Thermo> swThermo
(
    CO2Thermo::New(mesh)
);
CO2Thermo& thermo = swThermo();
//Info<<thermo<<endl;
//volScalarField& e = thermo.he();
//*********************************
//
//
volScalarField& rho = thermo.rho();
//
Info<< "Reading initField e, rho, p, T, U\n" << endl;
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

//volScalarField rho
//(
//    IOobject
//    (
//        "rho",
//        runTime.timeName(),
//        mesh,
//        IOobject::MUST_READ,
//        IOobject::AUTO_WRITE
//    ),
//    mesh
//);

volScalarField e
(
    IOobject
    (
        "e",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
//    dimensionedScalar("e",dimensionSet(0,2,-2,0,0,0,0),0)
);


volVectorField rhoU
(
    IOobject
    (
        "rhoU",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    rho*U
);

volScalarField rhoE
(
    IOobject
    (
        "rhoE",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    rho*(e + 0.5*magSqr(U))
);

//volScalarField press
//(
//    IOobject
//   (
//        "press",
//        runTime.timeName(),
//        mesh,
//        IOobject::NO_READ,
//        IOobject::AUTO_WRITE
//    ),
//    mesh,
//    dimensionedScalar("press",dimensionSet(1,-1,-2,0,0,0,0),0)
//);

//volScalarField temp
//(
//    IOobject
//    (
//        "temp",
//        runTime.timeName(),
//        mesh,
//        IOobject::NO_READ,
//        IOobject::AUTO_WRITE
//    ),
//    mesh,
//    dimensionedScalar("temp",dimensionSet(0,0,0,1,0,0, 0),0)
//);
volScalarField c
(
    IOobject
    (
        "c",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("c",dimensionSet(0,1,-1,0,0,0, 0),0)
);
//volVectorField gradT
//(
//    IOobject
//    (
//        "gradT",
//        runTime.timeName(),
//        mesh,
//        IOobject::READ_IF_PRESENT,
//        IOobject::AUTO_WRITE
//    ),
//    fvc::grad(T)
//    mesh,
//    dimensionedScalar("gradT",dimensionSet(0,-1,0,-1,0,0, 0),0)
//);
volScalarField qual
(
    IOobject
    (
        "qual",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("qual",dimensionSet(0,0,0,0,0,0, 0),0)
);
volScalarField volfrac
(
    IOobject
    (
        "volfrac",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("volfrac",dimensionSet(0,0,0,0,0,0, 0),0)
);

volScalarField cv
(
    IOobject
    (
        "cv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("cv",dimensionSet(0,2,-2,-1,0,0,0),0)
);
volScalarField cp
(
    IOobject
    (
        "cp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
//    dimensionedScalar("cp",dimensionSet(0,0,0,0,0,0, 0),0)
    dimensionedScalar("cp",dimensionSet(0,2,-2,-1,0,0, 0),0)
);
volScalarField lambda_l
(
    IOobject
    (
        "lambda_l",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("lambda_l",dimensionSet(1,1,-3,-1,0,0,0),0)
);
volScalarField lambda_t
(
    IOobject
    (
        "lambda_t",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("lambda_t",dimensionSet(1,1,-3,-1,0,0,0),0)
);
volScalarField lambda_eff
(
    IOobject
    (
        "lambda_eff",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("lambda_eff",dimensionSet(1,1,-3,-1,0,0,0),0)
);
volScalarField gamma
(
    IOobject
    (
        "gamma",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("gamma",dimensionSet(0,0,0,0,0,0, 0),0)
);

volScalarField dedr
(
    IOobject
    (
        "dedr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("dedr",dimensionSet(-1,5,-2,0,0,0,0),0)
);
//

surfaceScalarField pos
(
    IOobject
    (
        "pos",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("pos", dimless, 1.0)
);

surfaceScalarField neg
(
    IOobject
    (
        "neg",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("neg", dimless, -1.0)
);

surfaceScalarField phi("phi", fvc::flux(rhoU));

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);
//Info<<turbulence<<endl;
