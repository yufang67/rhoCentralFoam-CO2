        volScalarField yPlusTemp
        (
            IOobject
            (
                "yPlusTemp",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("yPlusTemp", dimless, 0.0)
        );

        volScalarField yPlus
        (
            IOobject
            (
                "yPlus",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("yPlus", dimless, 0.0)
        );

        volVectorField uPlus
        (
            IOobject
            (
                "uPlus",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector("uPlus", dimless, vector::zero)
        );

        volScalarField uTau
        (
            IOobject
            (
                "uTau",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("uTau", dimVelocity, 0.0)
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
                IOobject::NO_WRITE
            ),
            mesh
        );
	Info<< "Reading field rho\n" << endl;
	volScalarField rho
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

	Info<< "Reading field mu, nut,muEff\n" << endl;
  
        volScalarField mu
        (
            IOobject
            (
                "thermo:mu",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        volScalarField muEff
        (
            IOobject
            (
                "muEff",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
//volVectorField rhoU
//(
//    IOobject
//    (
//        "rhoU",
//        runTime.timeName(),
//        mesh,
//        IOobject::NO_READ,
//        IOobject::NO_WRITE
//    ),
//    rho*U
//);

//surfaceScalarField phi("phi", fvc::flux(rhoU));


//autoPtr<CO2Thermo> swThermo
//(
//    CO2Thermo::New(mesh)
//);
//CO2Thermo& thermo = swThermo();

//autoPtr<compressible::turbulenceModel> turbulence
//(
//    compressible::turbulenceModel::New
//    (
//        rho,
//        U,
//        phi,
//        thermo
//    )
//);

volScalarField::Boundary d = nearWallDist(mesh).y();
