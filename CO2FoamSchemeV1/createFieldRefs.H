volScalarField& p = thermo.p();
volScalarField& T = thermo.T();
//volScalarField& rho = thermo.rho();

volScalarField& mu= thermo.mu();
//const volScalarField& nu= thermo.nu();
//volScalarField& rPr= thermo.rPr();
//volScalarField& alpha= thermo.alphah();

//bool inviscid(true);
bool inviscid(false);
//if (max(mu.primitiveField()) > 0.0)
//if (gMax(mu) > 0.0)
//{
//    inviscid = false;
//}
//Info<<"invisid "<<inviscid<<endl;
