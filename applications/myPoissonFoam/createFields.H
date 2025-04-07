Info<< "Reading transportProperties\n" << endl;

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

dimensionedScalar alpha
(
    "alpha",                                  // Name of the scalar
    dimensionSet(0, 2, -1, 0, 0, 0, 0),       // Dimensions: [0 2 -1 0 0 0 0] (m^2/s for thermal diffusivity)
    transportProperties       // Lookup value from transportProperties
);

Info<< "Reading field Temperature\n" << endl;

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
