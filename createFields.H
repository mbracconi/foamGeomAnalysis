// Container of label of catchment basins 
    volScalarField lab
	(
		IOobject
		(
			"lab",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh,
		dimensionedScalar("lab", dimless, -3.0 )
	);
	
	volScalarField pore
	(
		IOobject
		(
			"pore",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh,
		dimensionedScalar("pore", dimless, -3.0 )
	);
	
	volScalarField pore_diam
	(
		IOobject
		(
			"pore_diam",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh,
		dimensionedScalar("pore_diam", dimless, -3.0 )
	);
	
	volScalarField cell_diam
	(
		IOobject
		(
			"cell_diam",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh,
		dimensionedScalar("cell_diam", dimless, -3.0 )
	);
	
	volScalarField cells
	(
		IOobject
		(
			"cells",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh,
		dimensionedScalar("cells", dimless, -3.0 )
	);
	
	volScalarField dist_map
	(
		IOobject
		(
			"dist_map",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh,
		dimensionedScalar("dist_map", dimless, -3.0 )
	);
	
	volScalarField cell_volume
	(
		IOobject
		(
			"cell_volume",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh,
		dimensionedScalar("cell_volume", dimless, 0.0 )
	);
	
	
