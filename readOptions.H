
	IOdictionary poreOptions
	(
		IOobject
		(
			"poreDict",
			runTime.time().constant(),
			mesh,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		)
	);
	

	string dump_(poreOptions.lookup("cellCenters"));
	word foamSurface_(poreOptions.lookup("foamSurface"));
	List<word>   inertBound_(poreOptions.lookup("inertBoundaries"));
	
	Switch createDistanceMap_(poreOptions.lookupOrDefault("createDistanceMap",false));

	Switch writeStats_(poreOptions.lookupOrDefault("writeStats",true));
	Switch writeStatsOnScreen_(poreOptions.lookupOrDefault("writeStatsOnScree",false));
