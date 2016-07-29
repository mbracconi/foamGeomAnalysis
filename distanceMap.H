	// Create distance map
	if(createDistanceMap_)
	{
		Info << " * Computing distance map .. " << endl;
		t1 = runTime.elapsedCpuTime();
		scalar tmp;
		forAll (dist_map, i)
		{	
			// Initialize distances for cell I
			dist_map[i] = HUGE;
		
			// Loop over the boundary cells, evaluate distances and select the lower one
			forAll(boundaryCell,celli) 
			{
				tmp = 0.;
				for (direction a=0; a<vector::nComponents; a++)
				{
					tmp += std::fabs(C[i].component(a)-C[boundaryCell[celli]].component(a));					
				}
				dist_map[i] = std::min(tmp,dist_map[i]);
			}
			if (i%(label(0.2*dist_map.size())+1) == 0)
				Info << "  * Accomplished: " << i << "/" << dist_map.size() << endl;
		}
	    
	    forAll(dist_map.boundaryField(), patchi)
		{
			  dist_map.boundaryField()[patchi] == dist_map.boundaryField()[patchi].patchInternalField();
		}
		
		Info << " * Distance map evaluated in " << runTime.elapsedCpuTime() - t1 << " s" << nl << endl;
		t1 = runTime.elapsedCpuTime();
	}
	
