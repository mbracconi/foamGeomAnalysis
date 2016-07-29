// Assign to each cell the mesh points and label them
	forAll(lab, i)
	{
		scalar dist = HUGE;
		label id = -1;
		forAll(center_list, p)
		{
			if ( dist > Foam::mag(C[i] - center_list[p]) )
			{
				dist = Foam::mag(C[i] - center_list[p]);
				id = p;
			}
		}
		lab[i] = scalar(id);
		
		// Initialize neib
		labelList neibI;
		forAll(mesh.cellCells()[i],s)
		{
			// Cell is neib if id > 0 && < size and not itself
			if( mesh.cellCells()[i][s] >= 0 && mesh.cellCells()[i][s] < lab.size() && mesh.cellCells()[i][s] != i)
			{
				neibI.append(mesh.cellCells()[i][s]);
			}
		}
		neib[i] = neibI;
	}
	
	lab.write();
	
	Info << " * Cell segmented in " << runTime.elapsedCpuTime() - t1 << " s" << nl << endl;
