	// Creation of a list of the present open-cells with the boundary cell identified
	t1 = runTime.elapsedCpuTime();
	List < foam_cell > cells_list;
	forAll(center_list, i)
	{
		foam_cell cc;
		cc.id = scalar(i);
		bool boundMeshCells = false;
		// Loop on the mesh cells
		forAll(lab, p)
		{
			// Check if this cell is labeleb with the proper id
			if( std::fabs(lab[p] - scalar(i) ) < 1.e-4)
			{
				bool bound = false;
				
				// Loop on the neighbouring cells
				forAll(neib[p], q)
				{
					// If the label of neighbouring differs this is a boundary cell
					// and the neighbours can be stored
					if( std::fabs(lab[p] - lab[neib[p][q]]) > 1.e-4)
					{
						cc.bound_cells.append(p);
						cc.neib_cells.append(lab[neib[p][q]]);
						bound = true;
						// Assign correct value to pore field
						pore[p] = 1.0;
						break;
					}
				}
				
				// Add to boundary cell all the cells which share a face with the surface of the foam
				forAll(boundaryCell, b)
				{
					if( p == boundaryCell[b] )
						cc.bound_cells.append(p); break;
				}
				
				// Check if the cell is owned by one the inert boundary
				forAll(inertBoundaryCells, b)
				{
					if( p == inertBoundaryCells[b] )
					{
						boundMeshCells = true; break;
					}
				}
				
				// Add the cell to the core of the open-cell only if it is not owned by the foam surface
				if(!bound)
					cc.cell_cells.append(p);
			}
		}
		// If the cell exists
		if(cc.cell_cells.size() > 0)
		{
			// Assign the internal cell
			labelHashSet p;
			p.insert(cc.neib_cells);
			cc.neib_cells = p.toc();
			
			// Assign the boundary cell
			labelHashSet h;
			h.insert(cc.bound_cells);
			cc.bound_cells = h.toc();
			
			// Define if the cell is on the border of the domain
			cc.boundCell = boundMeshCells;
			
			// Initialize maximum and minimum dimension of cell diameter
			cc.dcell_max = SMALL;
			cc.dcell_min = HUGE;

			// Add to the cell list
			cells_list.append(cc);
			
		}
	}
	
	 cells.write();
	 Info << " * Identified cell features in " << runTime.elapsedCpuTime() - t1 << " s" << nl << endl;
