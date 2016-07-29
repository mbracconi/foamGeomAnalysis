	t1 = runTime.elapsedCpuTime();
	
	List < scalar > pore_diam_list;
	List < scalar > cell_diam_list;
	List < label >  cell_neib_count_list;
	List < scalar > cell_sphericityAR_list;
	List < scalar > cell_sphericityJR_list;
	
	//Info << " * Cells summary: " << endl;
	forAll(cells_list, i)
	{	
		// Initialize cell volume and center of gravity
		scalar cellVoli = 0.0;
		vector cog = vector::zero;
		// Evaluate cell volume and center of gravity
		forAll(cells_list[i].cell_cells, g)
		{
			cog += C[cells_list[i].cell_cells[g]]*V[cells_list[i].cell_cells[g]];
			cellVoli += V[cells_list[i].cell_cells[g]];
			if(!cells_list[i].boundCell)
				cells[cells_list[i].cell_cells[g]] = cells_list[i].id;
		}
		cog /= cellVoli;
		
		cells_list[i].vcell = cellVoli;
		cells_list[i].cog = cog;
		
		// Assign volume of cell to related field
		forAll(cells_list[i].cell_cells, g)
		{
			cell_volume[cells_list[i].cell_cells[g]] = cellVoli;
		}
		
		// Found Center of gravity of the cell
		scalar dcell = 1.0;	
		// Evaluate pore diameter by computing the geometric mean of all the distances from the cell center of gravity
		scalar dist = 0.0;
		forAll(cells_list[i].bound_cells, g)
		{
			cell_volume[cells_list[i].bound_cells[g]] = cellVoli;
			dist = Foam::mag( C[cells_list[i].bound_cells[g]] - cog);
			dcell *= Foam::pow( dist ,1./scalar(cells_list[i].bound_cells.size()));
			cells_list[i].dcell_min = std::min(cells_list[i].dcell_min, dist);
			cells_list[i].dcell_max = std::max(cells_list[i].dcell_max, dist);
		}
		dcell *= 2.0; 
		cells_list[i].dcell = dcell;
		
		// Evaluate sphericity according to Al-Roush
		forAll(cells_list[i].bound_cells, g)
		{
			cells_list[i].sphericityAR += Foam::mag( C[cells_list[i].bound_cells[g]] - cog)/cells_list[i].dcell_min - Foam::mag( C[cells_list[i].bound_cells[g]] - cog)/cells_list[i].dcell_max;
		}
		cells_list[i].sphericityAR /= scalar(cells_list[i].bound_cells.size());
		
		if(!cells_list[i].boundCell)
			cell_sphericityAR_list.append(cells_list[i].sphericityAR);
			
		// Maximum and minimum cell diameter
		cells_list[i].dcell_min *= 2.0;
		cells_list[i].dcell_max *= 2.0;
		
		// Find the cell containing the pore center
		label cellOfCellCoG = mesh.findCell(cog);

		// Approximate cell diameter ads geometric mean to the center of gravity and assign to cog cell
		cell_diam[cellOfCellCoG] = dcell;
		
		// Add to cell diameter list for post processing if it is not on the bound of the domain
		if(!cells_list[i].boundCell)
			cell_diam_list.append(dcell);
	
		if(!cells_list[i].boundCell)
			cell_neib_count_list.append(cells_list[i].neib_cells.size());
			
		// Evaluate pore characteristic
		forAll(cells_list[i].neib_cells, g)
		{
			// Find id of neighbouring cell
			label idx;
			forAll(cells_list, k)
			{
				if(cells_list[k].id == cells_list[i].neib_cells[g] )
				{
					idx = k;
					break;
				}
			}
					
			// Initialize pore volume and pore diameter
			scalar poreVoli = 0.;
			vector cog_p = vector::zero;
			scalar dpore = 1.0;
			label cont = 0;
				
			// Evaluate pore center (pore center of gravity and pore volume)
			// Loop over the boundary cell of the cells
			forAll(cells_list[i].bound_cells, p)
			{
				// Loop on the neib of the boundary cell
				forAll(neib[cells_list[i].bound_cells[p]], q)
				{
					// If the neib is owned by idx cells -> cell p is of this pore
					if( lab[neib[cells_list[i].bound_cells[p]][q]] == cells_list[idx].id )
					{
						cog_p += C[cells_list[i].bound_cells[p]]*V[cells_list[i].bound_cells[p]];
						poreVoli += V[cells_list[i].bound_cells[p]];
						if(!cells_list[i].boundCell)
							cells[cells_list[i].bound_cells[p]] = cells_list[i].id;
					}
				}
			}
			if(poreVoli > 0.)
			{
				cog_p /= poreVoli;
				
				cells_list[i].cog_pore.append(cog_p);
					
			
				// Evaluate an approximation of the pore diameter
				// Loop on the boundary cell
				forAll(cells_list[i].bound_cells, p)
				{
					// Loop on the neib of the boundary cell
					forAll(neib[cells_list[i].bound_cells[p]], q)
					{
						// If the neib is owned by idx cells -> cell p is of this pore
						if( lab[neib[cells_list[i].bound_cells[p]][q]] == cells_list[idx].id )
						{
							// By considering the geometric mean of the distances from the centery of gravity
							// to the patch reacting wall, the pore diameter can be evaluated 
							// Just reactingWall since they sorround a pore
							forAll(boundaryCell, b)
							{
								if( cells_list[i].bound_cells[p] == boundaryCell[b] )
								{
									dpore *= Foam::mag(  C[cells_list[i].bound_cells[p]] - cog_p);
									cont++;
									// Because cell p can be only one of the reactingWall ones
									break;
								}
							}
							
						}
					}
				}
				
				if(cont > 0)
				{
					// Evaluate pore diam
					dpore = 2.0*Foam::pow(dpore,1./scalar(cont));
					
					// Consider for post provessing only the ones are of id larger than actual
					if(cells_list[idx].id > cells_list[i].id)
					{
						// Find the cell containing the pore center
						label cellOfPoreCoG = mesh.findCell(cog_p);
										
						// Approximate pore diameter ads geometric mean to the center of gravity
						pore_diam[cellOfPoreCoG] = dpore;
						
						cells_list[i].dpore.append(dpore);
						
						// Consider for post processing only pores of non-boundary cells
						if(!cells_list[i].boundCell)
						{
							// Append pore diameter to pore list for post processing
							pore_diam_list.append(dpore);
						}
					}
					else
					{
						cells_list[i].dpore.append(0.0);
					}
					Info << "  *  * pore : " << cells_list[i].id << " <-> " <<cells_list[idx].id << " : " << dpore << endl;
				}
			}

		}
		Info << nl << endl;
		
	}
	
	// Assing boundary cells
	forAll(cells.boundaryField(), patchi)
	{
		  cells.boundaryField()[patchi] == cells.boundaryField()[patchi].patchInternalField();
	}
	
	forAll(lab.boundaryField(), patchi)
	{
		  lab.boundaryField()[patchi] == lab.boundaryField()[patchi].patchInternalField();
	}
	
	forAll(cell_volume.boundaryField(), patchi)
	{
		  cell_volume.boundaryField()[patchi] == cell_volume.boundaryField()[patchi].patchInternalField();
	}
	Info << " * Properties measured in " << runTime.elapsedCpuTime() - t1 << " s" << nl << endl;
