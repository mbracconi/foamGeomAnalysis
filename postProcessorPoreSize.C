/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Author                                                               
    Mauro Bracconi <mauro.braccon@polimi.it>                             
    Department of Energy                                                  
    Politecnico di Milano                                                
    via La Masa, 34 - 20156 - Milano, Italy   
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    foamGeomAnalisys

Description
    Utility to evaluate the geometrical properties of an open-cell foam
    by segmenting the mesh. The cell are identified exploiting the Voronoi
    tessellation (requires voro++ library) and analysis on the pore and cell
    sizes, sphericity, neighbouring cells are performed
    

\*---------------------------------------------------------------------------*/

// Voro++ library
#include "voro++.hh"

// OpenFOAM
#include <iomanip>
#include <iostream>
#include <fstream>
#include "IOField.H"
#include "fvCFD.H"
#include "meshTools.H"
#include "fvMesh.H"
#include "CPCCellToCellStencil.H"
#include "IOobjectList.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

struct foam_cell
{
	// open-cell id
	scalar id;
	
	// computational cell owned by the open-cell
	labelList cell_cells;
	
	// computational cell owned by the open-cell boundary
	labelList bound_cells;
	
	// open-cell neighbouring
	labelList neib_cells;
	
	// cell diameter and maximinum and minimun size
	scalar dcell;
	scalar dcell_max;
	scalar dcell_min;
	
	// cell volume
	scalar vcell;
	
	// list of pore diameter
	List<scalar> dpore;
	// location of the cell center
	vector cog;
	// location of the pore centers
	List<vector> cog_pore;
	
	// Switch assessing if cell is located at the boundary of the domain
	Switch boundCell;
	
	// sphericityAR index
	scalar sphericityAR;
};

template <class T>
scalar average(List < T > l)
{
	scalar ave = 0.0;
	forAll(l, i)
	{
		ave += l[i];
	}
	return ave / scalar(l.size());
}

template <class T>
scalar dev_std(List < T > l, scalar mean)
{
	scalar dev = 0.0;
	forAll(l, i)
	{
		dev += Foam::pow(l[i]-mean, 2.0);
	}
	return Foam::sqrt(dev / scalar(l.size()));
}

int main(int argc, char *argv[])
{
	// Solver setup (folders, mesh, etc.)
	#include "setRootCase.H"
	#include "createTime.H"

	runTime.functionObjects().off();
	
	// Read mesh and prepare domain
	scalar t1 = runTime.elapsedCpuTime();
	
	#include "createMesh.H"
	#include "createFields.H"
	#include "readOptions.H"
	
	// Define cell centers and volumes
	const volVectorField& C = mesh.C();
	const scalarField& V = mesh.V();
	
	// Define patch reactingWall
	label patchID = mesh.boundaryMesh().findPatchID(foamSurface_);
	const labelList boundaryCell = mesh.boundaryMesh()[patchID].faceCells();
	
	// Define boundary inert cells 
	labelList inertBoundaryCells;
	labelHashSet bCells;
	forAll(inertBound_, i)
	{
		label patchIDB = mesh.boundaryMesh().findPatchID(inertBound_[i]);
		bCells.insert(mesh.boundaryMesh()[patchIDB].faceCells());
	}
	inertBoundaryCells = bCells.toc();
	
	// Define bounding boc
	const boundBox& boundBox = mesh.bounds();
	Vector<scalar> minBB = boundBox.min();
	Vector<scalar> maxBB = boundBox.max();
	
	// Initialize neibghrood list
	List < labelList > neib;
	neib.setSize(C.size());
	CPCCellToCellStencil wideStencil(mesh);
	
	Info << " * Mesh readed in " << runTime.elapsedCpuTime() - t1 << " s" << nl << endl;

	// Compute Voronoi 
	#include "voronoi.H"
	
	// Label all cells with the proper cell id
	#include "labeling.H"	
	
	// Assign the mesh cell to the proper open-cell
	#include "openCellAssignment.H"
	
	// Evaluate geometrical properties
	#include "geomEvaluation.H"
	
	// Evaluate if requested distance map
	#include "distanceMap.H"
	
	// Write stats 
	scalar mean, devstd;
	
	// On screen
	if(writeStatsOnScreen_)
	{

		Info << " * Pore size distribution   : "<<endl;
		forAll(pore_diam_list,i)
		{
			Info << "  * pore : " << i+1 << " : " << pore_diam_list[i] << " [m]" << endl;
		}
		mean = average(pore_diam_list);
		devstd = dev_std(pore_diam_list,mean);
		Info << " * Average pore size        : " << mean << " [m]" << endl;
		Info << " * Dev std pore size        : " << devstd << " [m]" << endl;
		Info << nl << endl;
		
		Info << " * Cell size distribution   : "<<endl;
		forAll(cell_diam_list,i)
		{
			Info << "  * cell : " << i+1 << " : " << cell_diam_list[i] << endl;
		}
		mean = average(cell_diam_list);
		devstd = dev_std(cell_diam_list,mean);
		Info << " * Average cell size        : " << mean << " [m]" << endl;
		Info << " * Dev std cell size        : " << devstd << " [m]" << endl;
		Info << nl << endl;
		
		Info << " * Neighbouring distribution: "<<endl;
		forAll(cell_neib_count_list,i)
		{
			Info << "  * cell : " << i+1 << " : " << cell_neib_count_list[i] << endl;
		}
		mean = average(cell_neib_count_list);
		devstd = dev_std(cell_neib_count_list,mean);
		Info << " * Average neighbouring     : " << mean << " [-]" << endl;
		Info << " * Dev std neighbouring     : " << devstd << " [-]" << endl;
		Info << nl << endl;
		
		Info << " * sphericityAR distribution  : "<<endl;
		forAll(cell_sphericityAR_list,i)
		{
			Info << "  * cell : " << i+1 << " : " << cell_sphericityAR_list[i] << endl;
		}
		mean = average(cell_sphericityAR_list);
		devstd = dev_std(cell_sphericityAR_list,mean);
		Info << " * Average sphericityAR       : " << mean << " [-]" << endl;
		Info << " * Dev std sphericityAR       : " << devstd << " [-]" << endl;
		Info << nl << endl; 
	}
	
	// On file
	if(writeStats_)
	{
		autoPtr<OFstream> fmassFlux_;

		fmassFlux_.reset( new OFstream("pore_distribution"));
		fmassFlux_() << "# Pore count \t diameter [m]" <<  endl;
		forAll(pore_diam_list,i)
		{
			fmassFlux_() << i+1 << "\t" << pore_diam_list[i] << endl;
		}
		
		
		fmassFlux_.reset( new OFstream("cell_distribution"));
		fmassFlux_() << "# cell count \t diameter [m]" << endl;
		forAll(cell_diam_list,i)
		{
			fmassFlux_() << i+1 << "\t" << cell_diam_list[i]  << endl;
		}

		fmassFlux_.reset( new OFstream("neighbouring_distribution"));
		fmassFlux_() << "# cell count \t neighbouring count [-]" << endl;
		forAll(cell_neib_count_list,i)
		{
			fmassFlux_() << i+1 << "\t" << cell_neib_count_list[i]  << endl;
		}
		
		fmassFlux_.reset( new OFstream("sphericityAR_distribution"));
		fmassFlux_() << "# cell count \t sphericityAR [-]" << endl;
		forAll(cell_sphericityAR_list,i)
		{
			fmassFlux_() << i+1 << "\t" << cell_sphericityAR_list[i]  << endl;
		}
		
		fmassFlux_.reset( new OFstream("cell_summary"));
		fmassFlux_() << "# cell summary" << endl;
		forAll(cells_list,i)
		{
			fmassFlux_() << " ################################################################################################# " << endl;
			fmassFlux_() << "  * cell i               : " << i + 1 << "/" << cells_list.size() <<  endl;
			fmassFlux_() << "  * id                   : " << cells_list[i].id << endl;
			fmassFlux_() << "  * center of gravity    : ( " << cells_list[i].cog[0] << " " << cells_list[i].cog[1] << " " << cells_list[i].cog[2] << " )" << endl;
			fmassFlux_() << "  * average cell diameter: " << cells_list[i].dcell << " [m]" << endl;
			fmassFlux_() << "  * minimum cell diameter: " << cells_list[i].dcell_min << " [m]" << endl;
			fmassFlux_() << "  * maximum cell diameter: " << cells_list[i].dcell_max << " [m]" << endl;
			fmassFlux_() << "  * cell volume          : " << cells_list[i].vcell << " [m3]" << endl;
			fmassFlux_() << "  * number conn cells    : " << cells_list[i].neib_cells.size() << endl;
			fmassFlux_() << "  * connected cells id   : (";
			forAll(cells_list[i].neib_cells, g)
			{
				fmassFlux_() << " " << cells_list[i].neib_cells[g];
			}
			fmassFlux_() << " )" << endl << endl;
			fmassFlux_() << "  * pore information     : " << endl;
			forAll(cells_list[i].dpore,g)
			{
				fmassFlux_() << "    * pore count     : " << g+1 << endl;
				fmassFlux_() << "    * pore cog       : ( " << cells_list[i].cog_pore[g][0] << " " << cells_list[i].cog_pore[g][1] << " " << cells_list[i].cog_pore[g][2] << " )" << endl;
				fmassFlux_() << "    * pore diameter  : " << cells_list[i].dpore[g] << " [m]" << endl << endl;
			}
		}
		
		fmassFlux_.reset( new OFstream("stats"));
		mean = average(pore_diam_list);
		devstd = dev_std(pore_diam_list,mean);
		fmassFlux_() << " * Average pore size        : " << mean << " [m]" << endl;
		fmassFlux_() << " * Dev std pore size        : " << devstd << " [m]" << endl;
		fmassFlux_() << " * Min pore size            : " << Foam::min(pore_diam_list) << " [m]" << endl;
		fmassFlux_() << " * Max pore size            : " << Foam::max(pore_diam_list) << " [m]" << endl << endl;
		mean = average(cell_diam_list);
		devstd = dev_std(cell_diam_list,mean);
		fmassFlux_() << " * Average cell size        : " << mean << " [m]" << endl;
		fmassFlux_() << " * Dev std cell size        : " << devstd << " [m]" << endl;
		fmassFlux_() << " * Min cell size            : " << Foam::min(cell_diam_list) << " [m]" << endl;
		fmassFlux_() << " * Max cell size            : " << Foam::max(cell_diam_list) << " [m]" << endl << endl;
		mean = average(cell_neib_count_list);
		devstd = dev_std(cell_neib_count_list,mean);
		fmassFlux_() << " * Average neighbouring     : " << mean << " [-]" << endl;
		fmassFlux_() << " * Dev std neighbouring     : " << devstd << " [-]" << endl;
		fmassFlux_() << " * Min neighbouring         : " << Foam::min(cell_neib_count_list) << " [-]" << endl;
		fmassFlux_() << " * Max neighbouring         : " << Foam::max(cell_neib_count_list) << " [-]" << endl << endl;
		mean = average(cell_sphericityAR_list);
		devstd = dev_std(cell_sphericityAR_list,mean);
		fmassFlux_() << " * Average sphericityAR       : " << mean << " [-]" << endl;
		fmassFlux_() << " * Dev std sphericityAR       : " << devstd << " [-]" << endl;
		fmassFlux_() << " * Min sphericityAR           : " << Foam::min(cell_sphericityAR_list) << " [-]" << endl;
		fmassFlux_() << " * Max sphericityAR           : " << Foam::max(cell_sphericityAR_list) << " [-]" << endl << endl;
	}

	
	// Write fields
	cells.write();
	cell_diam.write();
	pore_diam.write();
	cell_volume.write();
	pore.write();
	dist_map.write();
	
	forAll(cells_list, i)
	{
		cells_list[i].dpore.clear();
		cells_list[i].cog_pore.clear();
		cells_list[i].cell_cells.clear();
		cells_list[i].bound_cells.clear();
		cells_list[i].neib_cells.clear();
	}

	cells_list.clear();

	Info << "End" << endl;
	return 0;
}

// ************************************************************************* //
