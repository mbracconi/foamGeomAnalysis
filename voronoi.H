	// Initialize file reader for dump
	t1 = runTime.elapsedCpuTime();
    std::ifstream file(dump_.c_str());
    std::string line;
	
	// Initialize bounding box for Voronoi
    double maxX = 1.e-6; double maxY = 1.e-6; double maxZ = 1.e-6;
    double minX = 1.e6; double minY = 1.e6; double minZ = 1.e6;
   
    // Initialize vector of centers
    List < vector > center;
    vector centert;

    label numCell = 0;
    while(true) 
    {
        // Read new line
        std::getline(file, line);
        std::istringstream iss(line);
        
        // Create vector to memorize info contained in dump_ file
        double val[5];
        
        // Read from dump_ file
        if (!(iss >> val[0]>>val[1]>>val[2]>>val[3]>>val[4])) 
        { 
            break; 
        } 

        maxX = std::max(maxX, val[1]); maxY = std::max(maxY, val[2]); maxZ = std::max(maxZ, val[3]);
        minX = std::min(minX, val[1]); minY = std::min(minY, val[2]); minZ = std::min(minZ, val[3]);
        centert[0] = val[1]; centert[1] = val[2]; centert[2] = val[3];
        center.append(centert);
        numCell++;
    }

    // Initialize Voronoi container
    voro::container con(minX-10,maxX+10,minY-10,maxY+10,minZ-10,maxZ+10,8,8,64,false,false,false,16);
 
    // Cycle on the total number of cells and add to the container
    // Define real number of cells added to the container
    unsigned int numCellInCon = 0;
    forAll(center, celli)
    {        
        if (con.point_inside(center[celli][0],center[celli][1],center[celli][2])) 
        {
            con.put(numCellInCon,center[celli][0],center[celli][1],center[celli][2]);numCellInCon++;
        }
    }

	// Compute Voronoi
	t1 = runTime.elapsedCpuTime();
	List < vector > center_list;
	voro::voronoicell vorocell;
	voro::c_loop_all cl(con);
	int id = 0;
	if(cl.start()) do if(con.compute_cell(vorocell,cl)) 
	{
		// Gather information about the computed Voronoi cell
		cl.pos(centert[0],centert[1],centert[2]);
		
		if(std::fabs(centert[0]) < std::fabs(maxBB.component(0)*1.5) && std::fabs(centert[1]) < std::fabs(maxBB.component(1)*1.5))
		{
		center_list.append(centert);
		}
		
	} while (cl.inc());
	
	Info << " * Initialized Voronoi decomposition in " << runTime.elapsedCpuTime() - t1 << " s" << nl << endl;
