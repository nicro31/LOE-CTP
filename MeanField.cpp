#include "MeanField.h"

using namespace std;


double adjustConcentration(double &pTemp, Graph &g)
{
    // Adjust chemical potential to get desired carrier concentration
    double cP(0); //Chemical potential (eV)
    double aD,bD;
    double tol(1e-10); // Tolerance
    double errPTmp(0);
    double errP(-1e10);

    if( pTemp < g.m_concentration ) { aD = 0; bD = 20;}
    else { aD = -20, bD = 0;}
    while( (bD-aD) / 2 > tol )
    {
        pTemp = 0;
        cP = (aD + bD) / 2;

        for (int v(0) ; v < g.m_site.size() / g.m_nbrLevels / 8 ; v++)
        {
            for(int u(0) ; u < g.m_nbrLevels ; u++)
            {
                int i(u*g.m_site.size()/g.m_nbrLevels + g.m_orderMeanField[v]);
                pTemp += ( 1.0 / ( 1 + exp( (-g.m_site[i].m_energy - g.m_site[i].m_cP - cP) / (kB*g.m_temperature) ) ) ) ;

            }
        }
        pTemp /= g.m_site.size() / 8;
        //std::cout << pTemp << "     " << g.m_concentration << '\n';
        if (pTemp < g.m_concentration ) { aD = cP;}
        else { bD = cP;}

    }

    for(int i(0) ; i < g.m_site.size() ; i++)
    {
        g.m_site[i].m_cP += cP;
        g.m_site[i].m_occup = 1.0 / ( 1 + exp( (-g.m_site[i].m_energy - g.m_site[i].m_cP) / (kB*g.m_temperature) ) );
        //g.m_site[i].m_occupT = g.m_site[i].m_occup;
        errPTmp = fabs(g.m_site[i].m_occup - g.m_site[i].m_occupT)/g.m_site[i].m_occupT;
        if(errPTmp > errP)
        {
            errP = errPTmp;

        }


    }

    return errP;
}

double computeMobility(Graph &g, bool usePrevDistrib, std::string folder, int u)
{


    double relaxFactor(1.0);

    cout << "Computing mobility...";

    /// Pre-compute r.F to speed up calculation ///
    for (int i=0 ; i < g.m_site.size() ; i ++)
    {
        g.m_site[i].m_displacementField.resize(0);
        for (int j(0) ; j < g.m_site[i].m_neighbor.size() ; j ++)
        {
            g.m_site[i].m_displacementField.push_back(g.m_site[i].m_neighborDisplacement[j]*g.m_field);
        }
    }


    /// Init system sites occupancies (first time we compute mobility)///
    double pTemp(0);
    if ( !usePrevDistrib )
    {
        //cout << '\n' << "Init system properties...";

        for (int i(0) ; i < g.m_site.size() ; i ++)
        {
            // Set occupation probability
            g.m_site[i].m_occup = 1.0 / ( 1 + exp( (-g.m_site[i].m_energy ) / (kB*g.m_temperature) ) );

            g.m_site[i].m_occupT = g.m_site[i].m_occup;
            pTemp += g.m_site[i].m_occup / g.m_site.size();
        }

        //cout << "ok" << endl;
    }

    /// Init system sites occupancies (using previous occupancy distribution as the starting point)///
    if ( usePrevDistrib )
    {
        //cout << '\n' << "Init system properties...";

        for (int i(0) ; i < g.m_site.size() ; i ++)
        {
            // Set occupation probability
            g.m_site[i].m_occup = g.m_site[i].m_occupT;
            pTemp += g.m_site[i].m_occup / g.m_site.size();
        }

        //cout << "ok" << endl;
        //cout << "pTemp = " << pTemp << endl;
    }

    /// Adjust chemical potential to get the desired carrier concentration (first time we compute mobility)///
    if( !usePrevDistrib )
    {
        //cout << "Adjust chemical potential..." << endl;
        double cP(0); //Chemical potential (eV)
        double aD,bD;
        double tol(1e-10); // Tolerance

        if( pTemp < g.m_concentration ) { aD = 0; bD = 20;}
        else { aD = -20, bD = 0;}

        while( (bD-aD) / 2.0 > tol )
        {
            pTemp = 0;

            cP = (aD + bD) / 2;

            for (int i(0) ; i < g.m_site.size() ; i ++)
            {


                g.m_site[i].m_occup = 1.0 / ( 1 + exp( (-g.m_site[i].m_energy - cP) / (kB*g.m_temperature) ) );

                g.m_site[i].m_occupT = g.m_site[i].m_occup;
                g.m_site[i].m_cP = cP;
                pTemp += g.m_site[i].m_occup / g.m_site.size();
            }

            if (pTemp < g.m_concentration ) { aD = cP;}
            else { bD = cP;}



        }
    }


    /// Adjust chemical potential to get the desired carrier concentration (using previous occupancy distribution as the starting point)///
    if( usePrevDistrib )
    {
        double cP(0); //Chemical potential (eV)
        double aD,bD;
        double tol(g.m_concentration * 1e-6); // Tolerance

        if(fabs(g.m_concentration-pTemp) > tol)
        {
            if( pTemp < g.m_concentration ) { aD = 0; bD = 20;}
            else { aD = -20, bD = 0;}

            while( (bD-aD) / 2 > eps )
            {
                pTemp = 0;

                cP = (aD + bD) / 2;

                for (int i(0) ; i < g.m_site.size() ; i ++)
                {
                    pTemp += ( 1.0 / ( 1 + exp( (-g.m_site[i].m_energy - g.m_site[i].m_cP - cP) / (kB*g.m_temperature) ) ) ) / g.m_site.size();
                }

                if (pTemp < g.m_concentration ) { aD = cP;}
                else { bD = cP;}
            }

        }
        for(int i(0) ; i < g.m_site.size() ; i++)
        {
            g.m_site[i].m_cP += cP;
            g.m_site[i].m_occup = 1.0 / ( 1 + exp( (-g.m_site[i].m_energy - g.m_site[i].m_cP) / (kB*g.m_temperature) ) );
            g.m_site[i].m_occupT = g.m_site[i].m_occup;

        }
    }


    /// Solve equilibrium probability distribution ///
    //cout << "Solve master equation..." << endl;
    double accuracyMob(1e-9);
    double accuracyOcc(1e-8);
    double accuracyOccGlobal(1e-9);
    double err(1e6);
    double mobilityF(0);
    double errOccup(1e6);
    int nbIterTot(0);
    int nbIterUpdate(1e4); // Number of iterations after which mobility is recomputed
    int nbIterMax(1e8);
    string fName;
    if(u == 0) {fName = "X";}
    if(u == 1) {fName = "Y";}
    if(u == 2) {fName = "Z";}
    std::ofstream errFlux(folder + fName + "_convergence_distrib_F" + std::to_string(g.m_field.norm() * 1e8) + "_T" + std::to_string(g.m_temperature) + "_P" + std::to_string(g.m_concentration) + ".txt");
    errFlux << "Number of iterations" << "    " << "Error on mobility" << "     "  << "Error on occupancy" << "     "  << "Mobility" << "    " << "MeanOccup" << endl;

    do
    {


        pTemp = 0;
        double errPTmp(0);
        double errP(-1e10);

        // Solve master equation
        for (int v(0) ; v < g.m_site.size() / g.m_nbrLevels / 8 ; v++)
        {
            for(int u(0) ; u < g.m_nbrLevels ; u++)
            {
                int i(u*g.m_site.size()/g.m_nbrLevels + g.m_orderMeanField[v]);
                double term1(0);
                double term2(0);


                for (int j(0) ; j < g.m_site[i].m_neighbor.size() ; j ++)
                {
                        int k(g.m_site[i].m_neighborIdx[j]);

                        // Direct update of occupancy during iterations
                        term1 += g.m_site[k].m_occup * g.m_site[i].m_hopFrom[j];
                        term2 += (1-g.m_site[k].m_occup) * g.m_site[i].m_hopTo[j];


                        // No update during iterations
                        //term1 += g.m_site[k].m_occupT * g.m_site[i].m_hopFrom[j];
                        //term2 += (1-g.m_site[k].m_occupT) * g.m_site[i].m_hopTo[j];
                }

                g.m_site[i].m_occupT = g.m_site[i].m_occup;
                g.m_site[i].m_occup =  relaxFactor * (term1 / (term1 + term2)) + (1.0 - relaxFactor) * g.m_site[i].m_occupT;

                // Calculate the error on occupation
                errPTmp = fabs(g.m_site[i].m_occup - g.m_site[i].m_occupT)/g.m_site[i].m_occupT;
                if(errPTmp > errP) {errP = errPTmp;}
                g.m_site[i].m_occupT = g.m_site[i].m_occup;

                g.m_site[i].m_cP = -( log(1.0 / g.m_site[i].m_occup - 1) * kB * g.m_temperature - (-g.m_site[i].m_energy) );
                pTemp += g.m_site[i].m_occup;
                for(int w(1) ; w <8 ; w ++)
                {
                    int q(i+w*(g.m_site.size() / g.m_nbrLevels / 8));
                    g.m_site[q].m_occup = g.m_site[i].m_occup;
                    g.m_site[q].m_cP = g.m_site[i].m_cP;
                    g.m_site[q].m_occupT = g.m_site[q].m_occup;
                }
            }
        }
        errOccup = errP;
        pTemp /= g.m_site.size() / 8;

        /*for(int h(0) ; h < g.m_site.size() ; h++)
        {
            g.m_site[h].m_occupT = g.m_site[h].m_occup;
        }*/


        // Adjust concentration
        if(fabs(g.m_concentration-pTemp) > accuracyOccGlobal)
        {
            errOccup = adjustConcentration(pTemp,g);
        }


        nbIterTot ++;

        if(nbIterTot % nbIterUpdate == 0)
        {
            // Adjust concentration
            //adjustConcentration(pTemp,g);

            // Compute mobility
            double mobility(0);
            errPTmp = 0;
            errP = -1e10;
            for (int v(0) ; v < g.m_site.size() / g.m_nbrLevels / 8 ; v++)
            {
                for(int u(0) ; u < g.m_nbrLevels ; u++)
                {
                    int i(u*g.m_site.size()/g.m_nbrLevels + g.m_orderMeanField[v]);
                    for (int j(0) ; j < g.m_site[i].m_neighbor.size() ; j ++)
                    {
                            int k(g.m_site[i].m_neighborIdx[j]);
                            mobility += g.m_site[i].m_occup * ( 1 - g.m_site[k].m_occup ) * g.m_site[i].m_hopTo[j] * g.m_site[i].m_displacementField[j];
                    }

                    // Calculate the error on occupation
                    //errPTmp = fabs(g.m_site[i].m_occup - g.m_site[i].m_occupT)/g.m_site[i].m_occupT;
                    //if(errPTmp > errP) {errP = errPTmp;}

                }
            }
            //errOccup = errP;
            mobility /= pTemp * (g.m_site.size() / 8.0) * ( g.m_field.norm() );
            err = fabs(mobility - mobilityF) / fabs(mobilityF) / (double) nbIterUpdate;
            if (mobilityF == 0) {err = 1e10;}
            mobilityF = mobility;
            errFlux << nbIterTot << "   " << err  << "    " << errOccup << "      " << mobilityF << "       " << pTemp << endl ;

        }

    }while ( (err > accuracyMob || errOccup > accuracyOcc || fabs(g.m_concentration-pTemp) > accuracyOccGlobal) && nbIterTot <= nbIterMax );


    errFlux.close();
    cout << "Done (Mobility is " << mobilityF << ")" << endl;
    return mobilityF;
}

Graph createCubic(int dim, double latticeParam)
{
    Graph g;
    g.m_system.m_boxSize = (dim)*latticeParam;
    g.m_system.m_topo.m_chainLength = dim;
    g.m_system.m_topo.m_pedotNbr = dim*dim;
    g.m_system.m_topo.m_tosNbr = 0;


    int counter(0);

    /// Create original cell ///
    for (int i(0) ; i < dim ; i ++)
    {
        for (int j(0) ; j < dim ; j++)
        {
            Polymer p;
            p.m_chainLength = dim;

            for (int k(0) ; k < dim ; k++)
            {
                RigidFragment rf;
                Vect3 center;
                rf.m_type = edot;
                rf.m_molNbr = counter;
                rf.m_fragNbr = k;


                // Set position
                center.m_x = i*latticeParam;
                center.m_y = j*latticeParam;
                center.m_z = k*latticeParam;
                rf.m_center = center;

                p.fragment.push_back(rf);

            }

            counter ++;
            g.m_system.m_pedot.push_back(p);
        }
    }


    /// Extend cell by pbc ///
    for (int a(0) ; a < 2 ; a ++ )
    {
        for (int b(0) ; b < 2 ; b ++)
        {
            for (int c(0) ; c < 2 ; c ++)
            {
                if(a!= 0 || b != 0 || c != 0)
                {
                    for (int i(0) ; i < dim ; i ++)
                    {
                        for (int j(0) ; j < dim ; j++)
                        {
                            Polymer p;
                            p.m_chainLength = dim;

                            for (int k(0) ; k < dim ; k++)
                            {
                                RigidFragment rf;
                                Vect3 center;
                                rf.m_type = edot;
                                rf.m_molNbr = counter;
                                rf.m_fragNbr = k;


                                // Set position
                                center.m_x = i*latticeParam + a * g.m_system.m_boxSize;
                                center.m_y = j*latticeParam + b * g.m_system.m_boxSize;
                                center.m_z = k*latticeParam + c * g.m_system.m_boxSize;
                                rf.m_center = center;

                                p.fragment.push_back(rf);

                            }

                            counter ++;
                            g.m_system.m_pedot.push_back(p);
                        }
                    }

                }
            }
        }
    }
    g.m_system.m_boxSize = 2 * (dim)*latticeParam;
    g.m_system.m_topo.m_chainLength = dim;
    g.m_system.m_topo.m_pedotNbr = 8 * dim*dim;
    g.m_system.m_topo.m_tosNbr = 0;



    g.m_system.m_line1 = "TITLE     vmdmolecule0 t= 29000.00000" ;
    g.m_system.m_line2 = "REMARK    THIS IS A SIMULATION BOX" ;
    g.m_system.m_line3 = "       20       20  90.00  90.00  90.00 P 1           1" ;
    g.m_system.m_line4 = "MODEL        1" ;

    return g;

}


Vect3 findPercolationThreshold(Graph &g, std::string folder)
{
    // All the values in unit log10(H*1000) H in eV
    double threshMin(-30.0);
    double threshMax(5.0);
    double threshErr(0.01);
    Vect3 threshVect;
    Vect3 clusterIdVect;

    for (int i(0) ; i < 3 ; i++)
    {

        double aD(threshMin);
        double bD(threshMax);
        double thresh(0);
        int clusterId;
        bool isInf(false);

        while( (bD-aD) / 2.0 > threshErr || !isInf)
        {

            thresh = (bD+aD) / 2.0;
            std::cout << "threshold: " << thresh << '\n';
            double threshMod = (exp(thresh*log(10))/1000.0);

            g.cluster(threshMod);
            g.percolation();

            double clusterMax(0);

            for(int j(0) ; j < g.m_clusterNb ; j ++)
            {
                double cSize(0);
                switch(i)
                {
                    case 0 : cSize = g.m_clusterDim[j].m_x; break;
                    case 1 : cSize = g.m_clusterDim[j].m_y; break;
                    case 2 : cSize = g.m_clusterDim[j].m_z; break;

                }

                if(cSize > clusterMax)
                {
                    clusterMax = cSize;
                    clusterId = j;
                }
            }

            cout << clusterMax << '\n';

            if ( clusterMax == 1000 ) { aD = thresh; isInf = true;}
            else { bD = thresh; isInf = false;}
        }

        switch(i)
        {
            case 0 : threshVect.m_x = thresh; clusterIdVect.m_x = clusterId; break;
            case 1 : threshVect.m_y = thresh; clusterIdVect.m_y = clusterId; break;
            case 2 : threshVect.m_z = thresh; clusterIdVect.m_z = clusterId; break;

        }

        // Save result of the percolated graph
        string graphPdb(folder+"graph"+std::to_string(i) +".pdb");
        string graphPsf(folder+"graph"+std::to_string(i) +".psf");
        string systemClustered("system_cluster"+std::to_string(i)+".pdb");
        g.writeGraph(graphPdb,graphPsf,false);
        g.writeSystemClustered(systemClustered,folder);
        g.writeCurrent(folder + "occup"+std::to_string(i) +".pdb",folder + "current"+std::to_string(i) +".pdb",folder + "current"+std::to_string(i) +".psf",false);

    }

    // Write the result in a file
    ofstream os(folder + "percolationThreshold.txt");
    if(os.is_open())
    {
        os << "X" << "  " << "Y" << "   " << "Z" << '\n';
        os << setprecision(15) << threshVect.m_x << " " << threshVect.m_y << "    " << threshVect.m_z << '\n';
        os << clusterIdVect.m_x << "    " << clusterIdVect.m_y << "     " << clusterIdVect.m_z << '\n';
        os.close();
    }
    else
    {
        std::cerr << "Can't open percolation threshold file !" << '\n';
    }

    std::cout << "Percolation thresholds: " << threshVect.m_x << "  " << threshVect.m_y << "    " << threshVect.m_z << '\n';
    return threshVect;

}

double findPercolationThresholdN(Graph &g, std::string folder)
{
    // All the values in unit log10(H*1000) H in eV
    double threshMin(-30.0);
    double threshMax(5.0);
    double threshErr(0.01);

    double aD(threshMin);
    double bD(threshMax);
    double thresh(0);
    int clusterId;
    bool isInf(false);

    while( (bD-aD) / 2.0 > threshErr || !isInf)
    {

        thresh = (bD+aD) / 2.0;
        std::cout << "threshold: " << thresh << '\n';
        double threshMod = (exp(thresh*log(10))/1000.0);

        g.cluster(threshMod);
        double clusterMax(0);

        for(int j(0) ; j < g.m_clusterNb ; j ++)
        {
            if(g.m_clusterSiteNb[j] > clusterMax)
            {
                clusterMax = g.m_clusterSiteNb[j];
                clusterId = j;
            }
        }

        cout << clusterMax << '\n';

        if ( clusterMax >= g.m_site.size() / 2.0 ) { aD = thresh; isInf = true;}
        else { bD = thresh; isInf = false;}
    }


    // Save result of the percolated graph
    string graphPdb(folder+"graph" +".pdb");
    string graphPsf(folder+"graph" +".psf");
    string systemClustered("system_cluster.pdb");
    g.writeGraph(graphPdb,graphPsf,false);
    g.writeSystemClustered(systemClustered,folder);
    g.writeCurrent(folder + "occup" +".pdb",folder + "current"+".pdb",folder + "current"+".psf",false);

    // Write the result in a file
    ofstream os(folder + "percolationThreshold.txt");
    if(os.is_open())
    {
        os << "Threshold" << "  " << "ID" << '\n';
        os << setprecision(15) << thresh << " " << clusterId << '\n';
        os.close();
    }
    else
    {
        std::cerr << "Can't open percolation threshold file !" << '\n';
    }

    std::cout << "Percolation threshold: " << thresh << "   cluster Id: " << clusterId << '\n';
    return thresh;
}

void testCubic()
{
   Graph g = createCubic(5,1);

    /// Nodes segmentation parameters ///
    g.m_segmentLength = 1;

    /// Edges threshold parameters ///
    g.m_dThresholdInter = sqrt(3);
    g.m_dThresholdIntra = sqrt(3);
    g.m_dThresholdDopantInter = 0;
    g.m_dThresholdDopantIntra = 0;
    g.m_nbrLevels = 1;
    g.m_sphereInfluenceRad = -1;

    /// Segment system and create graph edges///
    g.segmentSystem();
    //g.m_system.writeSystem("partition.pdb");
    g.findNeighbors();
    //g.writeGraph("graph.pdb", "graph.psf",false);

    double sigma(0.1);
    g.gdm(sigma,0, (unsigned int) time(NULL));

    //std::vector<double> fValue = { 1e-3, 0.1 , 0.25 , 0.5 , 1 , 1.5 , 1.75, 2 , 3};
    //std::vector<double> fValue = { 3};
    //std::vector<double> fValue = { 1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1};
    std::vector<double> fValue = { 1e-6};

    //std::vector<double> tValue = {2,3 , 4, 5, 6};

    //std::vector<double> fValue = {1e-6};
    std::vector<double> tValue = {6};



    ofstream of("mobility.txt");

    for( int k(0) ; k < tValue.size() ; k ++)
    {

        for(int i(0) ; i < fValue.size() ; i++)
        {

            /// External parameters for mobility computation ///
            //g.m_temperature = 300;
            g.m_temperature = sigma /(tValue[k]*kB);
            //g.m_concentration = 0.05; // Number of sites occupied
            g.m_concentration = fValue[i];
            Vect3 field(1,0,0);
            //field = field * ( kB * g.m_temperature / 10.0 / g.m_system.m_boxSize);
            //field = field * fValue[i] * ( sigma / (1) );
            //field = field * 0.1 * ( sigma / (1) );
            field = field * 1e-3 * sigma ;
            std::cout << "Field is : " << field.m_x << " " << field.m_y << " " << field.m_z << " V/Angs" << '\n';
            g.m_field = field;

            /// Miller Abrahams parameters ///
            g.m_w0_intra = 1;
            g.m_w0_inter = 1;
            g.m_w0_dopantIntra = 1;
            g.m_w0_dopantInter = 1;
            g.m_loc_inter = 0.1;
            g.m_loc_dopantInter = 0.1;
            g.m_loc_intra = 0.1;
            g.m_loc_dopantIntra = 0.1;


            /// Compute mobility ///
            g.orderSites(0);



            /*for(int u(0) ; u < g.m_orderMeanField.size() ; u++)
            {
                std::cout << g.m_site[g.m_orderMeanField[u]].m_pos.m_x << " " << g.m_site[g.m_orderMeanField[u]].m_pos.m_y << " " << g.m_site[g.m_orderMeanField[u]].m_pos.m_z << '\n';
            }*/



            g.computeRate(maPheno);
            double mob = computeMobility(g, false, "", 0);
            //g.writeCurrent("occup.pdb", "current.pdb", "current.psf",true);
            g.writeGraph("graph.pdb","graph.psf",false);




        of << tValue[k] << "    " << fValue[i] << "    " <<  mob << '\n';


        ofstream occupFlux("occup.txt");
        for (int h(0) ; h < g.m_site.size() ; h++)
        {
            occupFlux << g.m_site[h].m_pos.m_x << " " << g.m_site[h].m_pos.m_y <<  " " << g.m_site[h].m_pos.m_z  << "   " << g.m_site[h].m_occup << '\n';
        }
        occupFlux.close();

        }

    of << '\n';

    }

    of.close();
}


void prepareForGorilla(std::string fileNameTop, std::string fileNamecoord, std::string folder)
{
    /// Extend and convert into gromacs a set of pdb file in a folder ///
    vector<string> coordFile;
    coordFile.push_back(fileNamecoord);


    //Determine number of pedot and tos
    ifstream is(fileNameTop);
    string line;
    int nbPdt(0);
    int nbTos(0);
    int nbAtTos(18);
    int nbAtPedot(0);
    while( line != "[ molecules ]")
    {
        getline(is, line);
    }
    getline(is, line);
    is >> line >> nbPdt >> line >> nbTos;
    is.close();
    int chainLength(round(600.0/nbPdt));
    std::cout << "CHAIN LENGTH !!!!!!!!!!   " << chainLength << '\n';
    nbAtPedot = (chainLength-2)*13+28;

    for (int i(0) ; i < coordFile.size() ; i++)
    {
        string filePdbExtended = coordFile[i].substr(0, coordFile[i].size()-4) + "_big.pdb";
        string filePdbExtendedParsed = coordFile[i].substr(0, coordFile[i].size()-4) + "_big_parsed.pdb";
        extendPdb(coordFile[i]);
        parseBigPdb(filePdbExtended, filePdbExtendedParsed, nbPdt, nbTos, nbAtPedot, nbAtTos);
        string topExtended =  fileNameTop.substr(0, fileNameTop.size()-4) + "_big.top";
        extendTop(fileNameTop);

        string fileNameGro( coordFile[i].substr(0, coordFile[i].size()-4) + ".gro");
        Graph g;
        g.m_fileNameTop = topExtended;
        g.m_fileNameCoord = filePdbExtendedParsed;
        g.m_segmentLength = 1;
        g.initSystem();
        g.segmentSystem();
        g.m_system.writeSystemGro(fileNameGro);

        // Remove all the unused intermediate file
        remove("extendPdb.tcl");
        remove(filePdbExtended.c_str());
        string filePdbExtendedT = filePdbExtended.substr(0, filePdbExtended.size()-4) + "T.pdb";
        remove(filePdbExtendedT.c_str());

        // Displace files in folder
        rename(filePdbExtendedParsed.c_str(), (folder + filePdbExtendedParsed).c_str());
        rename(topExtended.c_str(), (folder + topExtended).c_str());
        rename(fileNameGro.c_str(), (folder + fileNameGro).c_str());

    }
}


void plotPercolation(Graph &g, double threshMin, double threshMax, double threshStep, std::string folder, bool dimCluster)
{

    /// Threshold curve parameters ///
    // All the values in unit log10(H*1000) H in eV
    std::vector<double> pThreshold;
    double thresh(threshMin-threshStep);
    while(thresh < threshMax)
    {
        thresh += threshStep;
        pThreshold.push_back(exp(thresh*log(10))/1000.0);
    }


    /// Compute the threshold behavior ///
    std::ofstream oFlux(folder + "percolationPlot.txt");

    if(oFlux)
    {
        if(dimCluster)
        {
            oFlux << "Threshold" << " " << "clusterMaxNorm" << "   " << "clusterMaxXYZ" << " " << "clusterMaxX" << "   " << "clusterMaxY" << "   " << "clusterMaxZ" << "    " << "clusterMaxNb" << '\n';
        }
        else
        {
            oFlux << "Threshold" << "    " << "clusterMaxNb" << '\n';
        }


        for (int i(0) ; i < pThreshold.size() ; i ++)
        {
            g.cluster(pThreshold[i]);
            if(dimCluster)
            {
                g.percolation();
            }

            double clusterMaxNorm(0);
            double clusterMaxX(0);
            double clusterMaxY(0);
            double clusterMaxZ(0);
            double clusterMaxXYZ(0);
            double clusterMaxNb(0);

            for(int i(0) ; i < g.m_clusterNb ; i ++)
            {
                if(dimCluster)
                {
                    double cSizeXYZ = max(max(g.m_clusterDim[i].m_x,g.m_clusterDim[i].m_y) , g.m_clusterDim[i].m_z);
                    if(cSizeXYZ > clusterMaxXYZ)
                    {
                        clusterMaxXYZ = cSizeXYZ;
                    }
                    double cSizeX = g.m_clusterDim[i].m_x;
                    if(cSizeX > clusterMaxX)
                    {
                        clusterMaxX = cSizeX;
                    }
                    double cSizeY = g.m_clusterDim[i].m_y;
                    if(cSizeY > clusterMaxY)
                    {
                        clusterMaxY = cSizeY;
                    }
                    double cSizeZ = g.m_clusterDim[i].m_z;
                    if(cSizeZ > clusterMaxZ)
                    {
                        clusterMaxZ = cSizeZ;
                    }

                    double cSizeNorm = g.m_clusterDim[i].norm();
                    if(cSizeNorm > clusterMaxNorm)
                    {
                        clusterMaxNorm = cSizeNorm;
                    }
                }

                if(g.m_clusterSiteNb[i] > clusterMaxNb)
                {
                    clusterMaxNb = g.m_clusterSiteNb[i];
                }
            }

            if(dimCluster)
            {
                oFlux << log(pThreshold[i] * 1000.0)/log(10) << " " << clusterMaxNorm << "   " << clusterMaxXYZ << " " << clusterMaxX << "   " << clusterMaxY << "   " << clusterMaxZ << "  " << clusterMaxNb << '\n';
            }
            else
            {
                oFlux << log(pThreshold[i] * 1000.0)/log(10) << " " << clusterMaxNb << '\n';
            }
        }

    oFlux.close();
    }

    else
    {
        cerr << "Can't open file" << '\n';
    }
}


void computeMobilityDependance(Graph &g, std::vector<double> pField, std::vector<double> pTemperature, std::vector<double> pConcentration, std::string folder)
{
    /// Compute the mobility as a function of field///
    std::ofstream oFlux(folder + "Mobility_dependance.txt");
    if(oFlux)
    {

        for (int u(0) ; u < 3 ; u ++)
        {
            string fName;
            if(u == 0) {fName = "X";}
            if(u == 1) {fName = "Y";}
            if(u == 2) {fName = "Z";}

            oFlux << fName << " " << "field(V/cm)" << " " << "temperature (K)" << " " << "concentration (charge/site)" << " " << "Mobility(a.u.)" << '\n';

            for (int i(0) ; i < pField.size() ; i ++)
            {
                for(int j(0) ; j < pTemperature.size() ; j ++)
                {
                    for (int k(0) ; k < pConcentration.size() ; k ++)
                    {
                        // Set up the field
                        Vect3 field(0,0,0);
                        if(u == 0) {field.m_x = 1;}
                        if(u == 1) {field.m_y = 1;}
                        if(u == 2) {field.m_z = 1;}
                        field = field * pField[i] * 1e-8;
                        std::cout << "Field is : " << field.norm() * 1e8 << " V/cm" << '\n' << std::flush;
                        g.m_field = field;

                        // Set up the temperature and concentration
                        g.m_temperature = pTemperature[j];
                        g.m_concentration = pConcentration[k];

                        // Compute rate according to these parameters
                        g.computeRate(maTransfer);
                        //g.computeRate(maPheno);

                        // Compute mobility and save results
                        //bool usePrevDistrib(i!=0);
                        bool usePrevDistrib(false);
                        double mob = computeMobility(g,usePrevDistrib, folder, u);
                        oFlux << field.norm() * 1e8 << " " << g.m_temperature << "  " << g.m_concentration << " " << mob << '\n' ;
                        oFlux.flush();

                        // Save carrier distribution
                        ofstream occupFlux(folder + fName + "_distrib_F" + std::to_string(field.norm() * 1e8) + "_T" + std::to_string(g.m_temperature) + "_P" + std::to_string(g.m_concentration) + ".txt" );
                        for (int h(0) ; h < g.m_site.size() ; h++)
                        {
                            occupFlux << g.m_site[h].m_pos.m_x << " " << g.m_site[h].m_pos.m_y <<  " " << g.m_site[h].m_pos.m_z << "    " << g.m_site[h].m_level << "   " << g.m_site[h].m_occup << '\n';
                        }
                        occupFlux.close();

                    }
                }
            }

        }

        oFlux.close();

    }

    else
    {
        cerr << "Can't open file" << '\n';
    }

}


Graph initGraph(std::string fileNameTop, std::string fileNamecoord, std::string configFile, unsigned int seed, std::string folder, bool readPhysics)
{

    /// Create a graph and init its parameters ///
    Graph g;
    g.m_fileNameCoord = fileNamecoord;
    g.m_fileNameTop = fileNameTop;

    ifstream is(configFile);
    if(is.is_open())
    {
        string line;
        double param;
        getline(is,line);

        /// Graphic settings for VMD site orientation vectors ///
        is >> line >> param;
        g.m_vectLength = param;
        is >> line >> param;
        g.m_vectRadius = param;
        is >> line >> param;
        g.m_vectResol = param;
        getline(is,line);
        getline(is,line);
        getline(is,line);

        /// Nodes segmentation parameters ///
        is >> line >> param;
        g.m_segmentLength = param;
        getline(is,line);
        getline(is,line);
        getline(is,line);

        /// Edges threshold parameters ///
        is >> line >> param;
        g.m_dThresholdInter = param;
        is >> line >> param;
        g.m_dThresholdIntra = param;
        is >> line >> param;
        g.m_dThresholdDopantInter = param;
        is >> line >> param;
        g.m_dThresholdDopantIntra = param;
        is >> line >> param;
        g.m_sphereInfluenceRad = param;
        getline(is,line);
        getline(is,line);
        getline(is,line);

        /// Miller Abrahams parameters ///
        is >> line >> param;
        g.m_w0_intra = param;
        is >> line >> param;
        g.m_w0_inter = param;
        is >> line >> param;
        g.m_w0_dopantIntra = param;
        is >> line >> param;
        g.m_w0_dopantInter = param;
        is >> line >> param;
        g.m_loc_inter = param;
        is >> line >> param;
        g.m_loc_dopantInter = param;
        is >> line >> param;
        g.m_loc_intra = param;
        is >> line >> param;
        g.m_loc_dopantIntra = param;
        getline(is,line);
        getline(is,line);
        getline(is,line);

        /// MO levels information ///
        is >> line >> param;
        g.m_nbrLevels = param;
        is >> line >> param;
        g.m_homoID = param;
        is >> line >> param;
        g.m_nbrLevels2Read = param;
        is >> line >> param;
        g.m_dosBroadening = param;
        is >> line;
        for(int p(0) ; p < g.m_nbrLevels2Read ; p++)
        {
            is >> param;
            g.m_levelEnergies.push_back(param);
        }

    }
    else
    {
        cerr << "Can't open config file !" << '\n';
    }

    /// Input files to read and definition of some output files ///
    string energiesFile("energies.txt");
    vector<string> transferFiles;
    for (int i(0) ; i < g.m_nbrLevels ; i++)
    {
        for (int j(0) ; j < g.m_nbrLevels ; j ++ )
        {
                transferFiles.push_back("transfer_" + std::to_string(i+g.m_homoID) + "_" + std::to_string(j+g.m_homoID) + ".txt");
                std::cout << "transfer_" + std::to_string(i+g.m_homoID) + "_" + std::to_string(j+g.m_homoID) + ".txt" << '\n';
        }
    }
    string rawTransferFile("rawTransfer.txt");
    string rawEnergiesFile("rawEnergies.txt");
    string rawRateFile("rawRate.txt");
    string distAngleTranfer("3DTranfer.txt");

    /// Init the system ///
    g.initSystem();
    g.segmentSystem();
    g.findNeighbors();

    if(readPhysics)
    {

        // Output the some stats on the linking in the graph
        int nbLinkMax(-1e6);
        int nbLinkMin(1e6);
        double nbLinkmean(0);
        for (int i(0) ; i < g.m_site.size() ; i ++)
        {
            int nbLink(g.m_site[i].m_neighborIdx.size());
            if(nbLink > nbLinkMax) {nbLinkMax = nbLink;}
            if(nbLink < nbLinkMin) {nbLinkMin = nbLink;}
            nbLinkmean += nbLink;
        }
        nbLinkmean /= g.m_site.size();
        std::cout << "Max number of links: " << nbLinkMax << "  Min number of links: " << nbLinkMin << "    Mean number of links: " << nbLinkmean << '\n';

        g.importEnergies(energiesFile);
        g.importTransferIntegral(transferFiles);
        g.parseTransfer();

        /// Transfer integral min-max ///
        cout << "Min and Max transfer integral : " << '\n';
        int cTranfer(0);
        for (int x(0) ; x < g.m_nbrLevels ; x ++)
        {
            for(int y(0) ; y < g.m_nbrLevels ; y ++)
            {
                cout << "Level 1:  " << x << "  Level 2:  " << y << "  Min transfer integral:  " << g.m_minTransfer[cTranfer] << "  Max transfer integral:  " << g.m_maxTransfer[cTranfer] << '\n';
                cTranfer ++;
            }
        }


        /// Add molecular levels (and If we want to attribute energy levels from a gaussian distribution) ///
        g.orderSites(0);
        g.addMolLevels2Site();
        if(g.m_dosBroadening != 0)
        {
            g.gdm(g.m_dosBroadening,0, seed);
        }
        g.writeRawEnergies(folder+rawEnergiesFile);
        g.writeRawTranfer(folder+rawTransferFile);
        g.m_temperature = 300;
        g.computeRate(maTransfer);

        /// Check that we have full connectivity in the system ///
        cout << "Check full connectivity... ";
        g.cluster(1e-70);
        g.percolation();
        if(g.m_clusterNb != 1)
        {
            cout << "No full connectivity in the system... aborting" << '\n';
        }
        else
        {
            double xSize(g.m_clusterDim[0].m_x);
            double ySize(g.m_clusterDim[0].m_y);
            double zSize(g.m_clusterDim[0].m_z);

            if(xSize == 1e3 && ySize == 1e3 && zSize == 1e3)
            {
                cout << "Done." << '\n';
            }
            else
            {
                cout << "No full connectivity in the system... aborting" << '\n';
            }
        }

    }

    return g;

}

void readMobilityParam(std::string configFile, std::vector<double> &pField, std::vector<double> &pTemperature, std::vector<double> &pConcentration)
{
    ifstream is(configFile);
    if(is.is_open())
    {
        string line;
        double param;
        int paramN;
        for (int i(0) ; i < 33 ; i ++)
        {
            getline(is,line);
        }
        is >> paramN;
        for(int i(0) ; i < paramN ; i ++)
        {
            is >> param;
            pField.push_back(param);

        }

        getline(is,line);
        getline(is,line);
        getline(is,line);

        is >> paramN;
        for(int i(0) ; i < paramN ; i ++)
        {
            is >> param;
            pTemperature.push_back(param);
        }

        getline(is,line);
        getline(is,line);
        getline(is,line);
        is >> paramN;
        for(int i(0) ; i < paramN ; i ++)
        {
            is >> param;
            pConcentration.push_back(param);
        }


        is.close();

    }
    else
    {
        std::cerr << "Can't open config file!" << '\n';
    }

}



/// Coarse grain applications ///
/// ///////////////////////// ///
/// ///////////////////////// ///

void systemMorpho(std::string fileNameCoord)
{
    Coarsegrain morpho;
    std::cout << "here" << '\n';
    morpho.readSystem(fileNameCoord);
    morpho.initGeom();
    morpho.arrangeCrystalMol();
    morpho.writeFragmentOrient("orient.tcl");
    morpho.writeSystem("morphology.pdb");
    morpho.computeCrystalDist("distanceCryst.txt");
}
