#include "MeanField.h"

#ifdef _WIN32 || _WIN64
#include <direct.h>
#endif // _WIN32

#ifdef linux
#include <sys/stat.h>
#endif // linux

using namespace std;


int main(int argc, char *argv[])
{
    unsigned int seed(time(NULL));

    /// Parse the input arguments ///
    std::string mode(argv[1]);
    if( (argc < 5 && (mode != "prep" && mode != "coarse") ) || (argc < 4 && mode == "prep") || (argc < 4 && mode == "coarse") )
    {
        cerr << "Not enough arguments !!" << '\n';
        return 0;
    }


    if( ! (mode == "prep" || mode == "perco" || mode == "mob" || mode == "export" || mode == "check" || mode == "network" || mode == "perco_plot" || mode == "system_cluster" || mode == "transfer_geom" || mode == "perco_plot_N" || mode == "perco_N" || mode == "coarse" || mode == "bound_distrib" ) )
    {
        cerr << "Mode not available !!" << '\n';
        return 0;
    }
    std::string folder;
    if(argc == 5 || (mode == "prep" && argc == 4) || (mode == "coarse") ) {folder = "";}
    else
    {
        folder = argv[5];
        folder = folder + "_" + std::to_string(seed);
        #ifdef _WIN32 || _WIN64
        _mkdir((std::string(_getcwd(NULL,0)) + "\\"+folder).c_str());
        #endif // _WIN32
        #ifdef linux
        mkdir(folder.c_str(),0777);
        #endif // linux
        folder = folder + "/";

    }

    /// Get file names
    std::string fileNameTop;
    std::string fileNameCoord;
    if(mode != "coarse")
    {
        fileNameTop = argv[2];
        fileNameCoord = argv[3];
    }
    else
    {
        fileNameCoord = argv[2];
    }
    std::string configFile;
    if( (mode != "prep") && (mode != "coarse") ) {configFile = argv[4];}
    if( mode == "coarse" ) {configFile = argv[3];}

    /// Prepare file for gorilla ///
    if(mode == "prep")
    {
        prepareForGorilla(fileNameTop, fileNameCoord, folder);
    }

    /// Write bound distribution ///
    if(mode == "bound_distrib")
    {
        Graph g = initGraph(fileNameTop,fileNameCoord,configFile, seed, folder, false);
        g.writeboundDistribution("boundDistrib.txt");
    }


    /// Find percolation threshold ///
    if(mode == "perco")
    {
        Graph g = initGraph(fileNameTop,fileNameCoord,configFile, seed, folder);
        Vect3 thresholdPerco;
        thresholdPerco = findPercolationThreshold(g, folder);
    }

    /// Find percolation threshold ///
    if(mode == "perco_N")
    {
        Graph g = initGraph(fileNameTop,fileNameCoord,configFile, seed, folder);
        double thresholdPerco;
        thresholdPerco = findPercolationThresholdN(g, folder);
    }

    /// Export the percolation behavior curve ///
    if(mode == "perco_plot")
    {
        Graph g = initGraph(fileNameTop,fileNameCoord,configFile, seed, folder);
        double threshMin(-6);
        double threshMax(2);
        double threshStep(0.1);
        plotPercolation(g, threshMin, threshMax, threshStep, folder, true);
    }

    /// Export the percolation behavior curve ///
    if(mode == "perco_plot_N")
    {
        Graph g = initGraph(fileNameTop,fileNameCoord,configFile, seed, folder);
        double threshMin(-15);
        double threshMax(3);
        double threshStep(0.1);
        plotPercolation(g, threshMin, threshMax, threshStep, folder, false);
    }


    /// Compute mobility ///
    if(mode == "mob")
    {
        Graph g = initGraph(fileNameTop,fileNameCoord,configFile, seed, folder);
        std::vector<double> pField;
        std::vector<double> pTemperature;
        std::vector<double> pConcentration;
        readMobilityParam(configFile,pField,pTemperature,pConcentration);
        // Transform global carrier concentration in nb carrier per site
        for(int i(0) ; i < pConcentration.size() ; i++)
        {
            pConcentration[i] = pConcentration[i] / (double) g.m_site.size();
        }
        computeMobilityDependance(g, pField, pTemperature, pConcentration ,folder, maTransfer);
    }


    /// Export diagonal and off-diagonal disorder ///
    if(mode == "export")
    {
        Graph g = initGraph(fileNameTop,fileNameCoord,configFile, seed, folder);

    }


    /// Check system integrity ///
    if(mode == "check")
    {
        Graph g = initGraph(fileNameTop,fileNameCoord,configFile, seed, folder);

        cout << "Checking integrity..." << '\n';

        for(int i(0) ; i < g.m_site.size() ; i ++)
        {
            for(int j(0) ; j < g.m_site[i].m_neighborIdx.size() ; j++)
            {
                int idxNeigh(g.m_site[i].m_neighborIdx[j]);
                int id(-1);

                for(int k(0) ; k < g.m_site[idxNeigh].m_neighborIdx.size() ; k++)
                {
                    if(g.m_site[idxNeigh].m_neighborIdx[k] == i)
                    {
                        id = k;
                    }
                }

                assert(id != -1);
                if(id == -1)
                {
                    cout << "Integrity issue in the data !"  << '\n';
                    cout << "Missing neighbor"  << '\n';
                }
                else
                {
                    if( ((g.m_site[i].m_hopTo[j] - g.m_site[idxNeigh].m_hopFrom[id])/g.m_site[i].m_hopTo[j] > 1e-10) || ((g.m_site[i].m_hopFrom[j] -  g.m_site[idxNeigh].m_hopTo[id])/g.m_site[i].m_hopFrom[j] > 1e-10) )
                    {
                        cout << setprecision(20);
                        cout << "Integrity issue in the data !"  << '\n';
                        cout << "Molecules " << i << " and " << idxNeigh << '\n';
                        cout << g.m_site[i].m_hopTo[j]  << '\n';
                        cout << g.m_site[idxNeigh].m_hopFrom[id]  << '\n';
                        cout << g.m_site[i].m_hopFrom[j]  << '\n';
                        cout << g.m_site[idxNeigh].m_hopTo[id]  << '\n';
                        cout << '\n' << '\n';
                    }

                }
            }
        }

        cout << "Done." << '\n';



    }


    /// Export the transport network to visualize with VMD ///
    if(mode == "network")
    {
        string graphPdb(folder+"graph.pdb");
        string graphPsf(folder+"graph.psf");
        Graph g = initGraph(fileNameTop,fileNameCoord,configFile, seed, folder);

        double threshold(0);
        cout << "Value for clustering ? (666 is for no clustering) " << '\n';
        cin >> threshold;

        if(threshold != 666)
        {
            g.cluster(exp(threshold*log(10))/1000.0);
        }

        g.writeGraph(graphPdb,graphPsf,false);

        g.writeCurrent("occup.pdb","current.pdb","current.psf",false);
    }

    /// Export the pedot system to visualize with VMD, with molecule associated to clusterId ///
    if(mode == "system_cluster")
    {
        Graph g = initGraph(fileNameTop,fileNameCoord,configFile, seed, folder);

        double threshold(0);
        cout << "Value for clustering ? (666 is for no clustering) " << '\n';
        cin >> threshold;

        if(threshold != 666)
        {
            g.cluster(exp(threshold*log(10))/1000.0);
        }

        g.writeSystemClustered("system_cluster.pdb",folder);
    }

    /// Export transfer integrals and associated geometries ///
    if(mode == "transfer_geom")
    {
        Graph g = initGraph(fileNameTop,fileNameCoord,configFile, seed, folder);
        string transferGeomFile(folder+"transferGeom.txt");
        g.writeTransferDistanceAngle(transferGeomFile,0,0, maTransfer);
        g.writeFragmentVector("orientation.tcl");
    }

    /// Calculate morphological properties of a coarse grained MD ///
    if(mode == "coarse")
    {
        Graph g = initCoarseGraph(fileNameCoord, configFile);
        string transferGeomFile("transferGeom.txt");
        g.writeTransferDistanceAngle(transferGeomFile,0,0, maPheno);
        double threshMin(-15);
        double threshMax(3);
        double threshStep(0.05);
        plotPercolation(g, threshMin, threshMax, threshStep, "", false);

        /// Init site energy ///
        g.orderSites(0);
        if(g.m_dosBroadening != 0)
        {
            g.gdm(g.m_dosBroadening,0, seed, false);
        }

        std::vector<double> pField;
        std::vector<double> pTemperature;
        std::vector<double> pConcentration;
        readMobilityParam(configFile,pField,pTemperature,pConcentration);
        // Transform global carrier concentration in nb carrier per site
        for(int i(0) ; i < pConcentration.size() ; i++)
        {
            pConcentration[i] = pConcentration[i] / (double) g.m_site.size();
        }
        computeMobilityDependance(g, pField, pTemperature, pConcentration ,folder, maPheno);


    }

    return 0;

}


