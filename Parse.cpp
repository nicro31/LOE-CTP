#include "Parse.h"

using namespace std;

void parseBigPdb(string fileI, string fileO, int nbPdt, int nbTos, int nbAtPedot, int nbAtTos)
{
    cout << "Parsing the big pdb file... ";

    vector<string> header;
    vector<string> data;

    string line;

    ifstream is(fileI);
    ofstream os(fileO);

    // Read header
    for (int i(0) ; i < 5 ; i++)
    {
        getline(is, line);
        header.push_back(line);
        //cout << line << '\n';
    }

    // Read file
    while(!is.eof())
    {
        getline(is, line);
        if(line.size() >=13)
        {
            if(line.substr(13,2) != "16")
            {
                data.push_back(line);
            }
        }
        //cout << line << '\n';
    }
    data.pop_back();
    data.pop_back();

    is.close();

    /// Write in new file with correct molecules order (pedot and then tos)
    for (int i(0) ; i < header.size() ; i ++)
    {
        os << header[i] << '\n';
    }

    int counter(1);
    string aIdx;
    string pIdx;
    string tIdx;
    int counterRes(0);

    /// Write pedot
    for (int i(0) ; i < 8 ; i ++)
    {
        for (int j(0) ; j < nbAtPedot * nbPdt; j++)
        {
            int idx(i*(nbPdt*nbAtPedot + nbTos*nbAtTos) + j);

            // Correct number of atoms
            aIdx = std::to_string(counter);
            while(aIdx.length() < 5)
            {
                aIdx = " " + aIdx;
            }
            data[idx].replace(6,5,aIdx);

            // Correct number chains
            if((counter-1) % nbAtPedot == 0) {counterRes++;}
            pIdx = std::to_string(counterRes);
            while(pIdx.length() < 4)
            {
                pIdx = " " + pIdx;
            }
            data[idx].replace(22,4,pIdx);

            os << data[idx] << '\n';
            counter ++;
        }
    }


    /// Write tos
    int counterT(0);
    for (int i(0) ; i < 8 ; i ++)
    {
        for (int j(0) ; j < nbAtTos * nbTos; j++)
        {
            int idx(i*(nbPdt*nbAtPedot + nbTos*nbAtTos) + nbAtPedot*nbPdt + j);

            aIdx = std::to_string(counter);
            while(aIdx.length() < 5)
            {
                aIdx = " " + aIdx;
            }
            data[idx].replace(6,5,aIdx);

            // Correct number tos
            if(counterT % nbAtTos == 0) {counterRes++;}
            tIdx = std::to_string(counterRes);
            while(tIdx.length() < 4)
            {
                tIdx = " " + tIdx;
            }
            data[idx].replace(22,4,tIdx);

            os << data[idx] << '\n';
            counter ++;
            counterT ++;
        }
    }

    os << "END" << '\n';

    cout << "Done." << '\n';
}


void extendPdb(string fileI)
{
    cout << "Extending the pbd file with periodic boundary conditions... ";

    /// Create and execute the Tcl script that will perform the extension (Note VMD has to be in the PATH) ///
    ofstream script("extendPdb.tcl");
    string fileIShort = fileI.substr(0, fileI.size()-4);
    string filePdbT(fileIShort + "_bigT.pdb");
    string filePdb(fileIShort + "_big.pdb");
    if (script.is_open())
    {
        script << "# load required package" << '\n' << "package require topotools" << '\n' << "# load a molecule"
        << '\n' << "set mol [mol new " << fileI << " type pdb waitfor all]" << '\n' << "# do the magic"
        << '\n' << "set newmol [::TopoTools::replicatemol $mol 2 2 2 ]" << '\n' << "animate write pdb " << filePdbT
        << " $newmol" << '\n' << "quit" ;

        script.close();

        system("vmd -dispdev text -e extendPdb.tcl");
    }
    else
    {
        cout << "Can't open tcl script file to extend pdb volume." << '\n';
    }


    /// Open the new (extended) pdb file and write the header ///
    vector<string> header;
    string line;
    ifstream is(fileI);
    // Read header
    for (int i(0) ; i < 5 ; i++)
    {
        getline(is, line);
        header.push_back(line);
    }
    is.close();

    ifstream pdbT(filePdbT);
    ofstream pdb(filePdb);

    if (pdb.is_open() && pdbT.is_open())
    {
        pdb << header[0] << '\n' << header[1] << '\n' << header[2] << '\n';
        getline(pdbT, line);
        pdb << line << '\n' << header[4] << '\n';

        while(!pdbT.eof())
        {
            getline(pdbT, line);
            pdb << line << '\n';

        }

        pdb.close();
        pdbT.close();

    }
    else
    {
        cout << "Can't open extended pdb file to write header." << '\n';
    }

    cout << "Done." << '\n';
}


void extendTop(string fileI)
{
    cout << "Writing extended top file... ";

    string topExtended = fileI.substr(0, fileI.size()-4) + "_big.top";

    ifstream top(fileI);
    ofstream topEx(topExtended);
    string line("");

    if (top.is_open() && topEx.is_open())
    {
        while(line != "[ molecules ]" )
        {
            getline(top, line);
            topEx << line << '\n';
        }
        getline(top, line);
        topEx << line << '\n';
        int nbPdt(0);
        int nbTos(0);
        top >> line >> nbPdt >> line >> nbTos;
        topEx << "pedot    " << 8*nbPdt << '\n' << "tos    " << 8*nbTos << '\n';


        top.close();
        topEx.close();

    }
    else
    {
        cout << "Can't open top file to write." << '\n';
    }

    cout << "Done." << '\n';

}
