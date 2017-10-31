#include "System.h"

void System::readSystem(std::string fileNameTop, std::string fileNameCoord)
{
    /// Init topology ///
    std::cout << "Reading topology... ";
    m_topo.readTopo(fileNameTop);
    std::cout << "done" << '\n';
    std::cout << "Creating topology map... ";
    m_topo.createTopoMap();
    std::cout << "done" << '\n';
    std::cout << "Reading system coordinates... ";

    /// Read system ///
    Polymer p;
    RigidFragment monomer;
    monomer.m_type = edot;
    for ( int i(0) ; i < m_topo.m_chainLength ; i ++)
    {
        if( i == 0 || i == m_topo.m_chainLength -1) {monomer.m_atom.resize(14);}
        else { monomer.m_atom.resize(13); }
        p.fragment.push_back(monomer);
    }
    RigidFragment tosylate;
    tosylate.m_type = tos;
    tosylate.m_atom.resize(18);
    std::string line;
    double nbr;
    Atom a;

    std::ifstream inFlux(fileNameCoord.c_str());
    if (inFlux.is_open())
    {
        int counter(0);

        // With 1 supplementary line at beginning
        getline(inFlux,m_line1);

        getline(inFlux,m_line1);
        getline(inFlux,m_line2);
        inFlux >> line;
        inFlux >> nbr;
        m_boxSize.m_x = nbr;
        inFlux >> nbr;
        m_boxSize.m_y = nbr;
        inFlux >> nbr;
        m_boxSize.m_z = nbr;
        getline(inFlux,m_line3);
        getline(inFlux,m_line4);

        for (int i(0) ; i < m_topo.m_pedotNbr ; i ++ )
        {
            for (int j(0) ; j < (m_topo.m_chainLength-2)*13+28 ; j ++ )
            {
                inFlux >> line;
                inFlux >> nbr;
                inFlux >> nbr;
                a.m_name = nbr;
                line.resize(10);
                inFlux.read(&line[0],10);
                inFlux >> nbr;

                inFlux >> nbr;
                a.m_pos.m_x = nbr;
                inFlux >> nbr;
                a.m_pos.m_y = nbr;
                inFlux >> nbr;
                a.m_pos.m_z = nbr;
                inFlux >> nbr;
                inFlux >> nbr;
                a.m_type = m_topo.m_atomPedot[j];
                int id1(m_topo.m_topoMap[j].m_id1);
                int id2(m_topo.m_topoMap[j].m_id2);
                p.fragment[id1].m_atom[id2] = a;
                p.fragment[id1].m_molNbr = counter;

            }

            m_pedot.push_back(p);
            counter ++;

        }

        for (int i(0) ; i < m_topo.m_tosNbr ; i ++ )
        {
            for (int j(0) ; j < 18  ; j ++ )
            {
                inFlux >> line;
                inFlux >> nbr;
                inFlux >> nbr;
                line.resize(10);
                inFlux.read(&line[0],10);
                inFlux >> nbr;

                inFlux >> nbr;
                a.m_pos.m_x = nbr;
                inFlux >> nbr;
                a.m_pos.m_y = nbr;
                inFlux >> nbr;
                a.m_pos.m_z = nbr;
                inFlux >> nbr;
                inFlux >> nbr;
                a.m_type = m_topo.m_atomTos[j];
                tosylate.m_atom[j] = a;
                tosylate.m_molNbr = counter;
            }

            m_tos.push_back(tosylate);
            counter ++;
        }

        inFlux.close();
    }
    else
    {
        std::cerr << "Unable to open coordinates file" << std::endl;
    }
    std::cout << "done" << '\n';
}


void System::writeSystem(std::string fileNameCoord)
{
    std::cout << "Writing the system... ";
    int counter(1);
    int counterG(1);
    Atom a;

    std::ofstream outFlux(fileNameCoord.c_str());
    if (outFlux.is_open())
    {

        outFlux << m_line1;
        outFlux << '\n';
        outFlux << m_line2;
        outFlux << '\n';
        outFlux << "CRYST1";
        outFlux << std::setw(9) << m_boxSize.m_x;
        outFlux << std::setw(9) << m_boxSize.m_y;
        outFlux << std::setw(9) << m_boxSize.m_z;
        outFlux << m_line3;
        outFlux << '\n';
        outFlux << m_line4;
        outFlux << '\n';

        for (int i(0) ; i < m_topo.m_pedotNbr ; i ++ )
        {
            for (int j(0) ; j < (m_topo.m_chainLength-2)*13+28 ; j ++ )
            {

                int id1((m_topo.m_topoMap[j]).m_id1);
                int id2(m_topo.m_topoMap[j].m_id2);
                a = m_pedot[i].fragment[id1].m_atom[id2];

                outFlux << "ATOM";
                outFlux << std::setw(7);
                outFlux << counter;
                outFlux << std::setw(3);
                outFlux << a.m_type;
                outFlux << std::setw(8);
                outFlux << "RES A";
                outFlux << std::setw(4);
                outFlux << counterG;
                outFlux << std::fixed << std::setprecision(3) << std::setw(12) ;

                outFlux << a.m_pos.m_x;
                outFlux << std::setw(8);
                outFlux << a.m_pos.m_y;
                outFlux << std::setw(8);
                outFlux << a.m_pos.m_z;
                outFlux << std::fixed << std::setprecision(2) << std::setw(6);

                //outFlux << (double) id1;
                outFlux << m_pedot[i].fragment[id1].m_fragNbr ;
                outFlux << std::fixed << std::setprecision(2) << std::setw(6);
                outFlux << (double) m_pedot[i].fragment[id1].m_type;
                outFlux << '\n';

                counter ++;

            }

            counterG ++;

        }

        for (int i(0) ; i < m_topo.m_tosNbr ; i ++ )
        {
            for (int j(0) ; j < 18  ; j ++ )
            {
                a = m_tos[i].m_atom[j];

                outFlux << "ATOM";
                outFlux << std::setw(7);
                outFlux << counter;
                outFlux << std::setw(3);
                outFlux << a.m_type;
                outFlux << std::setw(8);
                outFlux << "RES A";
                outFlux << std::setw(4);
                outFlux << counterG;
                outFlux << std::fixed << std::setprecision(3) << std::setw(12) ;

                outFlux << a.m_pos.m_x;
                outFlux << std::setw(8);
                outFlux << a.m_pos.m_y;
                outFlux << std::setw(8);
                outFlux << a.m_pos.m_z;
                outFlux << std::fixed << std::setprecision(2) << std::setw(6);

                outFlux << 0.0;
                outFlux << std::fixed << std::setprecision(2) << std::setw(6);
                outFlux << (double) m_tos[i].m_type;
                outFlux << '\n';

                counter ++;
            }

            counterG ++;
        }

        outFlux << "TER" << '\n' << "ENDMDL" << '\n' ;

        outFlux.close();
    }
    else
    {
        std::cerr << "Unable to open coordinates file" << std::endl;
    }
    std::cout << "done" << '\n';
}


void System::writeSystemGro(std::string fileNameGro)
{
    std::cout << "Writing the system (Gromacs format)... ";
    int counter(1);
    Atom a;

    if(fileNameGro == "")
    {
        fileNameGro = "Pedot.gro";
    }

    std::ofstream outFlux(fileNameGro.c_str());
    if (outFlux.is_open())
    {

        outFlux << "PEDOT";
        outFlux << '\n';
        outFlux << std::setw(5) << m_topo.m_pedotNbr *  ((m_topo.m_chainLength-2)*13+28);
        outFlux << '\n';

        for (int i(0) ; i < m_topo.m_pedotNbr ; i ++ )
        {
            for (int j(0) ; j < (m_topo.m_chainLength-2)*13+28 ; j ++ )
            {

                int id1((m_topo.m_topoMap[j]).m_id1);
                int id2(m_topo.m_topoMap[j].m_id2);
                a = m_pedot[i].fragment[id1].m_atom[id2];

                outFlux << std::setw(5);
                outFlux << i+1;
                outFlux << std::setw(3);
                outFlux << "pdt";
                outFlux << std::setw(7);

                std::string line;
                if(a.m_type == C) {line = "C";}
                if(a.m_type == O) {line = "O";}
                if(a.m_type == H) {line = "H";}
                if(a.m_type == S) {line = "S";}

                outFlux << line;
                outFlux << std::setw(5);
                outFlux << counter;

                outFlux << std::fixed << std::setprecision(4) << std::setw(8) ;
                outFlux << a.m_pos.m_x * 0.1;
                outFlux << std::setw(8);
                outFlux << a.m_pos.m_y * 0.1;
                outFlux << std::setw(8);
                outFlux << a.m_pos.m_z * 0.1;
                /*outFlux << std::fixed << std::setprecision(4) << std::setw(8);
                outFlux << 0.0;
                outFlux << std::setw(8);
                outFlux << 0.0;
                outFlux << std::setw(8);
                outFlux << 0.0;*/

                outFlux << '\n';

                counter ++;

            }

        }

        outFlux << std::fixed << std::setprecision(4) << m_boxSize.m_x * 0.1 << " " << m_boxSize.m_y * 0.1 << " " << m_boxSize.m_z * 0.1 << '\n';

        outFlux.close();
    }
    else
    {
        std::cerr << "Unable to open coordinates file" << std::endl;
    }
    std::cout << "done" << '\n';
}


void System::initGeom()
{
    std::cout << "Initialize rigid segments properties... ";

    for ( int i(0) ; i < m_pedot.size() ; i++ )
    {
        for ( int j(0) ; j < m_pedot[i].fragment.size() ; j ++ )
        {
            bool invert(false);
            if( j % 2 == 0 ) {invert = true;}
            m_pedot[i].fragment[j].computeGeom(invert);
        }
    }

    for (int i(0) ; i < m_topo.m_tosNbr ; i ++ )
    {
        bool invert(false);
        m_tos[i].computeGeom(invert);
    }

    std::cout << "Done" << '\n';
}


int newBoxPbc(int box, const Vect3 &pbc)
{
    int newBoxId(box);

    if(pbc.m_x != 0)
    {
        newBoxId = lookupTabPbc(newBoxId,x);
    }
    if(pbc.m_y != 0)
    {
        newBoxId = lookupTabPbc(newBoxId,y);
    }
    if(pbc.m_z != 0)
    {
        newBoxId = lookupTabPbc(newBoxId,z);
    }

    return newBoxId;
}


int lookupTabPbc(int box, Direction dir)
{
    if(box==0)
    {
        if(dir == x) {return 4;}
        if(dir == y) {return 2;}
        if(dir == z) {return 1;}
    }
    if(box==1)
    {
        if(dir == x) {return 5;}
        if(dir == y) {return 3;}
        if(dir == z) {return 0;}
    }
    if(box==2)
    {
        if(dir == x) {return 6;}
        if(dir == y) {return 0;}
        if(dir == z) {return 3;}
    }
    if(box==3)
    {
        if(dir == x) {return 7;}
        if(dir == y) {return 1;}
        if(dir == z) {return 2;}
    }
    if(box==4)
    {
        if(dir == x) {return 0;}
        if(dir == y) {return 6;}
        if(dir == z) {return 5;}
    }
    if(box==5)
    {
        if(dir == x) {return 1;}
        if(dir == y) {return 7;}
        if(dir == z) {return 4;}
    }
    if(box==6)
    {
        if(dir == x) {return 2;}
        if(dir == y) {return 4;}
        if(dir == z) {return 7;}
    }
    if(box==7)
    {
        if(dir == x) {return 3;}
        if(dir == y) {return 5;}
        if(dir == z) {return 6;}
    }
}

void System::boundDistribution(std::string boundDistrib)
{
    std::cout << "Computing bound distribution... ";

    double dThresh(5);
    double anglePThresh(30);
    double angleNThresh(45);
    double angleThresh(2000);
    std::vector<int> nbBound;
    double dMinGlob(1e10);

    for(int i(0) ; i < m_pedot.size()-1 ; i++)
    {
        for(int j(i+1) ; j < m_pedot.size() ; j++)
        {

            int link(0);

            for(int k(0) ; k < m_pedot[i].fragment.size() ; k++)
            {

                double distanceMin(1e8);
                bool periodic(true);
                int id1;
                int id2;

                for(int l(0) ; l < m_pedot[j].fragment.size() ; l++)
                {
                    double d(Vect3::distancePbc(m_pedot[i].fragment[k].m_center,m_pedot[j].fragment[l].m_center,m_boxSize,periodic));
                    if( d < distanceMin)
                    {
                        distanceMin = d;
                        id1 = k;
                        id2 = l;
                    }
                }



                if(distanceMin < dMinGlob)
                {
                    dMinGlob = distanceMin;
                }

                if( distanceMin < dThresh)
                {
                    Vect3 orientPLoc1 = m_pedot[i].fragment[id1].m_orientP;
                    Vect3 orientPLoc2 = m_pedot[j].fragment[id2].m_orientP;
                    double danglePLoc(acos(orientPLoc1*orientPLoc2)*180.0/pi);
                    if(danglePLoc > 90.0) {danglePLoc = 180.0 - danglePLoc;}

                    Vect3 boundLoc = Vect3::displacementPbc(m_pedot[i].fragment[id1].m_center, m_pedot[j].fragment[id2].m_center,m_boxSize);
                    boundLoc = boundLoc / boundLoc.norm();
                    Vect3 orientNLoc1 = m_pedot[i].fragment[id1].m_orientN;
                    double angleNLoc(acos(orientNLoc1*boundLoc)*180.0/pi);
                    if(angleNLoc > 90.0) {angleNLoc = 180.0 - angleNLoc;}

                    //if(danglePLoc <= anglePThresh && angleNLoc <= angleNThresh)
                    if(danglePLoc * angleNLoc <= angleThresh)
                    {
                        link ++;

                    }

                }


            }

            if(link > 0)
            {
                nbBound.push_back(link);
                bool findMol(false);
                bool is1(false);
                bool is2(false);
                int counter(0);

                while( !findMol && counter < m_crystalMol.size() )
                {
                    for(int p(0) ; p < m_crystalMol[counter].size() ; p++)
                    {
                        if(m_crystalMol[counter][p] == i) {is1 = true;}
                        if(m_crystalMol[counter][p] == j) {is2 = true;}
                    }

                    if(is1 || is2)
                    {
                        findMol = true;
                        if(!is1) {m_crystalMol[counter].push_back(i);}
                        if(!is2) {m_crystalMol[counter].push_back(j);}

                    }

                    counter ++;
                }

                if(!findMol)
                {
                    std::vector<int> cryst;
                    cryst.push_back(i);
                    cryst.push_back(j);
                    m_crystalMol.push_back(cryst);
                }
            }

        }
    }



    m_crystalId.resize(m_pedot.size());
    for(int i(0) ; i < m_crystalId.size() ; i++)
    {
        m_crystalId[i] = -1;
    }

    int counterCryst(0);

    for(int i(0) ; i < m_crystalMol.size() ; i++)
    {
        bool present(false);
        std::vector<int> idCryst;

        for(int j(0) ; j < m_crystalMol[i].size() ; j++)
        {
            int iD(m_crystalMol[i][j]);
            if(m_crystalId[iD] != -1)
            {
                present = true;
                idCryst.push_back(m_crystalId[iD]);
            }
        }

        if(!present)
        {
            for(int j(0) ; j < m_crystalMol[i].size() ; j++)
            {
                int iD(m_crystalMol[i][j]);
                m_crystalId[iD] = i;
            }

            counterCryst ++;
        }
        else
        {
            for(int j(0) ; j < m_crystalMol[i].size() ; j++)
            {
                int iD(m_crystalMol[i][j]);
                m_crystalId[iD] = idCryst[0];
            }

            for(int j(1) ; j < idCryst.size() ; j++)
            {
                for(int k(0) ; k < m_crystalMol[idCryst[j]].size() ; k++)
                {
                    m_crystalId[m_crystalMol[idCryst[j]][k]] = idCryst[0];
                }
            }
        }
    }


    m_crystalMol.resize(0);
    std::vector<int> crystalMolNull;
    m_crystalMol.push_back(crystalMolNull);
    std::vector<int> idCryst;
    idCryst.push_back(-1);

    // Arrange m_crystalMol
    for(int i(0) ; i < m_crystalId.size() ; i++)
    {
        int id(m_crystalId[i]);
        int idxVector(-1);

        if(id != -1)
        {
            for(int j(0) ; j < idCryst.size() ; j++)
            {
                if(id==idCryst[j])
                {
                    idxVector = j;
                    m_crystalMol[j].push_back(i);
                }
            }

            if(idxVector == -1)
            {
                std::vector<int> crystalMolTemp;
                crystalMolTemp.push_back(i);
                m_crystalMol.push_back(crystalMolTemp);
                idCryst.push_back(id);
            }
        }
        else
        {
            m_crystalMol[0].push_back(i);
        }


    }

    // Then arrange m_crystalId to have consecutive numbering
    for(int i(0) ; i < m_crystalMol.size() ; i++)
    {
        for(int j(0) ; j < m_crystalMol[i].size() ; j++)
        {
            m_crystalId[m_crystalMol[i][j]] = i;
        }
    }


    //std::cout << "dMiniGLOB    " << dMinGlob << '\n' ;

    std::ofstream flux(boundDistrib, std::ofstream::app);

    if(flux)
    {
        for(int p(0) ; p < nbBound.size() ; p++)
        {
            flux << nbBound[p] << '\n';
        }
        flux.close();
    }
    else
    {
        std::cout << "Can't open file to write bound distribution !" << '\n';
    }


        std::cout << "Writing the system... ";
    int counter(1);

    std::ofstream outFlux("morphology.pdb");
    if (outFlux.is_open())
    {

        outFlux << "CRYST1";
        outFlux << std::setw(9) << m_boxSize.m_x;
        outFlux << std::setw(9) << m_boxSize.m_y;
        outFlux << std::setw(9) << m_boxSize.m_z;
        outFlux << m_line3;
        outFlux << '\n';

        for (int i(0) ; i < m_topo.m_pedotNbr ; i ++ )
        {
            for (int j(0) ; j < m_topo.m_chainLength ; j ++ )
            {

                outFlux << "ATOM";
                outFlux << std::setw(7);
                outFlux << counter;
                outFlux << std::setw(3);
                outFlux << "V";
                outFlux << std::setw(8);
                outFlux << "PDT X";
                outFlux << std::setw(4);
                outFlux << i + 1;
                outFlux << std::fixed << std::setprecision(3) << std::setw(12) ;

                outFlux << m_pedot[i].fragment[j].m_center.m_x;
                outFlux << std::setw(8);
                outFlux <<  m_pedot[i].fragment[j].m_center.m_y;
                outFlux << std::setw(8);
                outFlux <<  m_pedot[i].fragment[j].m_center.m_z;
                outFlux << std::fixed << std::setprecision(2) << std::setw(6);

                outFlux << m_crystalId[i];
                outFlux << std::fixed << std::setprecision(2) << std::setw(6);
                outFlux << 0.0;
                outFlux << '\n';

                counter ++;

            }

        }

        outFlux << "TER" << '\n' << "ENDMDL" << '\n' ;

        outFlux.close();
    }
    else
    {
        std::cerr << "Unable to open coordinates file" << std::endl;
    }



    std::cout << "Done" << '\n';
}



/// Coarse grain functions ///
/// ////////////////////// ///
/// ////////////////////// ///



void Coarsegrain::readSystem(std::string fileNameCoord)
{
    std::cout << "here" << '\n';

    /// Init topology ///
    m_topo.m_pedotNbr = 500;
    m_topo.m_chainLength = 12;
    m_topo.m_tosNbr = 0;
    std::cout << "Reading coordinates... " ;

    /// Read system ///
    Polymer p;
    p.m_chainLength = m_topo.m_chainLength;
    RigidFragment monomer;
    monomer.m_type = edot;
    for ( int i(0) ; i < m_topo.m_chainLength ; i ++)
    {
        p.fragment.push_back(monomer);
    }

    std::string line;
    double nbr;

    std::ifstream inFlux(fileNameCoord.c_str());
    if (inFlux.is_open())
    {
        int counter(0);

        inFlux >> line;
        inFlux >> nbr;
        m_boxSize.m_x = nbr;
        inFlux >> nbr;
        m_boxSize.m_y = nbr;
        inFlux >> nbr;
        m_boxSize.m_z = nbr;
        getline(inFlux,m_line3);

        for (int i(0) ; i < m_topo.m_pedotNbr ; i ++ )
        {
            for (int j(0) ; j < m_topo.m_chainLength ; j ++ )
            {
                monomer.m_molNbr = i;
                monomer.m_siteNbr = i;
                monomer.m_fragNbr = j;

                inFlux >> line;
                inFlux >> nbr;
                inFlux >> line;
                line.resize(10);
                inFlux.read(&line[0],10);
                inFlux >> nbr;

                inFlux >> nbr;
                monomer.m_center.m_x = nbr;
                inFlux >> nbr;
                monomer.m_center.m_y = nbr;
                inFlux >> nbr;
                monomer.m_center.m_z = nbr;
                inFlux >> nbr;
                inFlux >> nbr;

                p.fragment[j] = monomer;

            }

            m_pedot.push_back(p);
            counter ++;

        }

        inFlux.close();
    }
    else
    {
        std::cerr << "Unable to open coordinates file" << std::endl;
    }
    std::cout << "done" << '\n';
}


void Coarsegrain::writeSystem(std::string fileNameCoord)
{
    std::cout << "Writing the system... ";
    int counter(1);

    std::ofstream outFlux(fileNameCoord.c_str());
    if (outFlux.is_open())
    {

        outFlux << "CRYST1";
        outFlux << std::setw(9) << m_boxSize.m_x;
        outFlux << std::setw(9) << m_boxSize.m_y;
        outFlux << std::setw(9) << m_boxSize.m_z;
        outFlux << m_line3;
        outFlux << '\n';

        for (int i(0) ; i < m_topo.m_pedotNbr ; i ++ )
        {
            for (int j(0) ; j < m_topo.m_chainLength ; j ++ )
            {

                outFlux << "ATOM";
                outFlux << std::setw(7);
                outFlux << counter;
                outFlux << std::setw(3);
                outFlux << "V";
                outFlux << std::setw(8);
                outFlux << "PDT X";
                outFlux << std::setw(4);
                outFlux << i + 1;
                outFlux << std::fixed << std::setprecision(3) << std::setw(12) ;

                outFlux << m_pedot[i].fragment[j].m_center.m_x;
                outFlux << std::setw(8);
                outFlux <<  m_pedot[i].fragment[j].m_center.m_y;
                outFlux << std::setw(8);
                outFlux <<  m_pedot[i].fragment[j].m_center.m_z;
                outFlux << std::fixed << std::setprecision(2) << std::setw(6);

                outFlux << m_crystalId[i];
                outFlux << std::fixed << std::setprecision(2) << std::setw(6);
                outFlux << 0.0;
                outFlux << '\n';

                counter ++;

            }

        }

        outFlux << "TER" << '\n' << "ENDMDL" << '\n' ;

        outFlux.close();
    }
    else
    {
        std::cerr << "Unable to open coordinates file" << std::endl;
    }
    std::cout << "done" << '\n';





    // Write psf file to link rigid fragment of the same polymer
    std::string nameTmp(fileNameCoord);
    nameTmp.pop_back();
    nameTmp.pop_back();
    nameTmp.pop_back();
    nameTmp.pop_back();
    nameTmp = nameTmp + ".psf";
    int nbEDOT(m_topo.m_pedotNbr*m_topo.m_chainLength);

    std::ofstream outFlux2(nameTmp.c_str());
    if (outFlux2.is_open())
    {
        outFlux2 << "PSF" << '\n' << '\n' ;
        outFlux2 << "       1 !NTITLE" << '\n';
        outFlux2 << " REMARKS connectivity of the graph" << '\n' << '\n';
        outFlux2 << std::setw(8);
        outFlux2 << nbEDOT;
        outFlux2 << std::setw(7);
        outFlux2 << "!NATOM" << '\n' ;


        for (int i(0) ; i < nbEDOT ; i ++ )
        {
            outFlux2 << std::setw(8);
            outFlux2 << i+1;
            outFlux2 << std::setw(4);
            outFlux2 << "RES";
            outFlux2 << std::setw(3);
            int resId((int) (floor((float)i / (float) m_topo.m_chainLength) + 1));
            outFlux2 << 1;
            outFlux2 << std::setw(7);
            outFlux2 << resId;
            outFlux2 << std::setw(5);
            outFlux2 << "RES";
            outFlux2 << std::setw(5);
            outFlux2 << "RES";
            outFlux2 << std::setw(12) << std::fixed << std::setprecision(6);
            outFlux2 << 1.0;
            outFlux2 << std::setw(14) << std::fixed << std::setprecision(4);
            outFlux2 << 1.0;
            outFlux2 << std::setw(12);
            outFlux2 << 0;
            outFlux2 << '\n';

        }

        outFlux2 << '\n';
        outFlux2 << std::setw(8);
        int nbEdges(m_topo.m_pedotNbr*(m_topo.m_chainLength-1));
        outFlux2 << nbEdges;
        outFlux2 << std::setw(14);
        outFlux2 << "!NBOND: bonds" << '\n' ;


        int counter(0);
        int counterGlob(1);
        for (int i(0) ; i < m_topo.m_pedotNbr ; i ++ )
        {
            for (int j(0) ; j < m_topo.m_chainLength-1 ; j ++ )
            {
                    outFlux2 << std::setw(8);
                    outFlux2 << counterGlob;
                    outFlux2 << std::setw(8);
                    outFlux2 << counterGlob+1;
                    counter ++;

                    if (counter % 4 == 0) { outFlux2 << '\n'; }
                    counterGlob++;
                    if(j == m_topo.m_chainLength-2) {counterGlob++;}
            }
        }
        outFlux2.close();
    }
    else
    {
        std::cerr << "Unable to open psf file" << std::endl;
    }


}


void Coarsegrain::initGeom()
{
    std::cout << "Computing structural properties... ";
    double dThresh(6);
    double angleThresh(90);
    int linkThresh(5);
    m_crystalMol.resize(0);

    double dMinGlob(1e10);

    for ( int i(0) ; i < m_pedot.size() ; i++ )
    {
        for ( int j(0) ; j < m_pedot[i].fragment.size() ; j ++ )
        {
            Vect3 orientP;
            if( j != m_pedot[i].fragment.size() - 1)
            {
                orientP = Vect3::displacementPbc(m_pedot[i].fragment[j].m_center , m_pedot[i].fragment[j+1].m_center , m_boxSize);
                orientP = orientP / orientP.norm();
            }
            else
            {
                orientP = m_pedot[i].fragment[j-1].m_orientP;
            }
            m_pedot[i].fragment[j].m_orientP = orientP;
        }
    }


    for(int i(0) ; i < m_pedot.size()-1 ; i++)
    {
        for(int j(i+1) ; j < m_pedot.size() ; j++)
        {

            int link(0);

            for(int k(0) ; k < m_pedot[i].fragment.size() ; k++)
            {

                double distanceMin(1e8);
                bool periodic(true);
                int id1;
                int id2;

                for(int l(0) ; l < m_pedot[j].fragment.size() ; l++)
                {
                    double d(Vect3::distancePbc(m_pedot[i].fragment[k].m_center,m_pedot[j].fragment[l].m_center,m_boxSize,periodic));
                    if( d < distanceMin)
                    {
                        distanceMin = d;
                        id1 = k;
                        id2 = l;
                    }
                }



                if(distanceMin < dMinGlob)
                {
                    dMinGlob = distanceMin;
                }

                if( distanceMin < dThresh)
                {
                    Vect3 orientPLoc1 = m_pedot[i].fragment[id1].m_orientP;
                    Vect3 orientPLoc2 = m_pedot[j].fragment[id2].m_orientP;
                    double danglePLoc(acos(orientPLoc1*orientPLoc2)*180.0/pi);
                    if(danglePLoc > 90.0) {danglePLoc = 180.0 - danglePLoc;}

                    if(danglePLoc <= angleThresh)
                    {
                        link ++;
                    }
                }
            }

            if(link >= linkThresh)
            {
                bool findMol(false);
                bool is1(false);
                bool is2(false);
                int counter(0);

                while( !findMol && counter < m_crystalMol.size() )
                {
                    for(int p(0) ; p < m_crystalMol[counter].size() ; p++)
                    {
                        if(m_crystalMol[counter][p] == i) {is1 = true;}
                        if(m_crystalMol[counter][p] == j) {is2 = true;}
                    }

                    if(is1 || is2)
                    {
                        findMol = true;
                        if(!is1) {m_crystalMol[counter].push_back(i);}
                        if(!is2) {m_crystalMol[counter].push_back(j);}

                    }

                    counter ++;
                }

                if(!findMol)
                {
                    std::vector<int> cryst;
                    cryst.push_back(i);
                    cryst.push_back(j);
                    m_crystalMol.push_back(cryst);
                }

            }

        }
    }

    m_crystalId.resize(m_pedot.size());
    for(int i(0) ; i < m_crystalId.size() ; i++)
    {
        m_crystalId[i] = -1;
    }

    int counterCryst(0);

    for(int i(0) ; i < m_crystalMol.size() ; i++)
    {
        bool present(false);
        std::vector<int> idCryst;

        for(int j(0) ; j < m_crystalMol[i].size() ; j++)
        {
            int iD(m_crystalMol[i][j]);
            if(m_crystalId[iD] != -1)
            {
                present = true;
                idCryst.push_back(m_crystalId[iD]);
            }
        }

        if(!present)
        {
            for(int j(0) ; j < m_crystalMol[i].size() ; j++)
            {
                int iD(m_crystalMol[i][j]);
                m_crystalId[iD] = i;
            }

            counterCryst ++;
        }
        else
        {
            for(int j(0) ; j < m_crystalMol[i].size() ; j++)
            {
                int iD(m_crystalMol[i][j]);
                m_crystalId[iD] = idCryst[0];
            }

            for(int j(1) ; j < idCryst.size() ; j++)
            {
                for(int k(0) ; k < m_crystalMol[idCryst[j]].size() ; k++)
                {
                    m_crystalId[m_crystalMol[idCryst[j]][k]] = idCryst[0];
                }
            }
        }
    }

    //std::cout << "dmin : " << dMinGlob << '\n';

    std::cout << "Done" << '\n';
}


void Coarsegrain::writeFragmentOrient(std::string fileVector)
{
    std::ofstream of(fileVector.c_str());
    double m_vectLength(1.5);
    double m_vectRadius(0.2);
    double m_vectResol(12);

    if(of)
    {
        for(int i(0) ; i < m_pedot.size() ; i++)
        {
            for (int j(0) ; j < m_pedot[i].fragment.size() ; j++ )
            {
                Vect3 v1(m_pedot[i].fragment[j].m_center);
                Vect3 v2(m_pedot[i].fragment[j].m_orientN);
                Vect3 v3(m_pedot[i].fragment[j].m_orientP);

                /*
                // Write orient N
                of << "graphics top color 0" << '\n';
                of << "graphics top cylinder {" << v1.m_x << " " << v1.m_y << " " << v1.m_z << "} {" <<
                v1.m_x + v2.m_x*m_vectLength << " " << v1.m_y + v2.m_y*m_vectLength << " " << v1.m_z + v2.m_z*m_vectLength << "} radius " <<
                m_vectRadius << " resolution " << m_vectResol << " filled yes" << '\n';

                of << "graphics top color 0" << '\n';
                of << "graphics top cone {" <<
                v1.m_x + v2.m_x*m_vectLength << " " << v1.m_y + v2.m_y*m_vectLength << " " << v1.m_z + v2.m_z*m_vectLength << "} {" <<
                v1.m_x + v2.m_x*m_vectLength*1.20 << " " << v1.m_y + v2.m_y*m_vectLength*1.20 << " " << v1.m_z + v2.m_z*m_vectLength*1.20 << "} radius " <<
                2*m_vectRadius << " resolution " << m_vectResol << '\n';
                */

                // Write orient P
                of << "graphics top color 1" << '\n';
                of << "graphics top cylinder {" << v1.m_x << " " << v1.m_y << " " << v1.m_z << "} {" <<
                v1.m_x + v3.m_x*m_vectLength << " " << v1.m_y + v3.m_y*m_vectLength << " " << v1.m_z + v3.m_z*m_vectLength << "} radius " <<
                m_vectRadius << " resolution " << m_vectResol << " filled yes" << '\n';

                of << "graphics top color 1" << '\n';
                of << "graphics top cone {" <<
                v1.m_x + v3.m_x*m_vectLength << " " << v1.m_y + v3.m_y*m_vectLength << " " << v1.m_z + v3.m_z*m_vectLength << "} {" <<
                v1.m_x + v3.m_x*m_vectLength*1.20 << " " << v1.m_y + v3.m_y*m_vectLength*1.20 << " " << v1.m_z + v3.m_z*m_vectLength*1.20 << "} radius " <<
                2*m_vectRadius << " resolution " << m_vectResol << '\n';


            }
        }

        of.close();

    }
    else
    {
        std::cerr << "Can't open script file for writing rigid fragment orientations" << '\n';
    }
}


void Coarsegrain::arrangeCrystalMol()
{
    m_crystalMol.resize(0);
    std::vector<int> crystalMolNull;
    m_crystalMol.push_back(crystalMolNull);
    std::vector<int> idCryst;
    idCryst.push_back(-1);

    // Arrange m_crystalMol
    for(int i(0) ; i < m_crystalId.size() ; i++)
    {
        int id(m_crystalId[i]);
        int idxVector(-1);

        if(id != -1)
        {
            for(int j(0) ; j < idCryst.size() ; j++)
            {
                if(id==idCryst[j])
                {
                    idxVector = j;
                    m_crystalMol[j].push_back(i);
                }
            }

            if(idxVector == -1)
            {
                std::vector<int> crystalMolTemp;
                crystalMolTemp.push_back(i);
                m_crystalMol.push_back(crystalMolTemp);
                idCryst.push_back(id);
            }
        }
        else
        {
            m_crystalMol[0].push_back(i);
        }


    }

    // Then arrange m_crystalId to have consecutive numbering
    for(int i(0) ; i < m_crystalMol.size() ; i++)
    {
        for(int j(0) ; j < m_crystalMol[i].size() ; j++)
        {
            m_crystalId[m_crystalMol[i][j]] = i;
        }
    }
}


void Coarsegrain::computeCrystalDist(std::string fileNameDist)
{
    std::ofstream flux(fileNameDist);

    if(flux)
    {
        for(int i(1) ; i < m_crystalMol.size() - 1 ; i++)
        {
            for(int j(i+1) ; j < m_crystalMol.size() ; j++)
            {
                double distMean(0);
                int counter(0);
                int counter1(0);
                int counter2(0);
                bool periodic(true);
                Vect3 c1;
                Vect3 c2;

                for(int k(0) ; k < m_crystalMol[i].size() ; k++)
                {
                    for(int l(0) ; l < m_crystalMol[j].size() ; l++)
                    {


                        /*for(int m(0) ; m < m_pedot[m_crystalMol[i][k]].fragment.size() ; m++)
                        {
                            for(int n(0) ; n < m_pedot[m_crystalMol[j][l]].fragment.size() ; n++)
                            {

                                distMean += Vect3::distancePbc(m_pedot[m_crystalMol[i][k]].fragment[m].m_center , m_pedot[m_crystalMol[j][l]].fragment[n].m_center , m_boxSize , periodic );
                                counter ++;
                            }
                        }*/


                        for(int m(0) ; m < m_pedot[m_crystalMol[i][k]].fragment.size() ; m++)
                        {
                            c1 += m_pedot[m_crystalMol[i][k]].fragment[m].m_center;
                            counter1 ++;
                        }
                        for(int n(0) ; n < m_pedot[m_crystalMol[j][l]].fragment.size() ; n++)
                        {
                            c2 += m_pedot[m_crystalMol[j][l]].fragment[n].m_center;
                            counter2 ++;
                        }

                    }
                }

                c1 = c1 / counter1;
                c2 = c2 / counter2;
                distMean = Vect3::distancePbc(c1,c2,m_boxSize,periodic);

                //distMean /= counter;

                flux << distMean << '\n';
            }
        }

        flux.close();
    }
    else
    {
        std::cout << "Can't open file !" << '\n';
    }
}
