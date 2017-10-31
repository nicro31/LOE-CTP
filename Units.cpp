#include "Units.h"

AtomType mass2Type(float mass)
{
    if( mass == (float) 12.0100 ) { return C; }
    if( mass == (float) 16.0000 ) { return O; }
    if( mass == (float) 32.0600 ) { return S; }
    if( mass == (float) 1.0080 ) { return H; }
}

void RigidFragment::computeGeom(bool invert)
{
    if (m_type == edot)
    {
        Vect3 center;
        Vect3 orientN;
        Vect3 orientP;
        Vect3 orientP2;
        Vect3 orientT;

        center = (m_atom[0].m_pos + m_atom[1].m_pos + m_atom[2].m_pos + m_atom[3].m_pos + m_atom[4].m_pos) / 5.0;

        orientP = (m_atom[2].m_pos - m_atom[1].m_pos);
        orientP = orientP / orientP.norm();

        orientP2 = (m_atom[2].m_pos - m_atom[0].m_pos);
        orientP2 = orientP2 / orientP2.norm();

        orientN = orientP^orientP2;
        orientN = orientN / orientN.norm();

        orientT = m_atom[0].m_pos - center;
        orientT = orientT / orientT.norm();

        m_center = center;
        m_orientN = orientN;
        if(invert) {m_orientN = Vect3(0,0,0)-orientN;}
        m_orientP = orientP;
        m_orientT = orientT;


    }
    else
    {
        Vect3 center;
        Vect3 orientN;
        Vect3 orientP;
        Vect3 orientP2;

        center = (m_atom[4].m_pos + m_atom[5].m_pos + m_atom[6].m_pos + m_atom[7].m_pos + m_atom[8].m_pos + m_atom[9].m_pos) / 6.0;

        orientP = (m_atom[4].m_pos - m_atom[7].m_pos);
        orientP = orientP / orientP.norm();

        orientP2 = (m_atom[6].m_pos - m_atom[7].m_pos);
        orientP2 = orientP2 / orientP2.norm();

        orientN = orientP^orientP2;
        orientN = orientN / orientN.norm();

        m_center = center;
        m_orientN = orientN;
        m_orientP = orientP;

    }
}

void Topo::readTopo(std::string fileName)
{
    std::string line;
    double nbr;

    std::ifstream inFlux(fileName.c_str());
    if (inFlux.is_open())
    {
        /// Go to Pedot info ///
        while ( line != "[ atoms ]" )
        {
            getline(inFlux,line);
        }
        getline(inFlux,line);

        /// Reading Pedot chain atoms ///
        inFlux >> line;
        while( line != "[" )
        {
            inFlux >> line;
            inFlux >> nbr;
            inFlux >> line;

            inFlux >> nbr;
            inFlux >> nbr;
            inFlux >> nbr;
            inFlux >> nbr;
            inFlux >> line;

            m_atomPedot.push_back(mass2Type(nbr));
        }
        m_chainLength = (m_atomPedot.size() - 28) / 13 + 2;

        /// Reading Pedot chain bounds ///
        m_neighborPedot.resize(m_atomPedot.size());
        getline(inFlux,line);
        getline(inFlux,line);
        inFlux >> line;
        Bound b;
        while( line != "[" )
        {
            nbr = std::stoi(line.c_str());
            b.m_id1 = nbr - 1;
            inFlux >> nbr;
            b.m_id2 = nbr - 1;
            inFlux >> nbr;
            inFlux >> line;
            m_boundPedot.push_back(b);

            m_neighborPedot[b.m_id1].push_back(b.m_id2);
            m_neighborPedot[b.m_id2].push_back(b.m_id1);

        }

        /// Reading Pedot extremities ///
        getline(inFlux,line);
        getline(inFlux,line);
        inFlux >> line;
        Vect3 angle;
        while( line != "[" )
        {
            int cC(0);
            int cH(0);
            int cS(0);

            nbr = std::stoi(line.c_str());
            angle.m_x = nbr - 1;
            inFlux >> nbr;
            angle.m_y = nbr - 1;
            inFlux >> nbr;
            angle.m_z = nbr - 1;
            inFlux >> nbr;
            inFlux >> line;

            cC = ( (AtomType) m_atomPedot[angle.m_x] == C ) + ( (AtomType) m_atomPedot[angle.m_y] == C ) + ( (AtomType) m_atomPedot[angle.m_z] == C );
            cH = ( (AtomType) m_atomPedot[angle.m_x] == H ) + ( (AtomType) m_atomPedot[angle.m_y] == H ) + ( (AtomType) m_atomPedot[angle.m_z] == H );
            cS = ( (AtomType) m_atomPedot[angle.m_x] == S ) + ( (AtomType) m_atomPedot[angle.m_y] == S ) + ( (AtomType) m_atomPedot[angle.m_z] == S );

            if (cC == 1 && cH == 1 && cS == 1)
            {
                if( (AtomType) m_atomPedot[angle.m_x] == H ) {m_HTerm.push_back(angle.m_x);}
                if( (AtomType) m_atomPedot[angle.m_y] == H ) {m_HTerm.push_back(angle.m_y);}
                if( (AtomType) m_atomPedot[angle.m_z] == H ) {m_HTerm.push_back(angle.m_z);}
            }

        }



        /// Go to Tos info ///
        while ( line != "[ atoms ]" )
        {
            getline(inFlux,line);
        }
        getline(inFlux,line);

        /// Reading Tos chain atoms ///
        inFlux >> line;
        while( line != "[" )
        {
            inFlux >> line;
            inFlux >> nbr;
            inFlux >> line;
            inFlux >> nbr;
            inFlux >> nbr;
            inFlux >> nbr;
            inFlux >> nbr;
            inFlux >> line;

            m_atomTos.push_back(mass2Type(nbr));
        }

        /// Reading Tos chain bounds ///
        m_neighborTos.resize(m_atomTos.size());
        getline(inFlux,line);
        getline(inFlux,line);
        inFlux >> line;
        while( line != "[" )
        {
            nbr = std::stoi(line.c_str());
            b.m_id1 = nbr - 1;
            inFlux >> nbr;
            b.m_id2 = nbr - 1;
            inFlux >> nbr;
            inFlux >> line;
            m_boundTos.push_back(b);

            m_neighborTos[b.m_id1].push_back(b.m_id2);
            m_neighborTos[b.m_id2].push_back(b.m_id1);

        }

        /// Go to molecules info ///
        while ( line != "[ molecules ]" )
        {
            getline(inFlux,line);
        }
        getline(inFlux,line);
        inFlux >> line;
        inFlux >> nbr;
        m_pedotNbr = nbr;
        inFlux >> line;
        inFlux >> nbr;
        m_tosNbr = nbr;

        //std::cout << m_pedotNbr << "    " << m_tosNbr << '\n';

        inFlux.close();
    }
    else
    {
        std::cerr << "Unable to open topology file" << std::endl;
    }


}

void Topo::createTopoMap()
{
    Bound b;
    int id;
    int idC1;
    int idC2;
    int idC3;
    int idC4;
    int idC7;
    int idC8;
    int idO5;
    int idO6;
    int idH9;
    int idH10;
    int idH11;
    int idH12;
    int idS;
    m_topoMap.resize(m_atomPedot.size());

    ///// Pedot chain /////

    /// Add first H terminal atom of a Pedot chain ///
    id = m_HTerm[0];
    b.m_id1 = 0;
    b.m_id2 = 13;
    m_topoMap[id] = b;

    /// Find C number 1 ///
    idC1 = m_neighborPedot[m_HTerm[0]][0];

    for ( int fragNbr(0) ; fragNbr < m_chainLength ; fragNbr ++ )
    {

    /// Add C number 1 ///
    b.m_id1 = fragNbr;
    b.m_id2 = 1;
    m_topoMap[idC1] = b;

    /// Add S number 0 and C number 3 ///
    for (int i(0) ; i < m_neighborPedot[idC1].size() ; i ++)
    {
        if( m_atomPedot[m_neighborPedot[idC1][i]] == S ) { idS = m_neighborPedot[idC1][i];}
        if( m_atomPedot[m_neighborPedot[idC1][i]] == C && m_neighborPedot[idC1][i] != id ) { idC3 = m_neighborPedot[idC1][i];}
    }
    b.m_id1 = fragNbr;
    b.m_id2 = 0;
    m_topoMap[idS] = b;
    b.m_id1 = fragNbr;
    b.m_id2 = 3;
    m_topoMap[idC3] = b;

    /// Add C number 2 ///
    for (int i(0) ; i < m_neighborPedot[idS].size() ; i ++)
    {
        if( m_neighborPedot[idS][i] != idC1 ) { idC2 = m_neighborPedot[idS][i];}
    }
    b.m_id1 = fragNbr;
    b.m_id2 = 2;
    m_topoMap[idC2] = b;

    /// Add O number 5 ///
    for (int i(0) ; i < m_neighborPedot[idC3].size() ; i ++)
    {
        if( m_atomPedot[m_neighborPedot[idC3][i]] == O ) { idO5 = m_neighborPedot[idC3][i];}
    }
    b.m_id1 = fragNbr;
    b.m_id2 = 5;
    m_topoMap[idO5] = b;

    /// Add C number 7 ///
    for (int i(0) ; i < m_neighborPedot[idO5].size() ; i ++)
    {
        if( m_neighborPedot[idO5][i] != idC3 ) { idC7 = m_neighborPedot[idO5][i];}
    }
    b.m_id1 = fragNbr;
    b.m_id2 = 7;
    m_topoMap[idC7] = b;

    /// Add H number 9 and H number 10 and C number 8///
    bool first(true);
    for (int i(0) ; i < m_neighborPedot[idC7].size() ; i ++)
    {
        if(  m_atomPedot[m_neighborPedot[idC7][i]] == H )
        {
            if(first) { idH9 = m_neighborPedot[idC7][i]; first = false;}
            else{ idH10 = m_neighborPedot[idC7][i]; }

        }

        if(  m_atomPedot[m_neighborPedot[idC7][i]] == C )
        {
            idC8 = m_neighborPedot[idC7][i];
        }

    }
    b.m_id1 = fragNbr;
    b.m_id2 = 9;
    m_topoMap[idH9] = b;
    b.m_id1 = fragNbr;
    b.m_id2 = 10;
    m_topoMap[idH10] = b;
    b.m_id1 = fragNbr;
    b.m_id2 = 8;
    m_topoMap[idC8] = b;

    /// Add H number 11 and H number 12 and O number 6 ///
    first = true;
    for (int i(0) ; i < m_neighborPedot[idC8].size() ; i ++)
    {
        if(  m_atomPedot[m_neighborPedot[idC8][i]] == H )
        {
            if(first) { idH11 = m_neighborPedot[idC8][i]; first = false;}
            else{ idH12 = m_neighborPedot[idC8][i]; }

        }

        if(  m_atomPedot[m_neighborPedot[idC8][i]] == O )
        {
            idO6 = m_neighborPedot[idC8][i];
        }

    }
    b.m_id1 = fragNbr;
    b.m_id2 = 11;
    m_topoMap[idH11] = b;
    b.m_id1 = fragNbr;
    b.m_id2 = 12;
    m_topoMap[idH12] = b;
    b.m_id1 = fragNbr;
    b.m_id2 = 6;
    m_topoMap[idO6] = b;


    /// Add C number 4 ///
    for (int i(0) ; i < m_neighborPedot[idO6].size() ; i ++)
    {
        if( m_neighborPedot[idO6][i] != idC8 ) { idC4 = m_neighborPedot[idO6][i];}
    }
    b.m_id1 = fragNbr;
    b.m_id2 = 4;
    m_topoMap[idC4] = b;

    /// Find C number 1 ///
    id = idC2;
    for (int i(0) ; i < m_neighborPedot[idC2].size() ; i ++)
    {
        if( m_atomPedot[m_neighborPedot[idC2][i]] == C && m_neighborPedot[idC2][i] != idC4 ) { idC1 = m_neighborPedot[idC2][i];}
    }

    }

    /// Add Lat H terminal atom of a Pedot chain ///
    id = m_HTerm[1];
    b.m_id1 = m_chainLength - 1;
    b.m_id2 = 13;
    m_topoMap[id] = b;


}

bool crossingRaySphere(const Vect3 &rf1,const Vect3 &rf2,const Vect3 &rfInfluence,const double &sphereR, const Vect3 &boxSize)
{
    Vect3 vrf1(0,0,0);
    Vect3 vrf2 = Vect3::displacementPbc(rf1, rf2, boxSize);
    Vect3 vrfInfluence = Vect3::displacementPbc(rf1, rfInfluence, boxSize);

    // We have to check if the rfInfluence is in between the segment limit or not
    bool isInside(true);
    Vect3 v = Vect3::displacementPbc(rf2, rfInfluence, boxSize);
    double angle1( acos( (vrf2/vrf2.norm()) * (vrfInfluence/vrfInfluence.norm()) ) );
    double angle2( acos( (vrf2/vrf2.norm()) * (v/v.norm()) ) );
    if( (angle1 > pi / 2.0) || (angle2 < pi / 2.0) ) {isInside = false;}

    double d(0);

    if(isInside)
    {
        Vect3 v1((vrfInfluence-vrf1)^(vrfInfluence-vrf2));
        Vect3 v2(vrf2-vrf1);
        d = (v1.norm() / v2.norm());
    }
    else
    {
        d = std::min(vrfInfluence.norm(),v.norm());
    }

    return(d<=sphereR);
}


bool operator==(Vect3 const& v1, Vect3 const& v2)
{
    if(v1.m_x == v2.m_x && v1.m_y == v2.m_y && v1.m_z == v2.m_z)
    {
        return true;
    }
    else
    {
        return false;
    }
}
