#ifndef SYSTEM_H_INCLUDED
#define SYSTEM_H_INCLUDED
#include "Units.h"

struct System
{
    std::string m_line1, m_line2, m_line3, m_line4;
    double m_boxSize;
    std::vector<Polymer> m_pedot;
    std::vector<RigidFragment> m_tos;
    Topo m_topo;

    // Vector for crystallites structure
    std::vector<std::vector<int>> m_crystalMol;
    std::vector<int> m_crystalId;

    void readSystem(std::string fileNameTop, std::string fileNameCoord);
    void initGeom(); // Compute rigid fragment properties (center and orientation)
    void writeSystem(std::string fileNameCoord);
    void writeSystemGro(std::string fileNameGro = "");
    void boundDistribution(std::string boundDistrib);


};

struct Coarsegrain : System
{
    std::vector<std::vector<int>> m_crystalMol;
    std::vector<int> m_crystalId;

    void readSystem(std::string fileNameCoord);
    void initGeom(); // Compute rigid fragment properties (center and orientation)
    void writeSystem(std::string fileNameCoord);
    void writeFragmentOrient(std::string fileVector);
    void arrangeCrystalMol(); // Rearrange structure vector (crystalMol and crystalId)
    void computeCrystalDist(std::string fileNameDist); // For each pair of crystallite, compute all beads to beads distance
};

int newBoxPbc(int box, const Vect3 &pbc);

int lookupTabPbc(int box, Direction dir);


#endif // SYSTEM_H_INCLUDED
