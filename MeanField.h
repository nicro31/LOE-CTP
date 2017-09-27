#ifndef MEANFIELD_H_INCLUDED
#define MEANFIELD_H_INCLUDED
#include "Parse.h"

#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <fstream>
#include <time.h>
#include <omp.h>

double adjustConcentration(double &pTemp, Graph &g); // Set the chemical potential to get a desired carrier concentration
double computeMobility(Graph &g, bool usePrevDistrib, std::string folder, int u); // Compute graph mobility, for a graph with parameters already set (field, temperature, concentration)
Graph createCubic(int dim, double latticeParam); // Create graph on a cubic lattice (debug purpose)
Vect3 findPercolationThreshold(Graph &g, std::string folder); // Find the percolation threshold (in x, y, z) for a graph
double findPercolationThresholdN(Graph &g, std::string folder); // Find the percolation threshold in nb of molecules per cluster
void testCubic(); // Run calculations on a cubic graph (debug purpose)
void prepareForGorilla(std::string fileNameTop, std::string fileNamecoord, std::string folder); // Prepare pdb files for gorilla
void plotPercolation(Graph &g, double threshMin, double threshMax, double threshStep, std::string folder, bool dimCluster); // Compute the percolation behavior of the graph
void computeMobilityDependance(Graph &g, std::vector<double> pField, std::vector<double> pTemperature, std::vector<double> pConcentration, std::string folder); // Compute the mobility for different values of the parameters field, tempearture and concentration
Graph initGraph(std::string fileNameTop, std::string fileNameCoord, std::string configFile, unsigned int seed, std::string folder = "", bool readPhysics = true); // Prepare a graph from pbd file coordFile, with the top file fileNameTop, with parameters in configFile, and output results in folder
void readMobilityParam(std::string configFile, std::vector<double> &pField, std::vector<double> &pTemperature, std::vector<double> &pConcentration);

/// Coarse grain applications ///
void systemMorpho(std::string fileNameCoord); // Compute pi-pi stacking and lamellar properties of the MD coarse grained volume

#endif // MEANFIELD_H_INCLUDED
