#ifndef JOHNSON_H_INCLUDED
#define JOHNSON_H_INCLUDED
#include <iostream>
#include <fstream>
#include <vector>
#include <climits>
#include <exception>
#include <set>
#include <iomanip>
#include <cmath>

#define eps 1e-10

struct Edge
{
  int head;
  double cost;
};

using SingleSP = std::vector<double>;
using AllSP = std::vector<std::vector<double>>;
const double INF = 1e3;

std::vector<std::vector<Edge>> loadgraph(std::istream& is);

std::vector<std::vector<Edge>> addZeroEdge(std::vector<std::vector<Edge>> g);

SingleSP bellmanford(std::vector<std::vector<Edge>> &g, int s);

SingleSP djikstra(const std::vector<std::vector<Edge>>& g, int s);

AllSP johnson(std::vector<std::vector<Edge>> &g);


#endif // JOHNSON_H_INCLUDED
