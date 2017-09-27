#ifndef PARSE_H_INCLUDED
#define PARSE_H_INCLUDED

#include "Graph.h"

void parseBigPdb(std::string fileI, std::string fileO, int nbPdt, int nbTos, int nbAtPedot, int nbAtTos);

void extendPdb(std::string fileI);

void extendTop(std::string fileI);

#endif // PARSE_H_INCLUDED
