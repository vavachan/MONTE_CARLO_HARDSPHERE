#ifndef ATOM_INCLUDED_H
#define ATOM_INCLUDED_H
#include "vector.hpp"
class atom {
public:
    Vector pos;
    int cluster_index=-1;
    int connections=0;
    int index=0;
    int neigh_list[200]= {0};
    void update_neighbour(int neigh);
    void reset();
};
#endif
