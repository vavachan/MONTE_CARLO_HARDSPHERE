#ifndef ATOM_INCLUDED_H
#define ATOM_INCLUDED_H
#include "vector.hpp"
class atom {
public:
    Vector pos;
    Vector old_pos;
    double dist=0;
    int cluster_index=-1;
    int connections=0;
    int neighbours=0;
    int close_neighbours=0;
    int index=0;
    int close_index=0;
    int neigh_list[180]= {0};
    int close_neigh_list[50]= {0};
    void update_neighbour(int neigh);
    void close_update_neighbour(int neigh);
    void reset();
    void close_reset();
};
#endif
