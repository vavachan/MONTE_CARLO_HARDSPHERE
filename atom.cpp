#include "atom.hpp"
#include<iostream>
void atom::update_neighbour(int neigh)	{
     neigh_list[index]=neigh;
     index++;
}
void atom::close_update_neighbour(int neigh)	{
     close_neigh_list[close_index]=neigh;
     close_index++;
}
void atom::close_reset()	{
//	neighbours=0;
	close_neighbours=0;
////////for(int i=0;i<180;i++)
////////	neigh_list[i]=0;
        for(int i=0;i<50;i++)
        	close_neigh_list[i]=0;
//	index=0;
	close_index=0;
}
void atom::reset()	{
	neighbours=0;
//	close_neighbours=0;
        for(int i=0;i<180;i++)
        	neigh_list[i]=0;
////////for(int i=0;i<50;i++)
////////	close_neigh_list[i]=0;
	index=0;
//	close_index=0;
}
