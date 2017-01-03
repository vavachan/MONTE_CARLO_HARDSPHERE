#include "atom.hpp"
#include<iostream>
void atom::update_neighbour(int neigh)	{
     neigh_list[index]=neigh;
     index++;
}
void atom::reset()	{
	for(int i=0;i<15;i++)
		neigh_list[i]=0;
	index=0;
}
