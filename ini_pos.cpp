using namespace std;
#include "atom.hpp"
#include<iostream>
#include<math.h>
#include<iomanip>
int find_cube(int nAtoms) {
    for(int i=0;1; i++)
        if(i*i*i>=nAtoms) {
            return i;
            break;
        }
}
void inipos(atom * Atoms,int nAtoms,Vector box) {
    int cube=find_cube(nAtoms);
    int n=0;
    Vector grid;
    grid.x=2*box.x/(cube);
    grid.y=2*box.y/(cube);
    grid.z=2*box.z/(cube);
 // cout<<grid.x<<"\n";
 // cout<<nAtoms<<"\n";
    for(int i=0; i<cube; i++)
        for(int j=0; j<cube; j++)
            for(int k=0; k<cube; k++) {
                if(n>nAtoms-1)
                    break;
                Atoms[n].pos.x=-box.x+grid.x*i+0.5*grid.x;
                Atoms[n].pos.y=-box.y+grid.y*j+0.5*grid.y;
                Atoms[n].pos.z=-box.z+grid.z*k+0.5*grid.z;
                n++;
            }
}
Vector inipos_bcc(atom Atoms[],int nAtoms,Vector box,double a) {
    Vector v,X;
    int count=0;
    int n=4;
    box.x=(n)*a/2;
    box.y=(n)*a/2;
    box.z=(n)*a/2;
    cout<<"\n"<<box.x<<"\t"<<box.y<<"\t"<<box.z<<"\n";
    for(int i=0; i<n; i++)
        for(int j=0; j<n; j++)
            for(int k=0; k<n; k++)
                {
                    X.x=i*a;
                    X.y=j*a;
                    X.z=k*a;
                    Atoms[count].pos=v_add(X,v);
                    Atoms[count].pos=v_sub(Atoms[count].pos,box);
                    count=count+1;
                }
    v.x=a/2;
    v.y=a/2;
    v.z=a/2; 
    for(int i=0; i<n; i++)
        for(int j=0; j<n; j++)
            for(int k=0; k<n; k++)
                {
                    X.x=i*a;
                    X.y=j*a;
                    X.z=k*a;
                    Atoms[count].pos=v_add(X,v);
                    Atoms[count].pos=v_sub(Atoms[count].pos,box);
                    count=count+1;
                }
     cout<<"\nnatoms:"<<"\t"<<count<<"\n";
     return box;
}


