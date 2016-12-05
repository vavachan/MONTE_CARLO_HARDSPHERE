using namespace std;
#include <tuple>
#include "atom.hpp"
#include "bop.hpp"
#include<iostream>
#include<math.h>
long cluster(atom Atoms[],int nAtoms,int l,Vector box) {
    int* Nc;
    double shell_one=1.5;
    double rr,r;
    Vector r_ij;
    Nc= new (nothrow) int[nAtoms];
    if(Nc==nullptr) {
        cout<<"Memory Allocation Failed in cluster.cpp\n";
        return 0;
    }
    for(int n=0; n<nAtoms; n++) {
        Nc[n]=0;
    }
    for(int i=0; i<nAtoms; i++) {
        if(Atoms[i].cluster_index==-1) {
            Atoms[i].cluster_index=i;
        }
        for(int j=i+1; j<nAtoms; j++)
        {

            r_ij=v_sub(Atoms[i].pos,Atoms[j].pos);
            r_ij=VWrap(r_ij,box);
            rr=v_dot(r_ij,r_ij);
            r=sqrt(rr);
            if(r<shell_one) {
                if(Atoms[j].cluster_index == -1) {
                    if(check_bond(i,j,Atoms,nAtoms,l,box)) {
                        Atoms[j].cluster_index=Atoms[i].cluster_index;
                    }
                }
            }
        }
    }
    long max=0;
    for(int i=0; i<nAtoms; i++)
    {
        Nc[Atoms[i].cluster_index]++;
    }
    max=Nc[0];
    for(int i=0; i<nAtoms; i++)
        if(Nc[i]>max)
            max=Nc[i];
    delete[] Nc;
    return max;
}
