#ifndef BOP_INCLUDED_H
#define BOP_INCLUDED_H
double order_para_global(atom Atoms[],int nAtoms,int l,Vector box); 
double order_para_local(atom Atoms[],int nAtoms,int l,Vector box); 
void q_lm_i(atom Atoms[],int nAtoms,int i, Vector box,int l);
int check_bond(int i,int j, atom Atoms[],int nAtoms,int l,Vector box,double*,double*,double*,double*); 
long largest_cluster(atom Atoms[],int nAtoms,int l,Vector box,char label[]);
#endif


