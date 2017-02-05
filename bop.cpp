using namespace std;
#include "atom.hpp"
#include<tuple>
#include<iostream>
#include<math.h>
#include </home/group/boost/boost/math/special_functions/spherical_harmonic.hpp>
double shell_one=1.5;
double xi=10;
double order_para_global(atom Atoms[],int nAtoms,int l,Vector box) {
    Vector r_ij;
    double rr;
    long* Nb;
    Nb= new (nothrow) long [nAtoms];
    if(Nb==nullptr) {
        cout<<"Memory Allocation Failed in bop.hpp\n";
        return 0;
    }
    for(int n=0; n<nAtoms; n++) {
        Nb[n]=0;
    }
    long NB=0;
    double* Q_L_i;
    double* Q_Lj_i;
    double* Q_L_r;
    double* Q_Lj_r;
    Q_L_i= new (nothrow) double[2*l+1];
    Q_Lj_i= new (nothrow) double[2*l+1];
    Q_L_r= new (nothrow) double[2*l+1];
    Q_Lj_r= new (nothrow) double[2*l+1];
    double THETA,PHI,r;
    for(int i=0; i<nAtoms; i++) {
        for(int j=0; j<nAtoms; j++) {
            if(i != j) {
                r_ij=v_sub(Atoms[i].pos,Atoms[j].pos);
                r_ij=VWrap(r_ij,box);
                rr=v_dot(r_ij,r_ij);
                r=sqrt(rr);
                if(r<shell_one) {
                    Nb[i]=Nb[i]+1;
                    NB=NB+1;
                    THETA=acos(r_ij.z/r);
                    if(r_ij.x == 0. and r_ij.y == 0.) {
                        PHI=0;
                    }
                    else {
                        PHI=atan2(r_ij.y,r_ij.x);
                        if(PHI<0)
                            PHI=2*M_PI+PHI;
                    }
                    for(int m=0; m<2*l+1; m++) {
                        Q_Lj_i[m]=Q_Lj_i[m]+boost::math::spherical_harmonic_i <double , double > (l,m-l,THETA, PHI);
                        Q_Lj_r[m]=Q_Lj_r[m]+boost::math::spherical_harmonic_r <double , double > (l,m-l,THETA, PHI);
                    }
                }

            }
        }


        for(int m=0; m<2*l+1; m++) {
            Q_L_i[m]=Q_L_i[m]+Q_Lj_i[m];
            Q_Lj_i[m]=0;
            Q_L_r[m]=Q_L_r[m]+Q_Lj_r[m];
            Q_Lj_r[m]=0;
        }
    }
    double Q_L=0;
    for(int m=0; m<2*l+1; m++) {
        Q_L_i[m]=Q_L_i[m]/NB;
        Q_L_r[m]=Q_L_r[m]/NB;
        Q_L=Q_L+(Q_L_i[m]*Q_L_i[m]+Q_L_r[m]*Q_L_r[m]);
    }
    delete[] Nb;
    delete[] Q_L_r;
    delete[] Q_L_i;
    delete[] Q_Lj_r;
    delete[] Q_Lj_i;
    return sqrt(4*M_PI/(2*l+1)*Q_L);
}

double order_para_local(atom Atoms[],int nAtoms,int l,Vector box) {
    Vector r_ij;
    double rr;
    long* Nb;
    Nb= new (nothrow) long [nAtoms];
    if(Nb==nullptr) {
        cout<<"Memory Allocation Failed in bop.hpp\n";
        return 0;
    }
    for(int n=0; n<nAtoms; n++) {
        Nb[n]=0;
    }
    long NB=0;
    double Q_L=0;
    double* Q_Li_i;
    double* Q_Lj_i;
    double* Q_Li_r;
    double* Q_Lj_r;
    Q_Li_i= new (nothrow) double[2*l+1];
    Q_Lj_i= new (nothrow) double[2*l+1];
    Q_Li_r= new (nothrow) double[2*l+1];
    Q_Lj_r= new (nothrow) double[2*l+1];
    double THETA,PHI,r;
    for(int m=0; m<2*l+1; m++) {
        Q_Li_i[m]=0.0;
        Q_Li_r[m]=0.0;
    }
    for(int i=0; i<nAtoms; i++) {
        for(int m=0; m<2*l+1; m++) {
            Q_Lj_i[m]=0.0;
            Q_Lj_r[m]=0.0;
        }
        //   double Q_Lj_i[100]= {0.0};
        //   double Q_Lj_r[100]= {0.0};
        for(int j=0; j<nAtoms; j++) {
            if(i != j) {
                r_ij=v_sub(Atoms[i].pos,Atoms[j].pos);
                r_ij=VWrap(r_ij,box);
                rr=v_dot(r_ij,r_ij);
                r=sqrt(rr);
                if(r<shell_one) {
                    Nb[i]=Nb[i]+1;
                    NB=NB+1;
                    THETA=acos(r_ij.z/r);
                    if(r_ij.x == 0. and r_ij.y == 0.) {
                        PHI=0;
                    }
                    else {
                        PHI=atan2(r_ij.y,r_ij.x);
                        if(PHI<0)
                            PHI=2*M_PI+PHI;
                    }
                    for(int m=0; m<2*l+1; m++) {
                        Q_Lj_i[m]=Q_Lj_i[m]+boost::math::spherical_harmonic_i <double , double > (l,m-l,THETA, PHI);
                        Q_Lj_r[m]=Q_Lj_r[m]+boost::math::spherical_harmonic_r <double , double > (l,m-l,THETA, PHI);
                    }
                }
            }
        }
        double Q_Li=0;
        for(int m=0; m<2*l+1; m++) {
            Q_Li_i[m]=Q_Lj_i[m]/Nb[i];
            Q_Li_r[m]=Q_Lj_r[m]/Nb[i];
            Q_Li=Q_Li+(Q_Li_i[m]*Q_Li_i[m]+Q_Li_r[m]*Q_Li_r[m]);
        }
        Q_Li=4*M_PI/(2*l+1)*Q_Li;
        Q_L=Q_L+sqrt(Q_Li);
        // cout<<i<<"\t"<<sqrt(Q_Li)<<"\t"<<Nb[i]<<"\n";
    }
    delete[] Nb;
    delete[] Q_Li_i;
    delete[] Q_Lj_i;
    delete[] Q_Li_r;
    delete[] Q_Lj_r;

    return Q_L/nAtoms;
}


void q_lm_i(atom Atoms[],int nAtoms,int i, Vector box,int l,double* Q_Lj_i,double* Q_Lj_r)
{
    Vector r_ij;
    double rr;
    long NB=0;
    double Q_L=0;
    //double Q_Li_i[100]= {0.0},Q_Lj_i[100]= {0.0};
    //double Q_Li_r[100]= {0.0},Q_Lj_r[100]= {0.0};
    for(int m=0; m<2*l+1; m++) {
        Q_Lj_i[m]=0.0;
        Q_Lj_r[m]=0.0;
    }
    double THETA,PHI,r;
    for(int j=0; j<nAtoms; j++) {
        if(i != j) {
            r_ij=v_sub(Atoms[i].pos,Atoms[j].pos);
            r_ij=VWrap(r_ij,box);
            rr=v_dot(r_ij,r_ij);
            r=sqrt(rr);
            if(r<shell_one) {
                NB=NB+1;
                THETA=acos(r_ij.z/r);
                if(r_ij.x == 0. and r_ij.y == 0.) {
                    PHI=0;
                }
                else {
                    PHI=atan2(r_ij.y,r_ij.x);
                    if(PHI<0)
                        PHI=2*M_PI+PHI;
                }
                for(int m=0; m<2*l+1; m++) {
                    Q_Lj_i[m]=Q_Lj_i[m]+boost::math::spherical_harmonic_i <double , double > (l,m-l,THETA, PHI);
                    Q_Lj_r[m]=Q_Lj_r[m]+boost::math::spherical_harmonic_r <double , double > (l,m-l,THETA, PHI);
                }
            }
        }
    }
    for(int m=0; m<2*l+1; m++)
    {
        Q_Lj_i[m]=Q_Lj_i[m]/NB;
        Q_Lj_r[m]=Q_Lj_r[m]/NB;
    }
}
int check_bond(int i,int j, atom Atoms[],int nAtoms,int l,Vector box,double Q_lm_i_r[],double Q_lm_i_i[],double Q_lm_j_r[],double Q_lm_j_i[]) {
    double d_l_ij_i=0,d_l_ij_r=0;
    double d_mod_i=0;
    double d_mod_j=0;
    for(int m=0; m<2*l+1; m++)
    {
        d_mod_i=d_mod_i+(Q_lm_i_i[m]*Q_lm_i_i[m]+Q_lm_i_r[m]*Q_lm_i_r[m]);
        d_mod_j=d_mod_j+(Q_lm_j_i[m]*Q_lm_j_i[m]+Q_lm_j_r[m]*Q_lm_j_r[m]);

        d_l_ij_i=d_l_ij_i+(Q_lm_i_i[m]*Q_lm_j_r[m]-Q_lm_i_r[m]*Q_lm_j_i[m]);
        d_l_ij_r=d_l_ij_r+(Q_lm_i_r[m]*Q_lm_j_r[m]+Q_lm_i_i[m]*Q_lm_j_i[m]);
        //cout<<m<<"\t"<<(Q_lm_i_r[m]*Q_lm_j_r[m]+Q_lm_i_i[m]*Q_lm_j_i[m])<<"\t"<<(Q_lm_i_i[m]*Q_lm_j_r[m]-Q_lm_i_r[m]*Q_lm_j_i[m])<<"\n";
    }
    d_mod_i=sqrt(d_mod_i);
    d_mod_j=sqrt(d_mod_j);

    d_l_ij_i=d_l_ij_i/(d_mod_i*d_mod_j);
    d_l_ij_r=d_l_ij_r/(d_mod_i*d_mod_j);
    //  cout<<d_l_ij_r<<"\n";
    if(d_l_ij_r > 0.7 )
    {
        return 1;
    }
    else
    {
        return 0;
    }
}
long largest_cluster(atom Atoms[],int nAtoms,int l,Vector box) {
    Vector r_ij;
    long* Nb;
    double twob=2*box.x;
    double** Q_local_r;
    double** Q_local_i;
    double* Q_local;
    Q_local=new (nothrow) double [nAtoms];
    Nb= new (nothrow) long [nAtoms];
    Q_local_i= new (nothrow) double* [nAtoms];
    Q_local_r= new (nothrow) double* [nAtoms];
    for(int i=0; i<nAtoms; i++)
    {
        Q_local[i]=0;
        Q_local_i[i] = new (nothrow) double [2*l+1];
        Q_local_r[i] = new (nothrow) double [2*l+1];
    }
    int* Nc;
    int* old_label;
    double rr,r;
    Nc= new (nothrow) int[nAtoms];
    old_label = new (nothrow) int[nAtoms];
    if(Nc==nullptr) {
        cout<<"Memory Allocation Failed in cluster.cpp\n";
        return 0;
    }
    if(Nb==nullptr) {
        cout<<"Memory Allocation Failed in bop.hpp\n";
        return 0;
    }
    for(int n=0; n<nAtoms; n++) {
        Nb[n]=0;
        Nc[n]=0;
    }
    long NB=0;
    double Q_L=0;
    double* Q_Li_i;
    double* Q_Lj_i;
    double* Q_Li_r;
    double* Q_Lj_r;
    Q_Li_i= new (nothrow) double[2*l+1];
    Q_Lj_i= new (nothrow) double[2*l+1];
    Q_Li_r= new (nothrow) double[2*l+1];
    Q_Lj_r= new (nothrow) double[2*l+1];
    double THETA,PHI;
    for(int m=0; m<2*l+1; m++) {
        Q_Li_i[m]=0.0;
        Q_Li_r[m]=0.0;
    }
    for(int i=0; i<nAtoms; i++) {
        for(int m=0; m<2*l+1; m++) {
            Q_Lj_i[m]=0.0;
            Q_Lj_r[m]=0.0;
        }
        //   double Q_Lj_i[100]= {0.0};
        //   double Q_Lj_r[100]= {0.0};
        //

        for(int j=0; j<Atoms[i].neighbours; j++) {
            r_ij.x=Atoms[i].pos.x-Atoms[Atoms[i].neigh_list[j]].pos.x;
            r_ij.y=Atoms[i].pos.y-Atoms[Atoms[i].neigh_list[j]].pos.y;
            r_ij.z=Atoms[i].pos.z-Atoms[Atoms[i].neigh_list[j]].pos.z;
            r_ij.x=(r_ij.x-(twob*lround(r_ij.x/twob)));
            r_ij.y=(r_ij.y-(twob*lround(r_ij.y/twob)));
            r_ij.z=(r_ij.z-(twob*lround(r_ij.z/twob)));
            rr=r_ij.x*r_ij.x+r_ij.y*r_ij.y+r_ij.z*r_ij.z;
            //    for(int j=0; j<nAtoms; j++) {
            //      if(i != j) {
            //        r_ij=v_sub(Atoms[i].pos,Atoms[j].pos);
            //      r_ij=VWrap(r_ij,box);
            //    rr=v_dot(r_ij,r_ij);
            r=sqrt(rr);
            if(r<shell_one) {
                Nb[i]=Nb[i]+1;
                NB=NB+1;
                THETA=acos(r_ij.z/r);
                if(r_ij.x == 0. and r_ij.y == 0.) {
                    PHI=0;
                }
                else {
                    PHI=atan2(r_ij.y,r_ij.x);
                    if(PHI<0)
                        PHI=2*M_PI+PHI;
                }
                for(int m=0; m<2*l+1; m++) {
                    Q_Lj_i[m]=Q_Lj_i[m]+boost::math::spherical_harmonic_i <double , double > (l,m-l,THETA, PHI);
                    Q_Lj_r[m]=Q_Lj_r[m]+boost::math::spherical_harmonic_r <double , double > (l,m-l,THETA, PHI);
                }
            }
        }

        double Q_Li=0;
        for(int m=0; m<2*l+1; m++) {
            Q_Li_i[m]=Q_Lj_i[m]/Nb[i];
            Q_Li_r[m]=Q_Lj_r[m]/Nb[i];
            Q_Li=Q_Li+(Q_Li_i[m]*Q_Li_i[m]+Q_Li_r[m]*Q_Li_r[m]);
            Q_local_i[i][m]=Q_Li_i[m];
            Q_local_r[i][m]=Q_Li_r[m];
        }
        Q_Li=4*M_PI/(2*l+1)*Q_Li;
        Q_local[i]=sqrt(Q_Li);
        Q_L=Q_L+sqrt(Q_Li);
    }
    for(int n=0; n<nAtoms; n++) {
        Nc[n]=0;
    }
//  for(int i=0; i<nAtoms; i++) {
//      cout<<i<<"\t"<<Atoms[i].cluster_index<<"\n";
//  }
    for(int i=0; i<nAtoms; i++) {

        for(int j=0; j<Atoms[i].neighbours; j++) {
            r_ij.x=Atoms[i].pos.x-Atoms[Atoms[i].neigh_list[j]].pos.x;
            r_ij.y=Atoms[i].pos.y-Atoms[Atoms[i].neigh_list[j]].pos.y;
            r_ij.z=Atoms[i].pos.z-Atoms[Atoms[i].neigh_list[j]].pos.z;
            r_ij.x=(r_ij.x-(twob*lround(r_ij.x/twob)));
            r_ij.y=(r_ij.y-(twob*lround(r_ij.y/twob)));
            r_ij.z=(r_ij.z-(twob*lround(r_ij.z/twob)));
            rr=r_ij.x*r_ij.x+r_ij.y*r_ij.y+r_ij.z*r_ij.z;
//  for(int j=i+1; j<nAtoms; j++)
//  {

//      r_ij=v_sub(Atoms[i].pos,Atoms[j].pos);
//      r_ij=VWrap(r_ij,box);
//      rr=v_dot(r_ij,r_ij);
            r=sqrt(rr);
            if(r<shell_one) {
                //	cout<<i<<"\t"<<Atoms[i].cluster_index<<"\n";
                Atoms[i].close_neighbours+=1;
                Atoms[i].close_update_neighbour(Atoms[i].neigh_list[j]);
	//	cout<<i<<"\t"<<Atoms[i].neigh_list[j]<<"\n";
//		cout<<i<<"\t"<<j<<"\n";
                if(check_bond(i,Atoms[i].neigh_list[j],Atoms,nAtoms,l,box,Q_local_r[i],Q_local_i[i],Q_local_r[Atoms[i].neigh_list[j]],Q_local_i[Atoms[i].neigh_list[j]])) {
                    Atoms[i].connections+=1;
                    //      Atoms[j].connections+=1;
                }
            }
        }
    }
    for(int i=0; i<nAtoms; i++) {
//	cout<<i<<"\t"<<Atoms[i].connections<<"\n";
//                  cout<<i<<"\t"<<Q_local[i]<<"\t"<<Nb[i]<<"\n";
          if(Atoms[i].connections>xi) {
//        if(Nb[i]==4 and Q_local[i] >= 0.6) {
            Atoms[i].cluster_index=i;
           // cout<<i<<"\t"<<Q_local[i]<<"\t"<<Nb[i]<<"\t"<<Atoms[i].close_neighbours<<"\n";
        }
    }
    int min;
    int change=1;
    int flag=1;
    while(change) {
        for(int i=0; i<nAtoms; i++) {
            old_label[i]=Atoms[i].cluster_index;
        }
        for(int i=0; i<nAtoms; i++) {
              if(Atoms[i].connections>xi) {
          //  if(Nb[i]==4 and Q_local[i] >= 0.6) {
                min=Atoms[i].cluster_index;
                for(int n=0; n<Atoms[i].close_neighbours; n++) {
                       if(Atoms[Atoms[i].close_neigh_list[n]].connections>xi){
		   // cout<<i<<"\t"<<n<<"\t"<<Atoms[i].close_neigh_list[n]<<"\n";
            //        if(Nb[Atoms[i].close_neigh_list[n]]==4 and Q_local[Atoms[i].close_neigh_list[n]] >= 0.6) {
                        if(min>Atoms[Atoms[i].close_neigh_list[n]].cluster_index) {
                            min=Atoms[Atoms[i].close_neigh_list[n]].cluster_index;
                        }
                    }
                }
                for(int n=0; n<Atoms[i].close_neighbours; n++)
                {
                    if(Atoms[Atoms[i].close_neigh_list[n]].connections>xi)
                 //   if(Nb[Atoms[i].close_neigh_list[n]]==4 and Q_local[Atoms[i].close_neigh_list[n]] >= 0.6)
                    {
                        Atoms[Atoms[i].close_neigh_list[n]].cluster_index=min;
                    }
                }
            }
        }
        change=0;
        for(int i=0; i<nAtoms; i++) {
            flag=(old_label[i]!=Atoms[i].cluster_index);
            if(flag) {
                change=1;
                break;
            }
        }
    }
    long max=0;
    for(int i=0; i<nAtoms; i++)
    {  // cout<<i<<"\t"<<Atoms[i].cluster_index<<"\n";
        if(Atoms[i].cluster_index!=-1)
            Nc[Atoms[i].cluster_index]++;
    }
    max=Nc[0];
    for(int i=0; i<nAtoms; i++)
        if(Nc[i]>max)
            max=Nc[i];
    delete[] Nc;
    delete[] Nb;
    delete[] Q_Li_i;
    delete[] Q_Lj_i;
    delete[] Q_Li_r;
    delete[] Q_Lj_r;
    for(int i=0; i<nAtoms; i++)
    {
        delete [] Q_local_i[i];
        delete [] Q_local_r[i];
    }
    delete[] Q_local_i;
    delete[] Q_local_r;
    delete[] Q_local;
    delete[] old_label;
    return max;
}


