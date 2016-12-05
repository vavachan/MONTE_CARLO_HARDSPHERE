using namespace std;
#include "atom.hpp"
#include<tuple>
#include<iostream>
#include<math.h>
#include </home/varghese/Academitialismational/PHD1stsem/Monte_Carlo_fluid/boost_1_62_0/boost/math/special_functions/spherical_harmonic.hpp>

double order_para_global(atom Atoms[],int nAtoms,int l,Vector box) {
    Vector r_ij;
    double rr;
    double shell_one=1.5;
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
    double shell_one=1.5;
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
    double shell_one=1.5;
    long NB=0;
    double Q_L=0;
    //double Q_Li_i[100]= {0.0},Q_Lj_i[100]= {0.0};
    //double Q_Li_r[100]= {0.0},Q_Lj_r[100]= {0.0};
    for(int m=0;m<2*l+1;m++){
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
int check_bond(int i,int j, atom Atoms[],int nAtoms,int l,Vector box) {
    double* Q_lm_i_i;
    double* Q_lm_i_r;
    double* Q_lm_j_i;
    double* Q_lm_j_r;
    Q_lm_i_r = new (nothrow) double[2*l+1];
    Q_lm_j_r = new (nothrow) double[2*l+1];
    Q_lm_i_i = new (nothrow) double[2*l+1];
    Q_lm_j_i = new (nothrow) double[2*l+1];
    q_lm_i(Atoms,nAtoms,j,box,l,Q_lm_j_i,Q_lm_j_r);
    q_lm_i(Atoms,nAtoms,i,box,l,Q_lm_i_i,Q_lm_i_r);
    double d_l_ij_i=0,d_l_ij_r=0;
    double d_mod_i=0;
    double d_mod_j=0;
    for(int m=0; m<2*l+1; m++)
    {
        d_mod_i=d_mod_i+(Q_lm_i_i[m]*Q_lm_i_i[m]+Q_lm_i_r[m]*Q_lm_i_r[m]);
        d_mod_j=d_mod_j+(Q_lm_j_i[m]*Q_lm_j_i[m]+Q_lm_j_r[m]*Q_lm_j_r[m]);

        d_l_ij_i=d_l_ij_i+(Q_lm_i_i[m]*Q_lm_j_r[m]-Q_lm_i_r[m]*Q_lm_j_i[m]);
        d_l_ij_r=d_l_ij_r+(Q_lm_i_r[m]*Q_lm_j_r[m]+Q_lm_i_i[m]*Q_lm_j_i[m]);
    }
    d_mod_i=sqrt(d_mod_i);
    d_mod_j=sqrt(d_mod_j);

    d_l_ij_i=d_l_ij_i/(d_mod_i*d_mod_j);
    d_l_ij_r=d_l_ij_r/(d_mod_i*d_mod_j);
    if(d_l_ij_r > 0.5 )
    {
        delete[] Q_lm_i_i;
        delete[] Q_lm_j_i;
        delete[] Q_lm_i_r;
        delete[] Q_lm_j_r;
        return 1;
    }
    else
    {
        delete[] Q_lm_i_i;
        delete[] Q_lm_j_i;
        delete[] Q_lm_i_r;
        delete[] Q_lm_j_r;
        return 0;
    }
}

