using namespace std;
#include "atom.hpp"
#include<iostream>
#include<math.h>
#include </home/varghese/Academitialismational/PHD1stsem/Monte_Carlo_fluid/boost_1_62_0/boost/math/special_functions/spherical_harmonic.hpp>

double order_para(atom Atoms[],int nAtoms,int l,Vector box) {
    Vector r_ij;
    double rr;
    double shell_one=2.95;
    long Nb[5000],NB=0;
    double Q_L_i[100]= {0.0},Q_Lj_i[100]= {0.0};
    double Q_L_r[100]= {0.0},Q_Lj_r[100]= {0.0};
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
    return sqrt(4*M_PI/(2*l+1)*Q_L);
}



