using namespace std;
#include "atom.hpp"
#include<fstream>
#include<math.h>
#include<iostream>
void pair_correlation(atom* Atoms,int nAtoms,int END,Vector box,double temp) {
    double r=0,rr;
    static int g_r[1000]= {0};
    double bin_width=sqrt(v_dot(box,box))/1000.0;
    static int counter=0;
    static double prefact=4.0/3*3.14;
    double density=nAtoms/(2*box.x*2*box.y*2*box.z);
    char buffer[32];
    int bin_no=0;
    Vector dr;
    for(int i=0; i<nAtoms-1; i++)
        for(int j=i+1; j<nAtoms; j++) {
            dr=v_sub(Atoms[i].pos,Atoms[j].pos);
            dr=VWrap(dr,box);
            rr=v_dot(dr,dr);
            r=sqrt(rr);
            bin_no=int(r/bin_width);
            g_r[bin_no]=g_r[bin_no]+2;
        }
    counter=counter+1;
    if(END) {
	snprintf(buffer,sizeof(char)*32,"g_r_%f.dat",temp);
	std::ofstream G_R(buffer);
        double r_lower,r_upper,omega;
        for(int i=0; i<1000; i++) {
            r_lower=i*bin_width;
            r_upper=(i+1)*bin_width;
            omega=prefact*(pow(r_upper,3)-pow(r_lower,3));
            G_R<<(i)*bin_width<<"\t"<<g_r[i]/(counter*nAtoms*omega)*(1./density)<<"\n";
        }
        G_R.close();
    }

}
