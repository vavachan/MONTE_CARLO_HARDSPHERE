using namespace std;
#include "atom.hpp"
#include<fstream>
#include<math.h>
#include<iostream>
double pair_correlation(atom* Atoms,int nAtoms,int END,Vector box,double temp) {
    double r=0,rr;
    static int g_r[1000]= {0};
    double bin_width=0.01;//sqrt(v_dot(box,box))/1000.0;
    static int counter=0;
    static double prefact=4.0/3*3.14;
    double density=nAtoms/(2*box.x*2*box.y*2*box.z);
    char buffer[32];
    double R_CUT_HS=1.0;
    int bin_no=0;
    Vector dr;
    for(int i=0; i<nAtoms-1; i++)
        for(int j=i+1; j<nAtoms; j++) {
            dr=v_sub(Atoms[i].pos,Atoms[j].pos);
            dr=VWrap(dr,box);
            rr=v_dot(dr,dr);
            r=sqrt(rr);
	    if(r>9.99)
		continue;
            bin_no=int(r/bin_width);
            g_r[bin_no]=g_r[bin_no]+2;
        }
    counter=counter+1;
    double r_lower,r_upper,omega;
    if(END) {
	snprintf(buffer,sizeof(char)*32,"OUT/g_r.dat");
	std::ofstream G_R(buffer);
        for(int i=0; i<1000; i++) {
            r_lower=i*bin_width;
            r_upper=(i+1)*bin_width;
            omega=prefact*(pow(r_upper,3)-pow(r_lower,3));
            G_R<<(i)*bin_width<<"\t"<<g_r[i]/(counter*nAtoms*omega)*(1./density)<<"\n";
        }
        G_R.close();
    }
   int d=R_CUT_HS/bin_width;
   r_lower=d*bin_width;
   r_upper=(d+1)*bin_width;
   omega=prefact*(pow(r_upper,3)-pow(r_lower,3));
   return g_r[d]/(counter*nAtoms*omega)*(1./density);
    
}
