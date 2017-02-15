using namespace std;
#include "atom.hpp"
#include<iostream>
#include<math.h>
#include<iomanip>
#include<stdlib.h>
#include<random>
//#include "MC_header.hpp"
int find_cube(int nAtoms) {
    for(int i=0; 1; i++)
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
void random_ini(atom Atoms[],int nAtoms,Vector box,double R_CUT_HS)
{
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_real_distribution<long double> uni_d(0,1);
    long double random_pos;
    int flag=1;
    Vector dr;
    long double rr;
    long double twob=2*box.x;
    int count=0;
    for(int i=0; i<nAtoms;) {
        flag=1;
        Atoms[i].pos.x=(drand48()-0.5)*(twob);
        Atoms[i].pos.y=(drand48()-0.5)*(twob);
        Atoms[i].pos.z=(drand48()-0.5)*(twob);
        for(int j=0; j<i; j++)
        {
            dr.x=Atoms[i].pos.x-Atoms[j].pos.x;
            dr.y=Atoms[i].pos.y-Atoms[j].pos.y;
            dr.z=Atoms[i].pos.z-Atoms[j].pos.z;
            dr.x=(dr.x-(twob*lround(dr.x/twob)));
            dr.y=(dr.y-(twob*lround(dr.y/twob)));
            dr.z=(dr.z-(twob*lround(dr.z/twob)));
            rr=dr.x*dr.x+dr.y*dr.y+dr.z*dr.z;
            //   dr.x=fmod(dr.x,twob);
            //   dr.y=fmod(dr.y,twob);
            //   dr.z=fmod(dr.z,twob);
            //   if (dr.x >= box.x) {
            //       dr.x -=  twob;
            //   }
            //   else if (dr.x < -box.x) {
            //       dr.x += twob;
            //   }

            //   if (dr.y >= box.y) {
            //       dr.y -= twob;
            //   }
            //   else if (dr.y < -box.y) {
            //       dr.y += twob;
            //   }

            //   if (dr.z >= box.z) {
            //       dr.z -=twob;
            //   }
            //   else if (dr.z < -box.z) {
            //       dr.z += twob;
            //   }
            //   dr=v_sub(Atoms[i].pos,Atoms[j].pos);
            //   dr=VWrap(dr,box);
            //     rr=v_dot(dr,dr);
            if(rr<R_CUT_HS) {
                flag=0;
                count+=1;
                break;
            }

        }
        if(flag)
        {

            cout<<i<<"\t"<<Atoms[i].pos.x<<"\t"<<Atoms[i].pos.y<<"\t"<<Atoms[i].pos.z<<"\t"<<count<<"\n";
            i++;
            count=0;
        }

    }

}
Vector FCC(atom Atoms[],int nAtoms,Vector box,double a) {
    Vector v,X;
    int count=0;
    int n=4;
    box.x=(n)*a/2;
    box.y=(n)*a/2;
    box.z=(n)*a/2;
   // cout<<"\n"<<box.x<<"\t"<<box.y<<"\t"<<box.z<<"\n";
    for(int i=0; i<n; i++)
        for(int j=0; j<n; j++)
            for(int k=0; k<n; k++)
            {
                X.x=i*a;
                X.y=j*a;
                X.z=k*a;
		//cout<<count<<"\t"<<X.x<<"\t"<<X.y<<"\t"<<X.z<<"\n";
                Atoms[count].pos=v_add(X,v);
                Atoms[count].pos=v_sub(Atoms[count].pos,box);
                count=count+1;
            }
    v.x=0.;
    v.y=a/2.;
    v.z=a/2.;
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

    v.x=a/2.;
    v.y=0.;
    v.z=a/2.;
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
    v.x=a/2.;
    v.y=a/2.;
    v.z=0.;
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
