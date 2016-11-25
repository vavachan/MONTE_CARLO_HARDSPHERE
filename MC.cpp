using namespace std;
#include "MC_header.hpp"
#include "ini_pos.hpp"
#include "g_r.hpp"
double int_pot (double rr) {
    return 4.0*(1.0/pow(rr,6)-1.0/pow(rr,3));
}


double cal_pot (atom Atoms[],int nAtoms,int i) {
    Vector dr;
    double rr;
    double PE=0;
    static double P_cut=int_pot(pow(r_cut,2));
    for(int j=0; j<nAtoms; j++)
    {
        dr=v_sub(Atoms[j].pos,Atoms[i].pos);
        dr=VWrap(dr,box);
        rr=v_dot(dr,dr);
        if (rr<pow(r_cut,2) and rr) {
            PE=PE+4.0*(1.0/pow(rr,6)-1.0/pow(rr,3))-P_cut;
        }
    }
    return PE;
}
void print_pos (atom Atoms[],int nAtoms) {
    for (int i=0; i<nAtoms; i++)
        POSITION<<i<<"\t"<<Atoms[i].pos.x<<"\t"<<Atoms[i].pos.y<<"\t"<<Atoms[i].pos.z<<"\n";
}

void pot_energy(atom Atoms[],int nAtoms,int t) {
    for(int i=0; i<nAtoms; i++) {
        //PE=PE+0.5*v_dot(Atoms[i].vel,Atoms[i].vel);
        if(fmod(t,100)==0)
            MSD<<t<<"\t"<<Atoms[i].pos.x*sigma<<"\t"<<Atoms[i].pos.y*sigma<<"\t"<<Atoms[i].pos.z*sigma<<"\n";
    }
    //return PE;
}
void mcmove(atom Atoms[],int nAtoms) {
    double POTn,POTo,P,r;
    std::uniform_int_distribution<int> uni(0,nAtoms);
    auto random_integer = uni(rng);
    Vector random_dir;
    random_dir.x=uni_d(rng)-0.5;
    random_dir.y=uni_d(rng)-0.5;
    random_dir.z=uni_d(rng)-0.5;
    POTo=cal_pot(Atoms,nAtoms,random_integer);
    Atoms[random_integer].pos=v_Sadd(Atoms[random_integer].pos,random_dir,JUM);
    POTn=cal_pot(Atoms,nAtoms,random_integer);
    P=exp(-1.0/temp*(POTn-POTo));
    //cout<<POTo<<"\t"<<POTn<<"\t"<<random_integer<<"\n";
    r=uni_d(rng);
    if(P>1)
    {
        return;
    }
    else if (r<P)
    {
        return;
    }
    else
    {   Atoms[random_integer].pos=v_Sadd(Atoms[random_integer].pos,random_dir,-JUM);
        return;
    }

}
int main() {
    int nAtoms,N;
    atom Atoms[1000];
    cout<<"Enter the Box Dimenions:\n";
    cin>>box.x>>box.y>>box.z;
    cout<<"Enter the number of atoms:\n";
    cin>>nAtoms;
    cout<<"Enter how long you want to run the simulation:\n";
    cin>>N;
    cout<<"temp\n";
    cin>>temp;
    density=nAtoms/(2*box.x*2*box.y*2*box.z);
    cout<<"density:"<<nAtoms/(2*box.x*2*box.y*2*box.z);
    inipos(Atoms,nAtoms,box);
    clock_t begin=clock();
    int END=0;
    for(int i=0; i<N; i++) {
        for(int n=0; n<nAtoms; n++)
            mcmove(Atoms,nAtoms);
        pot_energy(Atoms,nAtoms,i);
        if(i>50000) {
            pair_correlation(Atoms,nAtoms,0,box);
            if(i==N-1)
                pair_correlation(Atoms,nAtoms,1,box);
        }
    }
    clock_t end =clock();
    double elapsed_time= double (end-begin)/CLOCKS_PER_SEC;
    cout<<"\n"<<elapsed_time;
    print_pos(Atoms,nAtoms);
    return 0;
}


