using namespace std;
#include "MC_header.hpp"
#include "ini_pos.hpp"
#include "g_r.hpp"
double int_pot (double rr) {
    return 4.0*epsilon*(pow(sigma2/rr,6)-pow(sigma2/rr,3));
}
double force(double rr) {
    return 24.0*epsilon*(2.0*pow(sigma2/rr,6)-pow(sigma2/rr,3));
}
double pertail(double RCUT) {
    return 1*16.0/3.0*M_PI*density*density*(2./3.*(pow(sigma/RCUT,9)-pow(sigma/RCUT,3)));
}
double cal_pot (atom Atoms[],int nAtoms) {
    Vector dr;
    double rr,r;
    double PE=0;
    static double P_cut=int_pot(pow(r_cut,2));
    vir_pre=0;
    for(int i=0; i<nAtoms-1; i++)
    {
        for(int j=i+1; j<nAtoms; j++)
        {
            dr=v_sub(Atoms[j].pos,Atoms[i].pos);
            dr=VWrap(dr,box);
            rr=v_dot(dr,dr);
            if (rr<RCUT_sq and rr) {
                PE+=4.0*(1.0/pow(rr,6)-1.0/pow(rr,3));
                vir_pre+=force(rr);
                //            cout<<i<<"\t"<<j<<"\t"<<rr<<"\n";
                //PE=PE+hd_pot(rr);
            }
        }
        //   cout<<i<<"\t"<<PE<<"\n";
    }
    return PE;
}
void print_pos (atom Atoms[],int nAtoms) {
    std::ofstream POSITION("config.dat");
    for (int i=0; i<nAtoms; i++)
        POSITION<<i<<"\t"<<Atoms[i].pos.x<<"\t"<<Atoms[i].pos.y<<"\t"<<Atoms[i].pos.z<<"\n";
    POSITION.close();
}
double mcmove(atom Atoms[],int nAtoms) {
    double POTn,POTo,P,r;
    POTo=PE;
    std::uniform_int_distribution<int> uni(0,nAtoms-1);
    auto random_integer = uni(rng);
    Vector random_dir;
    random_dir.x=uni_d(rng)-0.5;
    random_dir.y=uni_d(rng)-0.5;
    random_dir.z=uni_d(rng)-0.5;
    //  POTo=cal_pot(Atoms,nAtoms,random_integer);
    Atoms[random_integer].pos=v_Sadd(Atoms[random_integer].pos,random_dir,JUM);
    POTn=cal_pot(Atoms,nAtoms);
    P=exp(-1.0/(k_b*temp)*(POTn-POTo));
    //cout<<POTo<<"\t"<<POTn<<"\t"<<random_integer<<"\n";
    Iter++;
    r=uni_d(rng);
    if(P>1)
    {
        Nacc++;
        return POTn;
    }
    else if (r<P)
    {
        Nacc++;
        return POTn;
    }
    else
    {   Atoms[random_integer].pos=v_Sadd(Atoms[random_integer].pos,random_dir,-JUM);
        return POTo;
    }
    //free(Atoms);

}
double vmove(atom Atoms[],int nAtoms) {
    double old_vol,lnV0,lnV;
    double Lnew,f,r,P;
    double PE_new;
    old_vol=vol;
    lnV0=log(old_vol);
    lnV=lnV0+dlnV*(uni_d(rng)-0.5);
    vol=exp(lnV);
    Lnew=pow(vol,0.3333333)/2.0;
    f=Lnew/box.x;
    box.x *= f;
    box.y *= f;
    box.z *= f;
    for (int n=0; n<nAtoms; n++) {
        Atoms[n].pos.x *= f;
        Atoms[n].pos.y *= f;
        Atoms[n].pos.z *= f;
    }
    print_pos(Atoms,nAtoms);
    r_cut*=f;
    RCUT_sq*=f*f;
    PE_new=cal_pot(Atoms,nAtoms);
    P=exp(-1.0/(k_b*temp)*((PE_new-PE)+Press*(vol-old_vol)-(nAtoms+1)*k_b*temp*log(vol/old_vol)));
    r=uni_d(rng);
    Iter_v++;
    if(r<P)
    {
        Nacc_v++;
        return PE_new;
    }
    else
    {
        f=1.0/f;
        for (int n=0; n<nAtoms; n++) {
            Atoms[n].pos.x *= f;
            Atoms[n].pos.y *= f;
            Atoms[n].pos.z *= f;
        }

        box.x *= f;
        box.y *= f;
        box.z *= f;
        vol=old_vol;
        r_cut*=f;
        RCUT_sq*=f*f;
        return PE;
    }
}
int main() {
    int nAtoms,N;
    atom * Atoms;
    cout<<"Enter the Box Dimenions:\n";
    cin>>box.x>>box.y>>box.z;
    cout<<"Enter the number of atoms:\n";
    cin>>nAtoms;
    cout<<"Enter how long you want to run the simulation:\n";
    cin>>N;
    cout<<"temp\n";
    cin>>temp;
    cin>>Press;
    Atoms = new (nothrow) atom[nAtoms];
    cout<<sizeof(Atoms)<<"\n";
    if(Atoms==nullptr) {
        cout<<"Memory Allocation Failed\n";
        return 0;
    }
//    box=inipos_bcc(Atoms,nAtoms,box,2.5);
    inipos(Atoms,nAtoms,box);
    print_pos(Atoms,nAtoms);
    density=nAtoms/(2*box.x*2*box.y*2*box.z);
    cout<<"density:"<<nAtoms/(2*box.x*2*box.y*2*box.z)<<"\n";
    vol=2*box.x*2*box.y*2*box.z;
    PE=cal_pot(Atoms,nAtoms);
    cout<<PE<<"\n";
    double EqN=10000;
    int END=0;
    clock_t begin=clock();
    int rand;
    cout<<Press<<"\n";
    PE=vmove(Atoms,nAtoms);
    for(int i=0; i<EqN; i++) {
        std::uniform_int_distribution<int> uni(0,nAtoms);
        rand=uni(rng);
    //    cout<<i<<"\t"<<PE<<"\n";
        if(rand<nAtoms) {
            for(int n=0; n<nAtoms; n++)
                PE=mcmove(Atoms,nAtoms);
        }
        else
        {
            PE=vmove(Atoms,nAtoms);
            density=nAtoms/vol;
        }
        if(fmod(i,100)==0)
        {   cout<<i<<"\n";
	}
    }
    //        mcmove(Atoms,nAtoms);
    Nacc=0;
    Iter=0;
    Nacc_v=0;
    Iter_v=0;
    double P_sum=0;
    double den_sum=0;

    for(int i=0; i<N; i++) {

        std::uniform_int_distribution<int> uni(0,nAtoms);
        rand=uni(rng);
        if(rand<nAtoms) {
            for(int n=0; n<nAtoms; n++)
                PE=mcmove(Atoms,nAtoms);
        }
        else
      {
          PE=vmove(Atoms,nAtoms);
          density=nAtoms/vol;
      }
        P_sum+=vir_pre+pertail(r_cut);
        den_sum+=density;
        if(fmod(i,100)==0)
        {   cout<<i<<"\n";
            pair_correlation(Atoms,nAtoms,0,box,temp);
        }
        if(fmod(i,1000)==0);
        pair_correlation(Atoms,nAtoms,1,box,temp);
    }
    cout<<Press<<"\t"<<1.0/(3.0*vol)*P_sum/N+(den_sum/N)*temp<<"\t"<<den_sum/N<<"\n";
    //free(Atoms);
    cout<<"acceptance ratio\t"<<float(Nacc)/float(Iter)<<"\t"<<float(Nacc_v)/float(Iter_v)<<"\n";
    clock_t end =clock();
    double elapsed_time= double (end-begin)/CLOCKS_PER_SEC;
    cout<<"\n"<<elapsed_time;
    delete[] Atoms;
    return 0;
}


