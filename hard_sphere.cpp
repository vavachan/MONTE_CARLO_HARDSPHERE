using namespace std;
#include "MC_header.hpp"
#include "ini_pos.hpp"
#include "g_r.hpp"
#include "bop.hpp"
void print_pos (atom Atoms[],int nAtoms) {
    std::ofstream POSITION("config.dat");
    for (int i=0; i<nAtoms; i++)
        POSITION<<i<<"\t"<<Atoms[i].pos.x<<"\t"<<Atoms[i].pos.y<<"\t"<<Atoms[i].pos.z<<"\n";
    POSITION.close();
}

int hard_sphere(atom Atoms[],int nAtoms,int j) {
    Vector dr;
    double rr;
    for(int i=0; i<nAtoms; i++) {
        if(i!=j)
        {
            dr=v_sub(Atoms[i].pos,Atoms[j].pos);
            dr=VWrap(dr,box);
            rr=v_dot(dr,dr);
            if(rr<R_CUT_HS) {
                //cout<<i<<"\t"<<j<<"\t"<<rr<<"\n";
                return 0;
            }

        }
    }
    return 1;
}

void mcmove_hardsphere(atom Atoms[],int nAtoms) {
    double POTn,POTo,P,r;
    std::uniform_int_distribution<int> uni(0,nAtoms-1);
    auto random_integer = uni(rng);
    int flag=0;
    Vector random_dir;
    random_dir.x=uni_d(rng)-0.5;
    random_dir.y=uni_d(rng)-0.5;
    random_dir.z=uni_d(rng)-0.5;
    Atoms[random_integer].pos=v_Sadd(Atoms[random_integer].pos,random_dir,JUM);
    flag=hard_sphere(Atoms,nAtoms,random_integer);
    Iter++;
    if(flag)
    {
        Nacc++;
        return;
    }
    else
    {   Atoms[random_integer].pos=v_Sadd(Atoms[random_integer].pos,random_dir,-JUM);
        return;
    }

}
void vmove(atom Atoms[],int nAtoms,Vector box) {
    double old_vol,lnV0,lnV;
    double Lnew,f,r,P;
    old_vol=8*box.x*box.y*box.z;
    lnV0=log(old_vol);
    lnV=lnV0+dlnV*(uni_d(rng)-0.5);
    vol=exp(lnV);
    Lnew=pow(vol,0.3333333);
    f=Lnew/box.x;
    for (int n=0; n<nAtoms; n++) {
        Atoms[n].pos.x *= f;
        Atoms[n].pos.y *= f;
        Atoms[n].pos.z *= f;
    }
    P=exp(-1.0/temp*(Press*(old_vol-vol)-nAtoms*temp*log(vol/old_vol)));
    r=uni_d(rng);
    Iter_v++;
    if(r<P)
    {
        Nacc_v++;
        box.x *= f;
        box.y *= f;
        box.z *= f;
        return ;
    }
    else
    {
        f=1.0/f;
        for (int n=0; n<nAtoms; n++) {
            Atoms[n].pos.x *= f;
            Atoms[n].pos.y *= f;
            Atoms[n].pos.z *= f;
        }
        vol=old_vol;
        return ;
    }



}
void back_up(atom Atoms[],atom old_Atoms[],int nAtoms) {
    for(int n=0; n<nAtoms; n++) {
        old_Atoms[n].pos.x=Atoms[n].pos.x;
        old_Atoms[n].pos.y=Atoms[n].pos.y;
        old_Atoms[n].pos.z=Atoms[n].pos.z;
    }
}
void replace(atom Atoms[],atom old_Atoms[],int nAtoms) {
    for(int n=0; n<nAtoms; n++) {
        Atoms[n].pos.x=old_Atoms[n].pos.x;
        Atoms[n].pos.y=old_Atoms[n].pos.y;
        Atoms[n].pos.z=old_Atoms[n].pos.z;
    }
}
void reset(atom Atoms[],int nAtoms)
{
    for(int n=0; n<nAtoms; n++) {
        Atoms[n].cluster_index=-1;
        Atoms[n].connections=0;
    }


}
int move_accept(long nn,long no,long nc) {
    double P=0;
    double lambda=1.0;
    double r;
    P=exp(-1.0/temp*lambda*((nn-nc)*(nn-nc)-(no-nc)*(no-nc)));
    r=uni_d(rng);
    if(P>1)
    {
        return 1;
    }
    else if (r<P)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}
long umbrella(atom Atoms[],int nAtoms,int l,Vector box,long no,int flag,long HISTOGRAM[]) {
    long nn=0;
    atom* old_Atoms;
    old_Atoms = new (nothrow) atom[nAtoms];
    if(old_Atoms==nullptr) {
        cout<<"Memory Allocation Failed\n";
        return 0;
    }
    reset(Atoms,nAtoms);  // remove the cluster labels from the atoms
    back_up(Atoms,old_Atoms,nAtoms); //we need a copy of the config before the move.
    nn=largest_cluster(Atoms,nAtoms,l,box);  //calculate the new largest cluster in the system.
//      cout<<i<<"\t"<<nn<<"\t"<<no<<"\n";
    if(move_accept(nn,no,nc)) //move_Accept determines if the move should be accepted or note depending on the bias potential.
    {
        no=nn; // if it is accepted then no is the new nn.
        if(flag) {
            HISTOGRAM[nn]++;
        }
        return no;
    }
    else
    {
        if(flag)
        {
            HISTOGRAM[no]++;
        }
        replace(Atoms,old_Atoms,nAtoms); //if the move is rejected then reset the config to what it was before the move.
        return no;
    }
    delete[] old_Atoms;
}
int main() {
    int nAtoms,N;
    long* HISTOGRAM;
    long n;
    atom* Atoms;
    cout<<"Enter the number of atoms:\n";
    cin>>nAtoms;
    Atoms = new (nothrow) atom[nAtoms];
    HISTOGRAM = new (nothrow) long [nAtoms];
    if(Atoms==nullptr) {
        cout<<"Memory Allocation Failed\n";
        return 0;
    }
    for(int i=0; i<nAtoms; i++)
    {
        HISTOGRAM[i]=0;
    }
    cout<<"Enter how long you want to run the simulation:\n";
    cin>>N;
    cout<<"temp\n";
    cin>>temp;
    box=inipos_bcc(Atoms,nAtoms,box,1.3);
    density=nAtoms/(2*box.x*2*box.y*2*box.z);
    cout<<"density:"<<nAtoms/(2*box.x*2*box.y*2*box.z);
    clock_t begin=clock();
    double EqN=100;
    int END=0;
    for(int i=0; i<EqN; i++) {
        for(int n=0; n<nAtoms; n++)
            mcmove_hardsphere(Atoms,nAtoms);
        if(fmod(i,1000)==0)
        {   cout<<i<<"\n";
        }
    }
    cout<<"eq completed\n";
    Nacc=0;
    Iter=0;
    int flag=0;
    n=largest_cluster(Atoms,nAtoms,l,box); // calculate the largest cluster in the initial config.
    //cout<<n<<"\n";
    int rand=0;
    for(int i=0; i<N; i++) {
        for(int n=0; n<nAtoms; n++)     //MC_Sweep
            mcmove_hardsphere(Atoms,nAtoms);
        if(fmod(i,100)==0)
        {   cout<<i<<"\n";
            if(fmod(i,1000)==0)
            {
                pair_correlation(Atoms,nAtoms,1,box,temp);
            }
            else
                pair_correlation(Atoms,nAtoms,0,box,temp);
        }
        if(i>1000)
            flag=1;
    }
    for(int i=0; i<nAtoms; i++)
    {
    ;//    cout<<i<<"\t"<<HISTOGRAM[i]<<"\n";
    }
    print_pos(Atoms,nAtoms);
    cout<<"acceptance ratio\t"<<float(Nacc)/float(Iter)<<"\n";
    clock_t end =clock();
    double elapsed_time= double (end-begin)/CLOCKS_PER_SEC;
    cout<<"\n"<<elapsed_time;
    delete[] Atoms;
    delete[] HISTOGRAM;

    return 0;
}
