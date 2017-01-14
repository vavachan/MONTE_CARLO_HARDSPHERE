using namespace std;
#include "MC_header.hpp"
#include "ini_pos.hpp"
#include "g_r.hpp"
#include "bop.hpp"
void print_pos (atom Atoms[],int nAtoms) {
    char buffer[64];
    snprintf(buffer,sizeof(char)*64,"OUT/config_%f_%f_%f.dat",temp,float(nAtoms),Press);
    std::ofstream POSITION(buffer);
    POSITION<<nAtoms<<"\n";
    POSITION<<box.x<<"\n";
    Vector dr;
    for (int i=0; i<nAtoms; i++)
    {
        dr=VWrap(Atoms[i].pos,box);
        POSITION<<R_CUT_HS/2.0<<"\t"<<dr.x<<"\t"<<dr.y<<"\t"<<dr.z<<"\n";
    }
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
int check_overlap(atom Atoms[],int nAtoms)
{

    Vector dr;
    double rr,r;
    for(int i=0; i<nAtoms-1; i++) {
        for(int j=i+1; j<nAtoms; j++)
        {
            dr=v_sub(Atoms[i].pos,Atoms[j].pos);
            dr=VWrap(dr,box);
            rr=v_dot(dr,dr);
            r=sqrt(rr);
            if(r<R_CUT_HS) {
//                cout<<i<<"\t"<<j<<"\t"<<rr<<"\n";
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
        Atoms[n].reset();
    }


}
void vmove(atom Atoms[],int nAtoms) {
    double old_vol,lnV0,lnV;
    double Lnew,f,r,P;
    double PE_new;
    int flag=0;
    old_vol=vol;
    lnV0=log(old_vol);
    lnV=lnV0+dlnV*(uni_d(rng)-0.5);
    vol=exp(lnV);
    Lnew=pow(vol,1.0/3.0)/2.0;
    f=Lnew/box.x;
    box.x *= f;
    box.y *= f;
    box.z *= f;
    for (int n=0; n<nAtoms; n++) {
        Atoms[n].pos.x *= f;
        Atoms[n].pos.y *= f;
        Atoms[n].pos.z *= f;
    }
    flag=check_overlap(Atoms,nAtoms);
    P=exp(-1.0/(k_b*temp)*(Press*(vol-old_vol)-(nAtoms+1)*k_b*temp*log(vol/old_vol)));
    r=uni_d(rng);
    Iter_v++;
    if(r<P and flag)
    {
        Nacc_v++;
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

        box.x *= f;
        box.y *= f;
        box.z *= f;
        vol=old_vol;
        return ;
    }
}
int move_accept(long nn,long no,long nc) {
    double P=0;
    double lambda=0.1;
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
long umbrella(atom Atoms[],atom old_Atoms[],int nAtoms,int l,Vector box,long no,int flag,long HISTOGRAM[]) {
    long nn=0;
    reset(Atoms,nAtoms);  // remove the cluster labels from the atoms

    nn=largest_cluster(Atoms,nAtoms,l,box);  //calculate the new largest cluster in the system.
    // cout<<nn<<"\t"<<no<<"\n";
    if(move_accept(nn,no,nc)) //move_Accept determines if the move should be accepted or note depending on the bias potential.
    {
        no=nn; // if it is accepted then no is the new nn.
        // if(flag) {
        //     HISTOGRAM[nn]++;
        // }
        // return no;
    }
    else
    {
        //  if(flag)
        //  {
        //      HISTOGRAM[no]++;
        //  }
        replace(Atoms,old_Atoms,nAtoms); //if the move is rejected then reset the config to what it was before the move.
        return no;
    }
}
int main(int argc,char* argv[]) {
    int nAtoms,N;
    long* HISTOGRAM;
    long n;
    double g_d;
    atom* Atoms;
    bool restart;
    std::string action(argv[1]);
    if (action == "restart") {
        restart= true;
    } else if (action == "new") {
        restart= false;
    } else {
        cout<<"invalid argument\n\n\nexiting....";
        return 0;
    }
    bool bias;
    std::string action2(argv[2]);
    if (action2 == "bias") {
        bias= true;
    } else if (action2 == "no-bias") {
        bias= false;
    } else {
        cout<<"invalid argument\n\n\nexiting....";
        return 0;
    }
    if(!restart) {
        cin>>nAtoms;				//these too
        Atoms = new (nothrow) atom[nAtoms];
        cin>>density;		//Comment out these line
        vol=nAtoms/density;		//to change the program
        box.x=pow(vol,1.0/3.0)/2.;	//to start from a given
        box.y=pow(vol,1.0/3.0)/2.;	//initial parameters
        box.z=pow(vol,1.0/3.0)/2.;	//
        inipos(Atoms,nAtoms,box);	//
    }
    else
    {   std::ifstream infile(argv[3]);	//uncomment these lines to
        infile>>nAtoms;			//start the program
        infile>>box.x;			//from a inital
        box.z=box.y=box.x;
        Atoms = new (nothrow) atom[nAtoms];
        double a,b,c,d;				//uncomment the
        nAtoms=0;
        while(infile>>a>>b>>c>>d) { //>>e>>f>>g>>h) {	//start the pro
            Atoms[nAtoms].pos.x=b;			//from a inital
            Atoms[nAtoms].pos.y=c;			//configuraton
            Atoms[nAtoms].pos.z=d;			//
            nAtoms++;
        }						//
        density=nAtoms/(2*box.x*2*box.y*2*box.z);
    }
    HISTOGRAM = new (nothrow) long [nAtoms];
    if(Atoms==nullptr) {
        cout<<"Memory Allocation Failed\n";
        return 0;
    }
    for(int i=0; i<nAtoms; i++)
    {
        HISTOGRAM[i]=0;
    }
    cin>>N;
    cin>>temp;
    cin>>Press;
    cin>>nc;
    double EqN;
    cin>>EqN;
    vol=2*box.x*2*box.y*2*box.z;
    cout<<"no of Atoms:"<<nAtoms<<"\n";
    cout<<"N:"<<N<<"\n";
    cout<<"Pressure:"<<Press<<"\n";
    cout<<"temp:"<<temp<<"\n";
    cout<<"box:"<<box.x<<"\n";
    cout<<"density:"<<nAtoms/(2*box.x*2*box.y*2*box.z)<<"\n";
    cout<<"overlap:"<<check_overlap(Atoms,nAtoms)<<"\n";
    clock_t begin=clock();
    int END=0;
    int rand=0;
    char buffer[64];
    snprintf(buffer,sizeof(char)*32,"OUT/density_%f.dat",Press);
    std::ofstream DENSITY(buffer);
    cout<<"nc:"<<nc<<"\n";
    snprintf(buffer,sizeof(char)*32,"OUT/cluster_%f.dat",float(nc));
    std::ofstream CS(buffer);
    n=largest_cluster(Atoms,nAtoms,l,box); // calculate the largest cluster in the initial config.
    cout<<"n="<<n<<"\n";
    int flag=0;
    snprintf(buffer,sizeof(char)*64,"OUT/Histogram_%f_%f.dat",float(nc),Press);
    atom* old_Atoms;
    old_Atoms = new (nothrow) atom[nAtoms];
    int break_point=0;
    if(old_Atoms==nullptr) {
        cout<<"Memory Allocation Failed\n";
        return 0;
    }
    for(int i=0; i<EqN; i++) {
        std::uniform_int_distribution<int> uni(0,nAtoms);
        rand=uni(rng);
        back_up(Atoms,old_Atoms,nAtoms); //we need a copy of the config before the move.
        if(rand<nAtoms) {
            for(int n=0; n<nAtoms; n++)     //MC_Sweep
                mcmove_hardsphere(Atoms,nAtoms);
        }
        else {
            vmove(Atoms,nAtoms);
            density=nAtoms/vol;
            // cout<<box.x<<"\n";
            DENSITY<<i<<"\t"<<density<<"\n";
        }
        if(bias)
        {
            n=umbrella(Atoms,old_Atoms,nAtoms,l,box,n,flag,HISTOGRAM);
            CS<<i<<"\t"<<n<<"\n";
        }
        //    cout<<i<<"\n";
        //    reset(Atoms,nAtoms);
        //    n=largest_cluster(Atoms,nAtoms,l,box);
        //      cout<<i<<"\t"<<n<<"\n";
        //      reset(Atoms,nAtoms);
        if(fmod(i,1000)==0)
        {   cout<<i<<"\n";
            if(bias) {
                std::ofstream HIS(buffer);
                for(int n=0; n<nAtoms; n++)
                {
                    HIS<<n<<"\t"<<float(HISTOGRAM[n])/i<<"\n";
                }
                HIS.close();
            }
            print_pos(Atoms,nAtoms);
        }
        if(bias) {
            if(n==nc) {
                flag=1;
                break_point=i;
                break;
            }
        }
    }
    cout<<"eq completed\n";
    Nacc=0;
    Iter=0;
    Nacc_v=0;
    Iter_v=0;
    double den_sum=0;
    for(int i=0; i<N; i++) {
        std::uniform_int_distribution<int> uni(0,nAtoms);
        rand=uni(rng);
//   	back_up(Atoms,old_Atoms,nAtoms); //we need a copy of the config before the move.
        if(rand<nAtoms) {
            for(int n=0; n<nAtoms; n++)     //MC_Sweep
                mcmove_hardsphere(Atoms,nAtoms);
        }
        else {
            vmove(Atoms,nAtoms);
            density=nAtoms/vol;
            DENSITY<<i+break_point<<"\t"<<density<<"\n";
        }
        if(bias) {
            n=umbrella(Atoms,old_Atoms,nAtoms,l,box,n,flag,HISTOGRAM);
            if(flag)
                HISTOGRAM[n]++;
            CS<<i+break_point<<"\t"<<n<<"\n";
        }
        den_sum+=density;
        if(fmod(i,100)==0)
        {
            if(fmod(i,1000)==0)
            {   cout<<i<<"\n";
                g_d=pair_correlation(Atoms,nAtoms,1,box,temp);
                std::ofstream HIS(buffer);
                for(int n=0; n<nAtoms; n++)
                {
                    HIS<<n<<"\t"<<float(HISTOGRAM[n])<<"\n";
                }
                HIS.close();
            }
            else
                g_d=pair_correlation(Atoms,nAtoms,0,box,temp);
//      	cout<<i<<"\t"<<g_d<<"\n";
        }
    }
    cout<<temp*k_b*(den_sum/N)*(1+2.*M_PI*(den_sum/N)*pow(R_CUT_HS,3)*g_d/3.0)<<"\t"<<den_sum/N<<"\t"<<g_d<<"\n";
    cout<<"acceptance ratio\t"<<float(Nacc)/float(Iter)<<"\t"<<float(Nacc_v)/float(Iter_v)<<"\n";
    clock_t end =clock();
    double elapsed_time= double (end-begin)/CLOCKS_PER_SEC;
    cout<<"\n"<<elapsed_time/60.;
    delete[] Atoms;
    delete[] HISTOGRAM;
    delete[] old_Atoms;

    return 0;
}
