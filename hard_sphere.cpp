using namespace std;
#include "MC_header.hpp"
#include "ini_pos.hpp"
#include "g_r.hpp"
#include "bop.hpp"
void print_pos (atom Atoms[],int nAtoms) {
    std::ofstream POSITION("OUT/config.dat");
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
    double lambda=0.2;
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
    //cout<<nn<<"\t"<<no<<"\n";
    if(move_accept(nn,no,nc)) //move_Accept determines if the move should be accepted or note depending on the bias potential.
    {
        no=nn; // if it is accepted then no is the new nn.
        if(flag) {
            HISTOGRAM[nn]++;
        }
    	delete[] old_Atoms;
        return no;
    }
    else
    {
        if(flag)
        {
            HISTOGRAM[no]++;
        }
        replace(Atoms,old_Atoms,nAtoms); //if the move is rejected then reset the config to what it was before the move.
    	delete[] old_Atoms;
        return no;
    }
}
int main(int argc,char* argv[]) {
    int nAtoms,N;
    long* HISTOGRAM;
    long n;
    double g_d;
//  cout<<"Enter the number of atoms:\n";  //comment out these too
//  cin>>nAtoms;				//these too
    std::ifstream infile(argv[1]);	//uncomment these lines to
    infile>>nAtoms;			//start the program
    infile>>box.x;			//from a inital
    box.z=box.y=box.x;			//configuraton
    atom* Atoms;
    Atoms = new (nothrow) atom[nAtoms];
    double a,b,c,d;				//uncomment the
    nAtoms=0;
    while(infile>>a>>b>>c>>d) { //>>e>>f>>g>>h) {	//start the pro
        Atoms[nAtoms].pos.x=b;			//from a inital
        Atoms[nAtoms].pos.y=c;			//configuraton
        Atoms[nAtoms].pos.z=d;			//
        nAtoms++;
    }						//
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
    cin>>Press;
    cin>>nc;
//  cin>>density;		//Comment out these line
//  vol=nAtoms/density;		//to change the program
//  box.x=pow(vol,1.0/3.0)/2.;	//to start from a given
//  box.y=pow(vol,1.0/3.0)/2.;	//initial parameters
//  box.z=pow(vol,1.0/3.0)/2.;	//
//  inipos(Atoms,nAtoms,box);	//
    density=nAtoms/(2*box.x*2*box.y*2*box.z);
    vol=2*box.x*2*box.y*2*box.z;
    cout<<box.x<<"\n";
    cout<<"density:"<<nAtoms/(2*box.x*2*box.y*2*box.z)<<"\n";
    cout<<check_overlap(Atoms,nAtoms)<<"\n";
    clock_t begin=clock();
    double EqN=3000;
    int END=0;
    int rand=0;
    cout<<"Press="<<Press<<"\n";
    char buffer[32];
    snprintf(buffer,sizeof(char)*32,"OUT/density_%f.dat",Press);
    std::ofstream DENSITY(buffer);
    cout<<nc<<"\n";
    snprintf(buffer,sizeof(char)*32,"OUT/cluster_%f.dat",float(nc));
    std::ofstream CS(buffer);
    n=largest_cluster(Atoms,nAtoms,l,box); // calculate the largest cluster in the initial config.
    cout<<"n="<<n<<"\n";
    int flag=0;
    snprintf(buffer,sizeof(char)*32,"OUT/Histogram_%f.dat",float(nc));
    for(int i=0; i<EqN; i++) {
        std::uniform_int_distribution<int> uni(0,nAtoms);
        rand=uni(rng);
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

        n=umbrella(Atoms,nAtoms,l,box,n,flag,HISTOGRAM);
        CS<<i<<"\t"<<n<<"\n";
        //    cout<<i<<"\n";
//        reset(Atoms,nAtoms);
        //    n=largest_cluster(Atoms,nAtoms,l,box);
        //      cout<<i<<"\t"<<n<<"\n";
        //      reset(Atoms,nAtoms);
        if(fmod(i,1000)==0)
        {   cout<<i<<"\n";
       //   std::ofstream HIS(buffer);
       //   for(int n=0; n<nAtoms; n++)
       //   {
       //       HIS<<n<<"\t"<<float(HISTOGRAM[n])/i<<"\n";
       //   }
       //   HIS.close();
        }
	if(n==nc){
		flag=1;
		break;	
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
        if(rand<nAtoms) {
            for(int n=0; n<nAtoms; n++)     //MC_Sweep
                mcmove_hardsphere(Atoms,nAtoms);
        }
        else {
            vmove(Atoms,nAtoms);
            density=nAtoms/vol;
            DENSITY<<i+EqN<<"\t"<<density<<"\n";
        }
        n=umbrella(Atoms,nAtoms,l,box,n,flag,HISTOGRAM);
        CS<<i+EqN<<"\t"<<n<<"\n";
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
    print_pos(Atoms,nAtoms);
    cout<<temp*k_b*(den_sum/N)*(1+2.*M_PI*(den_sum/N)*pow(R_CUT_HS,3)*g_d/3.0)<<"\t"<<den_sum/N<<"\t"<<g_d<<"\n";
    cout<<"acceptance ratio\t"<<float(Nacc)/float(Iter)<<"\t"<<float(Nacc_v)/float(Iter_v)<<"\n";
    clock_t end =clock();
    double elapsed_time= double (end-begin)/CLOCKS_PER_SEC;
    cout<<"\n"<<elapsed_time;
    delete[] Atoms;
    delete[] HISTOGRAM;

    return 0;
}
