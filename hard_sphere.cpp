using namespace std;
#include "MC_header.hpp"
#include "ini_pos.hpp"
#include "g_r.hpp"
#include "bop.hpp"
#include <cstdlib>

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
int neigh_list_update(atom Atoms[],int nAtoms) {
    Vector dr;
    long double rr;
    long double twob=2*box.x;
//    cout<<"update\n";
    for(int n=0; n<nAtoms; n++) {
        Atoms[n].reset();
        Atoms[n].old_pos=Atoms[n].pos;
        Atoms[n].dist=0;
    }
    for(int i=0; i<nAtoms-1; i++) {
        for(int j=i+1; j<nAtoms; j++)
        {
            dr.x=Atoms[i].pos.x-Atoms[j].pos.x;
            dr.y=Atoms[i].pos.y-Atoms[j].pos.y;
            dr.z=Atoms[i].pos.z-Atoms[j].pos.z;
            dr.x=(dr.x-(twob*lround(dr.x/twob)));
            dr.y=(dr.y-(twob*lround(dr.y/twob)));
            dr.z=(dr.z-(twob*lround(dr.z/twob)));
            rr=dr.x*dr.x+dr.y*dr.y+dr.z*dr.z;
            //dr=v_sub(Atoms[i].pos,Atoms[j].pos);
            //dr=VWrap(dr,box);
            //rr=v_dot(dr,dr);
//	    cout<<i<<"\t"<<j<<"\t"<<rr<<"\n";
            if(rr<R_V_sq)
            {
                Atoms[i].neighbours+=1;
                Atoms[j].neighbours+=1;
                Atoms[i].update_neighbour(j);
                Atoms[j].update_neighbour(i);

            }

        }
    }
}

int verlet_list_HS(atom Atoms[],int nAtoms,int random_integer) {
    Vector dr;
    long double twob=2*box.x;
    long double rr;
    // cout<<"verlet\n";
    for(int i=0; i<Atoms[random_integer].neighbours; i++) {
        dr.x=Atoms[random_integer].pos.x-Atoms[Atoms[random_integer].neigh_list[i]].pos.x;
        dr.y=Atoms[random_integer].pos.y-Atoms[Atoms[random_integer].neigh_list[i]].pos.y;
        dr.z=Atoms[random_integer].pos.z-Atoms[Atoms[random_integer].neigh_list[i]].pos.z;
        dr.x=(dr.x-(twob*lround(dr.x/twob)));
        dr.y=(dr.y-(twob*lround(dr.y/twob)));
        dr.z=(dr.z-(twob*lround(dr.z/twob)));
        rr=dr.x*dr.x+dr.y*dr.y+dr.z*dr.z;
        //  dr=v_sub(Atoms[random_integer].pos,Atoms[Atoms[random_integer].neigh_list[i]].pos);
        //  dr=VWrap(dr,box);
        //  rr=v_dot(dr,dr);
        //          cout<<random_integer<<"\t"<<Atoms[random_integer].neigh_list[i]<<"\t"<<rr<<"\n";
        if(rr<R_CUT_HS_sq) {
            return 0;
        }
    }
    return 1;
}

int hard_sphere(atom Atoms[],int nAtoms,int j) {
    Vector dr;
    long double rr;
    long double twob=2*box.x;
    //cout<<"hard\n";
    for(int i=0; i<nAtoms; i++) {
        if(i!=j)
        {

            dr.x=Atoms[i].pos.x-Atoms[j].pos.x;
            dr.y=Atoms[i].pos.y-Atoms[j].pos.y;
            dr.z=Atoms[i].pos.z-Atoms[j].pos.z;
            dr.x=(dr.x-(twob*lround(dr.x/twob)));
            dr.y=(dr.y-(twob*lround(dr.y/twob)));
            dr.z=(dr.z-(twob*lround(dr.z/twob)));
            rr=dr.x*dr.x+dr.y*dr.y+dr.z*dr.z;
            //dr=v_sub(Atoms[i].pos,Atoms[j].pos);
            //dr=VWrap(dr,box);
            //rr=v_dot(dr,dr);
            if(rr<R_CUT_HS) {
//               cout<<j<<"\t"<<i<<"\t"<<rr<<"\n\n";
                return 0;
            }

        }
    }
    return 1;
}
int check_overlap(atom Atoms[],int nAtoms)
{
    Vector dr;
    long double rr,r;
    long double twob=2*box.x;
    for(int i=0; i<nAtoms; i++) {

        for(int j=0; j<Atoms[i].neighbours; j++) {
            dr.x=Atoms[i].pos.x-Atoms[Atoms[i].neigh_list[j]].pos.x;
            dr.y=Atoms[i].pos.y-Atoms[Atoms[i].neigh_list[j]].pos.y;
            dr.z=Atoms[i].pos.z-Atoms[Atoms[i].neigh_list[j]].pos.z;
            dr.x=(dr.x-(twob*lround(dr.x/twob)));
            dr.y=(dr.y-(twob*lround(dr.y/twob)));
            dr.z=(dr.z-(twob*lround(dr.z/twob)));
            rr=dr.x*dr.x+dr.y*dr.y+dr.z*dr.z;
            //  dr=v_sub(Atoms[random_integer].pos,Atoms[Atoms[random_integer].neigh_list[i]].pos);
            //  dr=VWrap(dr,box);
            //  rr=v_dot(dr,dr);
            //          cout<<random_integer<<"\t"<<Atoms[random_integer].neigh_list[i]<<"\t"<<rr<<"\n";
            if(rr<R_CUT_HS_sq) {
//		    cout<<i<<"\t"<<Atoms[i].neigh_list[j]<<"\t"<<rr<<"\n";
                return 0;
            }
        }
        //  for(int j=i+1; j<nAtoms; j++)
        //  {
        //      dr=v_sub(Atoms[i].pos,Atoms[j].pos);
        //      dr=VWrap(dr,box);
        //      rr=v_dot(dr,dr);
        //      r=sqrt(rr);
        //      if(r<R_CUT_HS) {
        //          return 0;
        //      }

        //  }

    }
    return 1;
}

void mcmove_hardsphere(atom Atoms[],int nAtoms) {
    Vector dr;
    long double rr;
    long double POTn,POTo,P,r;
    std::uniform_int_distribution<int> uni(0,nAtoms-1);
    auto random_integer = uni(rng);
    int flag=0;
    Vector random_dir;
    atom old_atom;
    old_atom=Atoms[random_integer];
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
    {   //   Atoms[random_integer].pos=v_Sadd(Atoms[random_integer].pos,random_dir,-JUM);
        Atoms[random_integer]=old_atom;
        return;
    }
}
void mcmove_ver_hardsphere(atom Atoms[],int nAtoms) {
    Vector dr;
    long double rr;
    long double POTn,POTo,P,r;
    std::uniform_int_distribution<int> uni(0,nAtoms-1);
    auto random_integer = uni(rng);
    long double twob=2*box.x;
    int flag=0;
    Vector random_dir;
    atom old_atom;
    old_atom=Atoms[random_integer];
    if(Atoms[random_integer].dist>gap)
    {
        nei_up+=1;
        cout<<"anybody\n";
        neigh_list_update(Atoms,nAtoms) ;
    }
    random_dir.x=uni_d(rng)-0.5;
    random_dir.y=uni_d(rng)-0.5;
    random_dir.z=uni_d(rng)-0.5;

    Atoms[random_integer].pos=v_Sadd(Atoms[random_integer].pos,random_dir,JUM);

    dr.x=Atoms[random_integer].pos.x-old_atom.pos.x;
    dr.y=Atoms[random_integer].pos.y-old_atom.pos.y;
    dr.z=Atoms[random_integer].pos.z-old_atom.pos.z;
    dr.x=(dr.x-(twob*lround(dr.x/twob)));
    dr.y=(dr.y-(twob*lround(dr.y/twob)));
    dr.z=(dr.z-(twob*lround(dr.z/twob)));
    rr=dr.x*dr.x+dr.y*dr.y+dr.z*dr.z;
    r=sqrt(rr);

    Atoms[random_integer].dist+=r;

    if(Atoms[random_integer].dist>gap)
    {
        nei_up+=1;
        neigh_list_update(Atoms,nAtoms) ;
    }
    flag=verlet_list_HS(Atoms,nAtoms,random_integer);
    Iter++;
    if(flag)
    {
        Nacc++;
        return;
    }
    else
    {
        Atoms[random_integer]=old_atom;
        return;
    }
}
void back_up(atom Atoms[],atom old_Atoms[],int nAtoms) {
    for(int n=0; n<nAtoms; n++) {
	    old_Atoms[n]=Atoms[n];
     // old_Atoms[n].pos.x=Atoms[n].pos.x;
     // old_Atoms[n].pos.y=Atoms[n].pos.y;
     // old_Atoms[n].pos.z=Atoms[n].pos.z;
    }
}
void replace(atom Atoms[],atom old_Atoms[],int nAtoms) {
    for(int n=0; n<nAtoms; n++) {
	    Atoms[n]=old_Atoms[n];
    //  Atoms[n].pos.x=old_Atoms[n].pos.x;
    //  Atoms[n].pos.y=old_Atoms[n].pos.y;
    //  Atoms[n].pos.z=old_Atoms[n].pos.z;
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
void close_reset(atom Atoms[],int nAtoms)
{
    for(int n=0; n<nAtoms; n++) {
        Atoms[n].cluster_index=-1;
        Atoms[n].connections=0;
        Atoms[n].close_reset();
    }
}
void vmove(atom Atoms[],int nAtoms) {
    long double old_vol,lnV0,lnV;
    long double Lnew,f,r,P;
    long double PE_new;
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
int move_accept(long nn,long no,long nc,int flag) {
    long double P=0;
    long double lambda=0.15;
    long double r;
// if(flag) {
//     if(abs(nn-nc)>5)
//         return 0;
//     else
//         return 1;
// }
//   else {
    P=exp(-1.0/temp*(lambda/2.0)*((nn-nc)*(nn-nc)-(no-nc)*(no-nc)));
    r=uni_d(rng);
//    cout<<r<<"\t"<<P<<"\n";
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
    // }
}
long umbrella(atom Atoms[],atom old_Atoms[],int nAtoms,int l,Vector box,long no,int flag,long HISTOGRAM[]) {
    long nn=0;
    close_reset(Atoms,nAtoms);  // remove the cluster labels from the atoms

//   reset(Atoms,nAtoms);  // remove the cluster labels from the atoms
    nn=largest_cluster(Atoms,nAtoms,l,box);  //calculate the new largest cluster in the system.
    cout<<nn<<"\t"<<no<<"\n";
    if(move_accept(nn,no,nc,flag)) //move_Accept determines if the move should be accepted or note depending on the bias potential.
    {
        no=nn; // if it is accepted then no is the new nn.
        return no;
    }
    else
    {
        replace(Atoms,old_Atoms,nAtoms); //if the move is rejected then reset the config to what it was before the move.
	box=old_box;
        return no;
    }
}
int main(int argc,char* argv[]) {
    int nAtoms,N;
    long* HISTOGRAM;
    long n;
    long double g_d;
    atom* Atoms;
    bool restart;
    bool bias;
//##############################################################################################################################
    std::string action(argv[1]);
    if (action == "restart") {
        restart= true;
    } else if (action == "new") {
        restart= false;
    } else {
        cout<<"invalid argument\n\n\nexiting....";
        return 0;
    }
    std::string action2(argv[2]);
    if (action2 == "bias") {
        bias= true;
    } else if (action2 == "no-bias") {
        bias= false;
    } else {
        cout<<"invalid argument\n\n\nexiting....";
        return 0;
    }
//##############################################################################################################################
    if(!restart) {
        cin>>nAtoms;				//these too
        Atoms = new (nothrow) atom[nAtoms];
        cin>>density;		//Comm
        vol=nAtoms/density;
        box.x=pow(vol,1.0/3.0)/2.;
        box.y=pow(vol,1.0/3.0)/2.;
        box.z=pow(vol,1.0/3.0)/2.;
        inipos(Atoms,nAtoms,box);
        //  random_ini(Atoms,nAtoms,box);

    }
    else
    {   std::ifstream infile(argv[3]);
        infile>>nAtoms;
        infile>>box.x;
        box.z=box.y=box.x;
        Atoms = new (nothrow) atom[nAtoms];
        long double a,b,c,d;				//uncomment the
        nAtoms=0;
        while(infile>>a>>b>>c>>d) { //>>e>>f>>g>>h)
            Atoms[nAtoms].pos.x=b;
            Atoms[nAtoms].pos.y=c;
            Atoms[nAtoms].pos.z=d;
            nAtoms++;
        }
        density=nAtoms/(2*box.x*2*box.y*2*box.z);
    }
//##############################################################################################################################
    neigh_list_update(Atoms,nAtoms);
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
    long double EqN;
    cin>>EqN;
    vol=2*box.x*2*box.y*2*box.z;
    char buffer[64];
    snprintf(buffer,sizeof(char)*64,"OUT/out_%d_%d_%f.dat",int(nAtoms),int(nc),Press);//_%d_%f.dat",int(nAtoms),Press);
    //freopen(buffer,"w",stdout);
    cout<<"no of Atoms:"<<nAtoms<<"\n"<<flush;
    cout<<"N:"<<N<<"\n"<<flush;
    cout<<"EqN:"<<EqN<<"\n"<<flush;
    cout<<"Pressure:"<<Press<<"\n"<<flush;
    cout<<"temp:"<<temp<<"\n"<<flush;
    cout<<"box:"<<box.x<<"\n"<<flush;
    cout<<"density:"<<nAtoms/(2*box.x*2*box.y*2*box.z)<<"\n"<<flush;
    cout<<"overlap:"<<check_overlap(Atoms,nAtoms)<<"\n"<<flush;
    int END=0;
    int rand=0;
    snprintf(buffer,sizeof(char)*64,"OUT/density_%d_%f.dat",int(nAtoms),Press);
    std::ofstream DENSITY(buffer);
    cout<<"nc:"<<nc<<"\n"<<flush;
    snprintf(buffer,sizeof(char)*64,"OUT/cluster_%d_%f.dat",int(nc),Press);
    std::ofstream CS(buffer);
    int flag=0;
    snprintf(buffer,sizeof(char)*64,"OUT/Histogram_%d_%d_%f.dat",int(nAtoms),int(nc),Press);
    atom* old_Atoms;
    old_Atoms = new (nothrow) atom[nAtoms];
    int break_point=0;
    if(old_Atoms==nullptr) {
        cout<<"Memory Allocation Failed\n";
        return 0;
    }
    nei_up=0;
    int N_iter=0;
    close_reset(Atoms,nAtoms);
    clock_t begin=clock();
    n=largest_cluster(Atoms,nAtoms,l,box); // calculate the largest cluster in the initial config.
    cout<<"n="<<n<<"\n"<<flush;
    int n1;
//#####################################################################################################################################3
    for(int i=0; i<EqN; i++) {
////////    if(!check_overlap(Atoms,nAtoms))
////////		  cout<<i<<"\t"<<"whoops\n";
//	cout<<n<<"\n";
     // close_reset(Atoms,nAtoms);
     // n1=largest_cluster(Atoms,nAtoms,l,box);
     // cout<<"\t"<<n1<<"\n";
    //  close_reset(old_Atoms,nAtoms);
    //  n1=largest_cluster(old_Atoms,nAtoms,l,box);
    //  cout<<n1<<"\n";
    //  close_reset(Atoms,nAtoms);
    //  n=largest_cluster(Atoms,nAtoms,l,box);
    //  cout<<n<<"\n";
        if(fmod(i,5)==0) {
            if(Nacc*(1.0/Iter)<0.5)
            {
                JUM=JUM*0.95;
            }
            else
            {
                JUM=JUM*1.05;
            }
            // cout<<i<<"\t"<<JUM<<"\n";
            Nacc=0;
            Iter=0;
        }
        if(fmod(i,nAtoms)==0) {
            if(Nacc_v*(1.0/Iter_v)<0.5)
            {
                dlnV=dlnV*0.95;
            }
            else
            {
                dlnV=dlnV*1.05;
            }
        //    cout<<Nacc_v<<"\t"<<Iter_v<<"\t"<<Nacc_v*(1.0/Iter_v)<<"\n";
            Nacc_v=0;
            Iter_v=0;

        }
        for(int n=0; n<nAtoms; n++)     //MC_Sweep
        {
            std::uniform_int_distribution<int> uni(0,nAtoms);
            rand=uni(rng);
            if(rand<nAtoms) {
                {
//              mcmove_hardsphere(Atoms,nAtoms);
                    mcmove_ver_hardsphere(Atoms,nAtoms);
                    N_iter++;
                }
            }
            else {
                vmove(Atoms,nAtoms);
                density=nAtoms/vol;
                // cout<<box.x<<"\n";
                DENSITY<<i<<"\t"<<density<<"\n"<<flush;
            }
        }
    //  for(int n=0;n<nAtoms;n++)
    //  {
    //  	;
    //  	cout<<n<<"\t"<<old_Atoms[n].pos.x-Atoms[n].pos.x<<"\t"<<old_Atoms[n].pos.y-Atoms[n].pos.y<<"\t"<<old_Atoms[n].pos.z-Atoms[n].pos.z<<"\n";//<<
    // // cout<<n<<"\t"<<Atoms[n].pos.x<<"\t"<<Atoms[n].pos.y<<"\t"<<Atoms[n].pos.z<<"\n";//<<
    //  }
        if(bias)// and (fmod(i,20)==0))
        {

            n=umbrella(Atoms,old_Atoms,nAtoms,l,box,n,flag,HISTOGRAM);
            back_up(Atoms,old_Atoms,nAtoms); //we need a copy of the config before the move.
	    old_box=box;
            // HISTOGRAM[n]+=1;
//	    if(fmod(i,1000)==0)
            CS<<i<<"\t"<<n<<"\n"<<flush;
        }

        //    cout<<i<<"\n";
        //  cout<<i<<"\t"<<n<<"\n";
        //      reset(Atoms,nAtoms);

        if(fmod(i,1000)==0)
        {   cout<<i*1.0/EqN<<"\n"<<flush;
            // if(bias) {
            //     std::ofstream HIS(buffer);
            //     for(int n=0; n<nAtoms; n++)
            //     {
            //         HIS<<n<<"\t"<<float(HISTOGRAM[n])/i<<"\n"<<flush;
            //     }
            //     HIS.close();
            //   }
            print_pos(Atoms,nAtoms);
        }
        if(bias) {
            if(n==nc) {
                flag=1;
                break_point=i;
                cout<<"yay\n";
                break;
            }
        }
    }
    if(!bias)
        break_point=EqN;
    cout<<"eq completed\n";
    cout<<nei_up*1.0/(N_iter)<<"\n";
    Nacc=0;
    Iter=0;
    Nacc_v=0;
    Iter_v=0;
    long double den_sum=0;
    for(int i=0; i<N; i++) {
////////if(bias and (fmod(i,10)==0)) {
////////    close_reset(old_Atoms,nAtoms);
////////    n=largest_cluster(old_Atoms,nAtoms,l,box);
////////    cout<<n<<"\n";
////////}
        if(fmod(i,5)==0) {
            if(Nacc*(1.0/Iter)<0.5)
            {
                JUM=JUM*0.95;
            }
            else
            {
                JUM=JUM*1.05;
            }
            // cout<<i<<"\t"<<JUM<<"\n";
            Nacc=0;
            Iter=0;
        }
        if(fmod(i,5*nAtoms)==0) {
            if(Nacc_v*(1.0/Iter_v)<0.5)
            {
                dlnV=dlnV*0.95;
            }
            else
            {
                dlnV=dlnV*1.05;
            }
            Nacc_v=0;
            Iter_v=0;
        }
        for(int n=0; n<nAtoms; n++)     //MC_Sweep
        {
            std::uniform_int_distribution<int> uni(0,nAtoms);
            rand=uni(rng);
            if(rand<nAtoms) {
                {
//              mcmove_hardsphere(Atoms,nAtoms);
                    mcmove_ver_hardsphere(Atoms,nAtoms);
                    N_iter++;
                }
            }
            else {
                vmove(Atoms,nAtoms);
                density=nAtoms/vol;
                // cout<<box.x<<"\n";
                DENSITY<<i+break_point<<"\t"<<density<<"\n"<<flush;
            }
        }
        if(bias and (fmod(i,10)==0)) {
            n=umbrella(Atoms,old_Atoms,nAtoms,l,box,n,flag,HISTOGRAM);
            back_up(Atoms,old_Atoms,nAtoms); //we need a copy of the config before the move.
	    old_box=box;
            if(flag)
                HISTOGRAM[n]++;
//	    if(fmod(i,1000)==0)
            CS<<i+break_point<<"\t"<<n<<"\n"<<flush;
        }
        den_sum+=density;
        if(fmod(i,100)==0)
        {
            if(fmod(i,1000)==0)
            {   cout<<i*1.0/N<<"\n"<<flush;
                g_d=pair_correlation(Atoms,nAtoms,1,box,temp,Press);
                std::ofstream HIS(buffer);
                for(int n=0; n<nAtoms; n++)
                {
                    HIS<<n<<"\t"<<float(HISTOGRAM[n])<<"\n"<<flush;
                }
                HIS.close();
                print_pos(Atoms,nAtoms);
            }
            else
                g_d=pair_correlation(Atoms,nAtoms,0,box,temp,Press);
//      	cout<<i<<"\t"<<g_d<<"\n";
        }
    }
    cout<<temp*k_b*(den_sum/N)*(1+2.*M_PI*(den_sum/N)*pow(R_CUT_HS,3)*g_d/3.0)<<"\t"<<den_sum/N<<"\t"<<(M_PI/6.0)*den_sum/N<<"\t"<<g_d<<"\n";
    cout<<"acceptance ratio\t"<<float(Nacc)/float(Iter)<<"\t"<<float(Nacc_v)/float(Iter_v)<<"\n";
    clock_t end =clock();
    double elapsed_time= double (end-begin)/CLOCKS_PER_SEC;
    cout<<elapsed_time/60.<<"\n";
    delete[] Atoms;
    delete[] HISTOGRAM;
    delete[] old_Atoms;

    return 0;
}
