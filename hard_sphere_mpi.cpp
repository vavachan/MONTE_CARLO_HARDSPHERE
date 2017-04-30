using namespace std;
#include "MC_header.hpp"
#include "ini_pos.hpp"
#include "g_r.hpp"
#include "bop.hpp"
#include <cstdlib>
#include <string.h>
#include <iomanip>
#include <mpi.h>
void print_pos (atom Atoms[],int nAtoms,int i,char label[],int n) {
    char buffer[64];
    snprintf(buffer,sizeof(char)*64,"OUT/configs/n_config/config_%.2f_%d_%.2f_%d_%s.dat",temp,nAtoms,Press,int(n),label);
    std::ofstream POSITION(buffer,std::ios_base::app);
    snprintf(buffer,sizeof(char)*64,"OUT/configs/n_config/restart_config_%.2f_%d_%.2f_%d_%s.dat",temp,nAtoms,Press,int(nc),label);
    std::ofstream R_POSITION(buffer);
    POSITION<<nAtoms<<"\n";
    POSITION<<box.x<<"\n";
    POSITION<<i<<"\n";
    R_POSITION<<nAtoms<<"\n";
    R_POSITION<<box.x<<"\n";
    R_POSITION<<i<<"\n";
    Vector dr;
    for (int i=0; i<nAtoms; i++)
    {
        dr=VWrap(Atoms[i].pos,box);
        POSITION<<R_CUT_HS/2.0<<setprecision(15)<<"\t"<<dr.x<<"\t"<<dr.y<<"\t"<<dr.z<<"\n";
        R_POSITION<<R_CUT_HS/2.0<<setprecision(15)<<"\t"<<dr.x<<"\t"<<dr.y<<"\t"<<dr.z<<"\n";
    }
    POSITION.close();
    R_POSITION.close();

}
int neigh_list_update(atom Atoms[],int nAtoms) {
    Vector dr;
    long double rr;
    long double twob=2*box.x;
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
int check_overlap_mc(atom Atoms[],int nAtoms)
{
    Vector dr;
    long double rr,r;
    long double twob=2*box.x;
    for(int i=0; i<nAtoms; i++) {
        for(int j=i+1; j<nAtoms; j++)
        {
            dr=v_sub(Atoms[i].pos,Atoms[j].pos);
            dr=VWrap(dr,box);
            rr=v_dot(dr,dr);
            r=sqrt(rr);
            if(r<R_CUT_HS) {
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
            if(rr<R_CUT_HS_sq) {
//		cout<<i<<"\t"<<Atoms[i].neigh_list[j]<<"\t"<<rr<<"\n";
                return 0;
            }
        }

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
    //cout<<nn<<"\t"<<no<<"\t"<<P<<"\n";
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
    // }
}
long umbrella(atom Atoms[],atom old_Atoms[],int nAtoms,int l,Vector box,long no,int flag,long HISTOGRAM[],char label[],long nc_temp) {
    long nn=0;
    close_reset(Atoms,nAtoms);  // remove the cluster labels from the atoms

//   reset(Atoms,nAtoms);  // remove the cluster labels from the atoms
    nn=largest_cluster(Atoms,nAtoms,l,box,label);  //calculate the new largest cluster in the system.
    //cout<<nn<<"\t"<<no<<"\t";
    if(move_accept(nn,no,nc_temp,flag)) //move_Accept determines if the move should be accepted or note depending on the bias potential.
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
    long tempnc;
    long double g_d;
    atom* Atoms;
    bool restart;
    bool bias;
    int START=0;
    int index=0;
    char label[20];
    bool equilibriate=true;
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
        //   inipos(Atoms,nAtoms,box);
        //  random_ini(Atoms,nAtoms,box,R_CUT_HS);
        box=FCC(Atoms,nAtoms,box,density);
    }
    else
    {   std::ifstream infile(argv[3]);

        infile>>nAtoms;
        infile>>box.x;
        infile>>START;
        box.z=box.y=box.x;
        Atoms = new (nothrow) atom[nAtoms];
        long double a,b,c,d;				//uncomment the
        nAtoms=0;
        while(infile>>a>>b>>c>>d) { //>>e>>f>>g>>h)
            Atoms[nAtoms].pos.x=b;
            Atoms[nAtoms].pos.y=c;
            Atoms[nAtoms].pos.z=d;
            //    cout<<nAtoms<<"\t"<<Atoms[nAtoms].pos.x<<"\t"<<Atoms[nAtoms].pos.y<<"\t"<<Atoms[nAtoms].pos.z<<"\n";
            nAtoms++;
        }
        cout<<"hellow\n";
        density=nAtoms/(2*box.x*2*box.y*2*box.z);
    }
//##############################################################################################################################
    neigh_list_update(Atoms,nAtoms);
    HISTOGRAM = new (nothrow) long [nAtoms+1];
    if(Atoms==nullptr) {
        cout<<"Memory Allocation Failed\n";
        return 0;
    }
    for(int i=0; i<nAtoms+1; i++)
    {
        HISTOGRAM[i]=0;
    }
    //cin>>N;
    //cin>>temp;
    //cin>>Press;
    //cin>>nc;
    N=2000000;
    temp=1.0;
    if(restart)
    {
        nc=atoi(argv[4]);
//    index=atoi(argv[5]);
        strncpy(label,argv[5],sizeof(label)-1);
        Press=stof(argv[6]);
    }
    else
    {
        nc=atoi(argv[3]);
//index=atoi(argv[4]);
        strncpy(label,argv[5],sizeof(label)-1);
        Press=stof(argv[5]);
    }
    double EqN;
    if(equilibriate)
        EqN=200000;
    else
        EqN=0;
    vol=2*box.x*2*box.y*2*box.z;
    char buffer[64];
    snprintf(buffer,sizeof(char)*64,"OUT/out_%d_%d_%.2f_%d.dat",int(nAtoms),int(nc),Press,label);//_%d_%f.dat",int(nAtoms),Press);
    //freopen(buffer,"w",stdout);
    cout<<"no of Atoms:"<<nAtoms<<"\n"<<flush;
    cout<<"N:"<<N<<"\n"<<flush;
    cout<<"EqN:"<<EqN<<"\n"<<flush;
    cout<<"Pressure:"<<Press<<"\n"<<flush;
    cout<<"temp:"<<temp<<"\n"<<flush;
    cout<<"box:"<<box.x<<"\n"<<flush;
    cout<<"density:"<<nAtoms/(2*box.x*2*box.y*2*box.z)<<"\n"<<flush;
    cout<<"overlap:"<<check_overlap(Atoms,nAtoms)<<"\n"<<flush;
//############################################################################################
//MPI STUFF
    int my_id,numproc;
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_id);
    MPI_Comm_size (MPI_COMM_WORLD, &numproc);
    nc=10+my_id*10;
    int END=0;
    int rand=0;
    snprintf(buffer,sizeof(char)*64,"OUT/density/density_%d_%d_%.2f_%s.dat",int(nAtoms),int(nc),Press,label);//_%d_%f.dat",int(nAtoms),Press);
    std::ofstream DENSITY(buffer);
    cout<<"nc:"<<nc<<"\n"<<flush;
    snprintf(buffer,sizeof(char)*64,"OUT/clusters/cluster_%d_%d_%.2f_%s.dat",int(nAtoms),int(nc),Press,label);//);//_%d_%f.dat",int(nc),Press);
    std::ofstream CS(buffer);
    int flag=0;
    snprintf(buffer,sizeof(char)*64,"OUT/Histogram/Histogram_%d_%d_%.2f_%s.dat",int(nAtoms),int(nc),Press,label);//);
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
    n=largest_cluster(Atoms,nAtoms,l,box,label); // calculate the largest cluster in the initial config.
    cout<<"n="<<n<<"\n"<<flush;
    int n1;
    Vector dr;
//#####################################################################################################################################3
    back_up(Atoms,old_Atoms,nAtoms);
    if(n<nc)
        tempnc=int(n/10)*10+10;
    if(n>nc)
        tempnc=nc;
    cout<<my_id<<"\t"<<tempnc<<"\n";
    //exit(0);
    for(int i=START; i<(START+EqN); i++) {
        if(fmod(i,20)==0 and !bias) {
            close_reset(Atoms,nAtoms);
            n1=largest_cluster(Atoms,nAtoms,l,box,label);
            HISTOGRAM[n1]++;
            std::ofstream HIS(buffer);
            for(int n=0; n<nAtoms+1; n++)
            {
                HIS<<n<<"\t"<<float(HISTOGRAM[n])<<"\n"<<flush;
            }
            CS<<i<<"\t"<<n1<<"\n"<<flush;
            HIS.close();
        }
        //cout<<i<<"\n";
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
//			cout<<i<<"\t"<<n<<"\t"<<rand<<"\n";
                    mcmove_ver_hardsphere(Atoms,nAtoms);
                    N_iter++;
                }
            }
            else {
                vmove(Atoms,nAtoms);
                density=nAtoms/vol;
                // cout<<box.x<<"\n";
            }
        }
        DENSITY<<i<<"\t"<<M_PI/6.0*density<<"\n"<<flush;
        if(bias)// and (fmod(i,20)==0))
        {
            n=umbrella(Atoms,old_Atoms,nAtoms,l,box,n,flag,HISTOGRAM,label,tempnc);

            back_up(Atoms,old_Atoms,nAtoms); //we need a copy of the config before the move.
            old_box=box;
            // HISTOGRAM[n]+=1;
//	    if(fmod(i,1000)==0)
            CS<<i<<"\t"<<n<<"\n"<<flush;
        }


        if(fmod(i,100)==0)
        {   cout<<i*1.0/(START+EqN)<<"\n"<<flush;
            // if(bias) {
            //     std::ofstream HIS(buffer);
            //     for(int n=0; n<nAtoms; n++)
            //     {
            //         HIS<<n<<"\t"<<float(HISTOGRAM[n])/i<<"\n"<<flush;
            //     }
            //     HIS.close();
            //   }
            print_pos(Atoms,nAtoms,i,label,n);
        }
        if(bias) {
            if(n==tempnc)
	    {
    		cout<<my_id<<"\t"<<tempnc<<"\n";
                tempnc=tempnc+10;
	    //	cout<<"tenmpnc\t"<<tempnc<<"\n";
	    }
            if(n==nc) {
                flag=1;
                break_point=i;
                cout<<"yay\n";
                break;
            }
        }
    }
    if(!bias and equilibriate)
    {
        break_point=START+EqN;
        flag=1;
    }
    flag=1;
    cout<<"eq completed\n";
    cout<<nei_up*1.0/(N_iter)<<"\n";
    Nacc=0;
    Iter=0;
    Nacc_v=0;
    Iter_v=0;
    int count=0;
    long double den_sum=0;
    MPI_Barrier(MPI_COMM_WORLD);
    //#############################################################################################
    //sampling !!!!
    for(int i=break_point; i<(N+break_point); i++) {
        if(fmod(i,20)==0 and !bias) {
            close_reset(Atoms,nAtoms);
            n1=largest_cluster(Atoms,nAtoms,l,box,label);
            HISTOGRAM[n1]++;
            std::ofstream HIS(buffer);
            for(int n=0; n<nAtoms+1; n++)
            {
                HIS<<n<<"\t"<<float(HISTOGRAM[n])<<"\n"<<flush;
            }
            CS<<i<<"\t"<<n1<<"\n"<<flush;
            HIS.close();
        }
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
            Nacc_v=0;
            Iter_v=0;
        }
        if(fmod(i,100)==0)
        {
            int* local;
            int* all;
            local = new (nothrow) int[2*numproc];
            all = new (nothrow) int[2*numproc];
            for(int k=0; k<2*numproc; k++)
            {
                local[k]=0;
                all[k]=0;
            }
            local[my_id+0*numproc]=n;
            local[my_id+1*numproc]=nc;
            //  cout<<my_id<<"\t"<<n<<"\t"<<nc<<"\n";
            MPI_Reduce(&local[0], &all[0], 2*numproc, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
            //  cout<<"first_all\n";
            //  cout<<my_id<<"\t";
            //  for(int k=0; k<2*numproc; k++)
            //      cout<<all[k]<<"\t";
            //  cout<<"\n";
            local[my_id+0*numproc]=0;
            local[my_id+1*numproc]=0;
            if (my_id==0)
            {
                int* CUR_CLU;
                int* BIAS;
                CUR_CLU = new (nothrow) int[numproc];
                BIAS = new (nothrow) int[numproc];
                for(int k=0; k<numproc; k++)
                {
                    CUR_CLU[k]=all[k+0*numproc];
                    BIAS[k]=all[k+1*numproc];
                    //	    cout<<my_id<<"\t"<<CUR_CLU[k]<<"bias\t"<<BIAS[k]<<"\n";
                }
                int direction = 1;
                int start = 0;
                int end = numproc-1;
                if (drand48() > 0.5) {
                    direction = -1;
                    start = 1;
                    end = numproc;
                }
                double lambda=0.15;
                for(int k=start ; k<end ; k++)
                {
                    int nbr=k+direction;
                    double w_o=0.5*lambda*(CUR_CLU[k]-BIAS[k])*(CUR_CLU[k]-BIAS[k])+0.5*lambda*(CUR_CLU[nbr]-BIAS[nbr])*(CUR_CLU[nbr]-BIAS[nbr]);
                    double w_n=0.5*lambda*(CUR_CLU[k]-BIAS[nbr])*(CUR_CLU[k]-BIAS[nbr])+0.5*lambda*(CUR_CLU[nbr]-BIAS[k])*(CUR_CLU[nbr]-BIAS[k]);
                    int swap=0;
                    if(w_n <= w_o)
                        swap=1;
                    else if (exp(-1.*(w_n-w_o)) < drand48()) {
                        swap=0;
                    }
                    else {
                        swap=1;
                    }
                    if(swap)
                    {
                        //cout<<"swap\t"<<k<;
                        double swapn=BIAS[nbr];
                        double swapk=BIAS[k];
                        //double swn=CUR_CLU[nbr];
                        //double swk=CUR_CLU[k];
                        //	cout<<"inswap\t"<<k<<"\t"<<nbr<<"\n";
                        local[k+1*numproc]=swapn;
                        BIAS[k]=swapn;
                        local[k+0*numproc]=1.0;
                        local[nbr+1*numproc]=swapk;
                        BIAS[nbr]=swapk;
                        local[nbr+0*numproc]=1.0;

                    }
                    // cout<<k<<"\t";
                    // for(int k=0; k<2*numproc; k++)
                    //     cout<<local[k]<<"\t";
                    // cout<<"\n";


                }



            }
            // cout<<my_id<<"\t";
            // for(int k=0; k<2*numproc; k++)
            //     cout<<local[k]<<"\t";
            // cout<<"\n";

            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Allreduce (&local[0], &all[0], 2*numproc, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            //  cout<<my_id<<"\t"<<"clu\t"<<all[my_id+0*numproc]<<"\n";
            //  cout<<my_id<<"\t"<<"bias\t"<<all[my_id+1*numproc]<<"\n";
            if(all[my_id+0*numproc])
                nc=all[my_id+1*numproc];
            //  cout<<my_id<<"\t"<<n<<"\t"<<nc<<"\n";
            //cout<<my_id<<"\t";
            //for(int k=0; k<2*numproc; k++)
            //    cout<<all[k]<<"\t";
            //cout<<"\n";

            delete[] local;
            delete[] all;
        }

        for(int n=0; n<nAtoms; n++)     //MC_Sweep
        {
            std::uniform_int_distribution<int> uni(0,nAtoms);
            rand=uni(rng);
            if(rand<nAtoms) {
                {
                    mcmove_ver_hardsphere(Atoms,nAtoms);
                    N_iter++;
                }
            }
            else {
                vmove(Atoms,nAtoms);
                density=nAtoms/vol;
                // cout<<box.x<<"\n";
            }
        }
        DENSITY<<i<<"\t"<<M_PI/6.0*density<<"\n"<<flush;
        if(bias and (fmod(i,20)==0)) {
            n=umbrella(Atoms,old_Atoms,nAtoms,l,box,n,flag,HISTOGRAM,label,nc);
            back_up(Atoms,old_Atoms,nAtoms); //we need a copy of the config before the move.
            old_box=box;
//	    cout<<"yo\t"<<flag<<"\n";
            if(flag)
            {
                HISTOGRAM[n]++;
                count++;
                //            cout<<i<<"\t"<<n<<"\n";
            }
            if(fmod(i,500)==0)
                print_pos(Atoms,nAtoms,i,label,n);
            CS<<i<<"\t"<<n<<"\n"<<flush;
        }
        den_sum+=density;
        if(fmod(i,100)==0)
        {
            if(!bias)
                print_pos(Atoms,nAtoms,i,label,n);
            if(fmod(i,1000)==0)
            {   cout<<(i-break_point)*1.0/N<<"\n"<<flush;
                g_d=pair_correlation(Atoms,nAtoms,1,box,temp,Press);
                std::ofstream HIS(buffer);
                for(int n=0; n<nAtoms+1; n++)
                {
                    HIS<<n<<"\t"<<float(HISTOGRAM[n])<<"\n"<<flush;
                }
                HIS.close();
            }
            else
                g_d=pair_correlation(Atoms,nAtoms,0,box,temp,Press);
//      	cout<<i<<"\t"<<g_d<<"\n";
        }
    }
    cout<<count<<"\n";
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
