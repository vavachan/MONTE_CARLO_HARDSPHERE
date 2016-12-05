#include<iostream>
#include<math.h>
#include<stdlib.h>
#include<fstream>
#include<time.h>
#include<ctime>
#include<random>
#include "vector.hpp"
#include "atom.hpp"
const double sigma=1;
const double epsilon=1;
const double m=1;
const double k_b=1;
const double unit_time=pow(m*sigma*sigma/epsilon,0.5);
double temp;
double JUM=.15;
Vector box;
double density=0;
long nAtoms;
std::ofstream MSD("position_x.dat");
std::ofstream POSITION("config.dat");
std::random_device rd;
std::mt19937 rng(rd());
std::uniform_real_distribution<double> uni_d(0.0,1.0);
double DELTA_T=.001,r_cut=2.5;//pow(2.0,1.0/6);
double R_CUT=1;
double R_CUT_HS=1;
double RCUT_sq=r_cut*r_cut;
long Nacc=0,Iter=0;
