#ifndef MC_HEADER_INCLUDED_H
#define MC_HEADER_INCLUDED_H
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
const double sigma2=sigma*sigma;
const double epsilon=1;
const double m=1;
const double k_b=1;
const double unit_time=pow(m*sigma*sigma/epsilon,0.5);
double temp;
double JUM=0.05;
Vector box;
Vector old_box;
double density=0;
long nAtoms;
std::ofstream MSD("position_x.dat");
std::random_device rd;
std::mt19937 rng(rd());
std::uniform_real_distribution<double> uni_d(0.0,1.0);
///////////////////////////////////////////////////////
double DELTA_T=.001,r_cut=2.5*sigma;//pow(2.0,1.0/6);
double RCUT_sq=r_cut*r_cut;
///////////////////////////////////////////////////////
double R_CUT_HS=sigma;
double R_CUT_HS_sq=R_CUT_HS*R_CUT_HS;

double R_V=1.5;
double R_V_sq=1.5*1.5;

double gap=(R_V-R_CUT_HS)/2.0;
double gap_sq=(gap)*(gap);

double PE;
double vol;
double dlnV=0.002;
double vir_pre=0;
double percor;
double Press=3.0;
long Nacc=0,Iter=0;
long Nacc_v=0,Iter_v=0;
long nc=80;
int l=6;
int nei_up=0;
#endif
