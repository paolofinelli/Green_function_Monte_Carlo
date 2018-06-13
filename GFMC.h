/*
 *  GFMC.h
 *  This version uses reconfiguration
 *
 *  Created by Yifei Shi on 4/16/16.
 *  Copyright 2016 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef GFMC_H
#define GFMC_H

#include<fstream>
#include<iostream>
#include<cmath>
#include<stdlib.h>
#include<vector>
#include <sstream>
#include "MersenneTwister.h"
#include "sorting.h"
#include "/Users/yifeishi/software/eigen/Eigen/Dense"

const int Dim = 10;

class GFMC {
private:
	int L_;		//NUmber of reconfig steps
	int M_;		//Number of walkers
	std::vector<int> Walker;
	std::vector<double> Wj;
	std::vector<double> Bj;
	std::vector<double> Ghis;
	int EldestG;
	int NG;
		
	Eigen::MatrixXd H;
	
	MTRand ran;
	
public:
	GFMC(int, int);
	void Diffuse(bool);
	void ReConfig(bool);
	void ApplyO();
	double ExpO();
	double ExpW();
	double ShowWav(int);
	double ShowEav(int);
	double ExpValue(int);
	
	double ShowE(int);
	
	void Init();
	
	

};

GFMC::GFMC(int L, int M)
:	L_(L)
,	M_(M)
,	EldestG(0)
,	NG(0)
{
	Wj.resize(M);
	Bj.resize(M);
	Ghis.resize(L);
	Walker.resize(M);
	
	H = Eigen::MatrixXd::Zero(Dim,Dim);

	for(int i = 0; i < Dim -1; i ++) {
		H(i,i) = 0.4;
		H(i,i+1) = sqrt(i+1)/3.;
		H(i+1,i) = sqrt(i+1)/3.; }
	H(Dim-1,Dim-1) = 0.4;
	
	H(0,Dim-1) = 1;
	H(Dim-1,0) = 1;
	H(2,7) = 0.5;
	H(7,2) = 0.5;
	H(6,3) = 0.5;
	H(3,6) = 0.5;
	
	for(int k = 0; k < M; k ++) {
		Walker[k] = int(ran.randDblExc()*Dim);
		Wj[k] = 1.;
		Bj[k] = 0.;
	}
	
	Eigen::MatrixXd O = Eigen::MatrixXd::Zero(Dim,Dim);
	
	O(Dim-1,Dim-1) = 1;
	
	std::cout << O << "\n";

	
	std::cout << H << "\n";
	
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esm(H);
	
	std::cout << "Exact value: " << esm.eigenvalues() << "\n";
	
	std::cout << "Expectation of O:  " << (esm.eigenvectors().col(Dim-1)).transpose() * O * esm.eigenvectors().col(Dim-1) << "\n";
}

void GFMC::Diffuse(bool ms) {
	for(int mi = 0; mi < M_; mi ++) {
		
		double Bx = 0., Bxp = 0.;
		
		for(int j = 0; j < Dim; j ++) {
			Bx += H(j,Walker[mi]); 
		}
		
		//Wj[mi] *= Bx;
		
		if(Bx == 0.) std::cout << "ZZ!!!\n";
		
		double p = ran.randDblExc();
		
		double where = 0.;
		for(int j = 0; j < Dim; j ++) {
			where += H(j,Walker[mi]);
			if(where/Bx > p) {
				Walker[mi] = j;
				break;
			}
		}
		
		for(int i = 0; i < Dim; i ++) {
			Bxp += H(i,Walker[mi]); }
		
		Wj[mi] *= Bx;
		Bj[mi] = Bxp;
				
	} //end for mi
}

void GFMC::ReConfig(bool method) {
	if(method) {
		double W = 0.;
		for(int mi = 0; mi < M_; mi ++) {
			W += Wj[mi]; }
		
		int* nconfig = (int*)malloc(sizeof(int)*M_);
		
		for(int mi = 0; mi < M_; mi ++) {
			double r = ran.randDblExc();
			double p = 0.;
			for(int ni = 0; ni < M_; ni ++) {
				p += Wj[ni]/W;
				if(r < p) {
					nconfig[mi] = Walker[ni];
					break;
				}
			}
		}
		
		for(int ni = 0; ni < M_; ni ++) {
			Walker[ni] = nconfig[ni];
		}
		
		for(int mi = 0; mi < M_; mi ++)
			//Wj[mi] = W/double(M_);	
			Wj[mi] = 1.;
		
		if(NG < L_) {
			Ghis[NG] = W;
			NG ++;
		}
		else {
			Ghis[EldestG] = W;
			EldestG = EldestG < L_-1 ? EldestG+1 : 0;
		}
	}
	
	else {
		std::vector<double> randarray(M_);
		
		for(int i = 0; i < M_; i ++) 
			randarray[i] = ran.randDblExc();
		
		std::vector<int> idx = MergeSort(randarray, M_);
		
		double W = 0.;
		for(int mi = 0; mi < M_; mi ++) 
			W += Wj[mi];
		
		int* nconfig = (int*)malloc(sizeof(int)*M_);
		
		double cp = Wj[0]/W;
		int j = 0;
		int currentidx = 0;
	
		while(currentidx < M_) {
			while(cp < randarray[idx[currentidx]]) {
				j ++;
				cp += Wj[j]/W; }
			nconfig[currentidx] = Walker[j];
			currentidx ++;
		}
		
		for(int ni = 0; ni < M_; ni ++) {
			Walker[ni] = nconfig[ni];
		}
		
		for(int mi = 0; mi < M_; mi ++)
			//Wj[mi] = W/double(M_);	
			Wj[mi] = 1.;
		
		if(NG < L_) {
			Ghis[NG] = W;
			NG ++;
		}
		else {
			Ghis[EldestG] = W;
			EldestG = EldestG < L_-1 ? EldestG+1 : 0;
		}
	}
}

void GFMC::ApplyO() {
	for(int mi = 0; mi < M_; mi ++) {
		Walker[mi] = Walker[mi] < Dim/2 ? Walker[mi] + Dim/2 : Walker[mi] - Dim/2;
	}
}

//a dumb operator, for test
double GFMC::ExpValue(int k) {
			
	if(NG < L_) {
		std::cout << "Not Thermalized!\n";
		return Walker[k]==9? Wj[k] : 0.;
	}
	else {
		double G = 1.;
		for(int i = 0; i < L_; i ++)
			G *= Ghis[i];
		
		return  Walker[k] == 9 ? Wj[k]*G : 0.;
	}
	
}

double GFMC::ExpW() {
	double wt = 0.;
	double G = 1.;
	for(int i = 0; i < L_; i ++)
		G *= Ghis[i];
	
	for(int mi = 0; mi < M_; mi ++) {
		wt += Wj[mi]; }
	return wt*G;
}

double GFMC::ExpO() {
/*	double owt = 0.;
	for(int mi = 1; mi < M_; mi ++) {
		for(int ni = 0; ni < mi; ni ++) {
			if(abs(Walker[mi]-Walker[ni]) == Dim/2)
				owt += sqrt(Wj[mi]*Wj[ni]);
		}
	}
	return owt; */
	double G = 1.;
	for(int i = 0; i < L_; i ++)
		G *= Ghis[i];
	
	double wt = 0.;
	for(int mi = 0; mi < M_; mi ++) {
		if(Walker[mi] == 9)
			wt += Wj[mi]; }
	return wt*G;
	
}

double GFMC::ShowEav(int k) {
	if(NG < L_) {
		std::cout << "Not Thermalized!\n";
		return Wj[k]*Bj[k];
	}
	else {
		double G = 1.;
		for(int i = 0; i < L_; i ++)
			G *= Ghis[i];
		
		return Wj[k]*Bj[k]*G;
	}
}

double GFMC::ShowWav(int k) {
	if(NG < L_) {
		std::cout << "Not Thermalized!\n";
		return Wj[k];
	}
	else {
		double G = 1.;
		for(int i = 0; i < L_; i ++)
			G *= Ghis[i];
		
		return Wj[k]*G;
	}
}

void GFMC::Init(){
	for(int mi = 0; mi < M_; mi ++) {
		Walker[mi] = int(ran.randDblExc()*Dim);
		Wj[mi] = 1.;
		Bj[mi] = 0.;
		NG = 0;
		EldestG = 0;
	}
}
	

#endif
