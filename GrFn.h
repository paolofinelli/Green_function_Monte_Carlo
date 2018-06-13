/*
 *  GrFn.h
 *  
 *
 *  Created by Yifei Shi on 4/16/16.
 *  Copyright 2016 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef GRFN_H
#define GRFN_H

#include<fstream>
#include<iostream>
#include<cmath>
#include<stdlib.h>
#include<vector>
#include <sstream>
#include "MersenneTwister.h"
#include "/Users/yifeishi/software/eigen/Eigen/Dense"

const int Dim = 9;

class GrFn {
private:
	int L_;		//size of Bhis
	double Wav;
	double Eav;		//average of w and E;
	double* Bhis;	//Bx history
	
	int eldest;	//position eldest element in Bx
	
	int Walker;
	
	int BxCt; //size of Bhist
	
	Eigen::MatrixXd H;
	
	MTRand ran;
	
public:
	GrFn(int);
	void Diffuse(bool);
	double ShowWav();
	double ShowEav();
	
	double ShowE();
	
	void Init();
	
	

};

GrFn::GrFn(int L):
L_(L),
Wav(0.),
Eav(0.),
eldest(0),
BxCt(0)
{
	H = Eigen::MatrixXd::Zero(Dim,Dim);
	
	Bhis = (double*) malloc(sizeof(double)*(L+1));
	
	for(int i = 0; i < 9; i ++) {
		H(i, i + 1) = sqrt(i+1)/3.;
		H(i + 1, i) = sqrt(i+1)/3.;
	}
	
	H(9,0) = 1.;
	H(0,9) = 1.;
	
	H(8,2) = 0.5;
	H(2,8) = 0.5;
		
	Walker = int(ran.randDblExc()*Dim);
	
	//Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esm(H);
	
	//std::cout << "Exact value: " << esm.eigenvalues()[9] << "\n";
}

void GrFn::Diffuse(bool ms) {
	double Bx = 0., Bxp = 0.;
	for(int i = 0; i < Dim; i ++) {
		Bx += H(i,Walker); 
		
	}
	if(Bx ==0) std::cout << "ZZ!!!\n";

	if(BxCt < L_+1) {
		BxCt ++;
		Bhis[BxCt-1] = Bx;
		 }
	else {
		Bhis[eldest] = Bx;
		eldest = eldest < L_ ? eldest +1 : 0;
	}
	
	double p = ran.randDblExc();
	
	double where = 0.;
	for(int j = 0; j < Dim; j ++) {
		where += H(j,Walker);
		if(where/Bx > p) {
			Walker = j;
			break;
		}
	}
	
	for(int i = 0; i < Dim; i ++) {
		Bxp += H(i,Walker); }
	
	if(ms) {
		double py = 1;
		for(int j = 0; j < L_+1; j ++) {
			py *= Bhis[j]; }
		
		Eav += py;
		Wav += py/Bx;
		if(Eav ==0 || Wav ==0)
			std::cout << "W!!\n";
	}
}

double GrFn::ShowE() {
	return Eav/Wav;
}

double GrFn::ShowEav() {
	return Eav;
}

double GrFn::ShowWav() {
	return Wav;
}

void GrFn::Init(){
	BxCt = 0;
	Walker = int(ran.randDblExc()*Dim);
	eldest = 0;
	Wav = 0.;
	Eav = 0.;
}
	

#endif
