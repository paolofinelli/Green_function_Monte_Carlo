/*
 *  GFMCO.cpp
 *  This version tries to calculate expectation value of 
 *
 *  Created by Yifei Shi on 4/16/16.
 *  Copyright 2016 __MyCompanyName__. All rights reserved.
 *
 */

#include "GFMC.h"

int main(int argc, char* argv[]) {
		
	int kp, Nm, L, M;
	bool apl;
	int Nk;
	
	sscanf(argv[1], "%d", &kp);
	sscanf(argv[2], "%d", &Nm);
	sscanf(argv[3], "%d", &L);
	sscanf(argv[4], "%d", &M);
	sscanf(argv[5], "%d", &Nk);
	sscanf(argv[6], "%d", &apl);
	
	double Eav = 0., Wav = 0.;
	
	GFMC a(L, M);

	for(int k = 0; k < Nk; k ++) {
	
		for(int l = 0; l < Nm+20; l ++) {
			for(int j = 0; j < kp; j ++) {
				a.Diffuse(0);				
			}
			a.ReConfig(0);			
		}
		
		if(apl)
			a.ApplyO();
		
		for(int i = 0; i < Nm; i ++) {
			for(int j = 0; j < kp; j ++) {
				a.Diffuse(0);				
			}
			a.ReConfig(0);
		}
		
		
		for(int mi = 0; mi < M; mi ++) {
		Eav = Eav*double(k)/double(k+1) + a.ShowEav(mi)/double(M)/double(k+1);
		
		Wav = Wav*double(k)/double(k+1) + a.ShowWav(mi)/double(M)/double(k+1);
		}
		
		
		a.Init();
		
	}
			

	
	
		
	std::cout << "Finally: Eav: " << Eav << "\t Wav: " << Wav << "\t E: "<< Eav/Wav << "\n";
	
	/*Eigen::MatrixXd H = Eigen::MatrixXd::Zero(Dim,Dim);
		
	H = Eigen::MatrixXd::Zero(Dim,Dim);
	
	for(int i = 0; i < Dim -1; i ++) {
		H(i,i+1) = sqrt(i+1)/3.;
		H(i+1,i) = sqrt(i+1)/3.; }
	
	H(0,Dim-1) = 1;
	H(Dim-1,0) = 1;
	H(1,7) = 0.5;
	H(7,1) = 0.5;
	H(3,6) = 0.5;
	H(6,3) = 0.5;
	
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esm(H);
	
	std::cout << "Exact value: " << esm.eigenvalues() << "\n";
	
	std::cout << H << "\n"; */
	
	return 1;
	
	
}
