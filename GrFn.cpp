/*
 *  GrFn.cpp
 *  
 *
 *  Created by Yifei Shi on 4/16/16.
 *  Copyright 2016 __MyCompanyName__. All rights reserved.
 *
 */

#include "GrFn.h"

int main(int argc, char* argv[]) {
		
	int Nt, Nm, L;
	int Nk;
	
	sscanf(argv[1], "%d", &Nt);
	sscanf(argv[2], "%d", &Nm);
	sscanf(argv[3], "%d", &L);
	sscanf(argv[4], "%d", &Nk);
	
	double E = 0., Eav = 0., Wav = 0.;
	
	GrFn a(L);

	for(int k = 0; k < Nk; k ++) {
		for(int j = 0; j < Nt; j ++) {
			a.Diffuse(0);
			
		}
		
		for(int j = 0; j < Nm; j ++) {
			a.Diffuse(1);
			
		} 
		Eav = Eav*double(k)/double(k+1) + a.ShowEav()/double(k+1);
		Wav = Wav*double(k)/double(k+1) + a.ShowWav()/double(k+1);
		
		a.Init();
	}
		
	std::cout << "Finally: Eav" << Eav << "\tWav: " << Wav << "\t E: "<< Eav/Wav << "\n";
	
	Eigen::MatrixXd H = Eigen::MatrixXd::Zero(10,10);
		
	for(int i = 0; i < 9; i ++) {
		H(i, i + 1) = sqrt(i+1)/3.;
		H(i + 1, i) = sqrt(i+1)/3.;
	}
	H(9,0) = 1.;
	H(0,9) = 1.;
	
	H(8,2) = 0.5;
	H(2,8) = 0.5;
		
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esm(H);
	
	std::cout << "Exact value: " << esm.eigenvalues() << "\n";
	
	std::cout << H << "\n";
	
	return 1;
	
	
}
