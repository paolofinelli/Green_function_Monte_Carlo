/*
 *  sorttest.cpp
 *  
 *
 *  Created by Yifei Shi on 5/12/16.
 *  Copyright 2016 __MyCompanyName__. All rights reserved.
 *
 */


#include<iostream>
#include<stdio.h>
#include "sorting.h"
#include "MersenneTwister.h"

int main(int argc, char* argv[]) {
	int N;
	MTRand ran;
	
	sscanf(argv[1], "%d", &N);
	
	std::vector<double> a(N);
	for(int i = 0; i < N; i ++) {
		a[i] = ran.randDblExc();}
	
	std::vector<int> idx = MergeSort(a, N);
	
	
	for(int i = 0; i < N; i ++)
		std::cout << a[i] << ", ";
	std::cout << "\n";

	for(int i = 0; i < N; i ++)
		std::cout << a[idx[i]] << ", ";
		
	std::cout << "\n";
		
		
		
}