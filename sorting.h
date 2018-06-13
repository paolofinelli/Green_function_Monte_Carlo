/*
 *  sorting.h
 *  
 *
 *  Created by Yifei Shi on 5/12/16.
 *  Copyright 2016 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef SORTING_H
#define SORTING_H


#include<cmath>
#include<vector>
#include<stdlib.h>

//return the new location array; changes Npiece to new value
void ArrayMerge(const std::vector<double> &usortarray, std::vector<int> &idx, std::vector<int> &location, int &Npiece, const int &N) {
	int cp1 = 0, cp2 = 1;
	
	std::vector<int> Newlocation;
	
	Newlocation.push_back(0);
	
	int l1 , l2;			//location of the two small arrays
	
	double N1, N2;		//the two numbers to be compared
	
	while(cp1 < Npiece) {
		
		l1 = location[cp1];
		l2 = location[cp2];
		
		int endloc = cp2 < Npiece ? location[cp2+1]-1 : N-1;
		
		int* ipart = (int*)malloc(sizeof(int)*(endloc-location[cp1]));
		
		int Npsize = 0;		//size of the combined array
		
		int reachend = 0;
		
		while(l1 < location[cp2] || l2 < endloc+1) {
			
			if(reachend == 0) {
				N1 = usortarray[idx[l1]];
				N2 = usortarray[idx[l2]];
				
				
				if(N1 < N2) {
					ipart[Npsize] = idx[l1];
					l1 ++;
					if(l1 >= location[cp2])
						reachend = 1;
				}
				else {
					ipart[Npsize] = idx[l2];
					l2 ++; 
					if(l2 >= endloc+1)
						reachend = 2;
				}
				
				Npsize ++;
			}
			
			else if(reachend == 1) {
				ipart[Npsize] = idx[l2];
				l2 ++; 
				Npsize ++;
			}
			else {
				ipart[Npsize] = idx[l1];
				l1 ++;
				Npsize ++;
			}
			
		}		//combine 2 arrays
		
		for(int i = 0; i < Npsize; i ++) {
			idx[location[cp1]+i] = ipart[i];
		}
		
		if(cp2 < Npiece) 
			Newlocation.push_back(l2);
		cp1 += 2;
		cp2 += 2;
	}		//a whole sweep
	
	Npiece = Npiece%2 == 1 ? (Npiece+1)/2-1 : (Npiece+1)/2;
	location.resize(Npiece);
	location = Newlocation;
	
	
}

//input an unsorted array, and the number of elements. output an int array, j(i) gives the index of the i th smallest element

std::vector<int> MergeSort(const std::vector<double> &usortarray, const int &N) {
	int Npiece = N-1;		//number of small arrays
	
	std::vector<int> result;
	
	std::vector<int> location;		//location of the beginning of small arrays
	
	result.resize(N);
	
	location.resize(N);
	
	for(int i = 0; i < N; i ++) {
		result[i] = i;
		location[i] = i;
	}
	
	while(Npiece > 0) {
		ArrayMerge(usortarray, result, location, Npiece, N);
	}
	
	return result;
}



#endif