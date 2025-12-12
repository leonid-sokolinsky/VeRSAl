/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: VeRSAl (Vertex Retrieve by Simplex Algorithm) No MPI
Module: Problem-bsfTypes.h (Predefined BSF Problem Types)
Prefix: PT_bsf
Author: Leonid B. Sokolinsky
This source code is a part of BSF Skeleton
==============================================================================*/
#pragma once
#include "Problem-Types.h"		// Problem Types 
//=========================== BSF Types =========================
struct PT_bsf_parameter_T {		// Type of Parameter for workers (current approximation)
	PT_vector_T v;				// Current vertex
	PT_vector_i_T basis_v;		// Current basis
};

struct PT_bsf_mapElem_T {		// Type of map-list elements
	int entryNumber;
};

struct PT_bsf_reduceElem_T {	// Type of reduce-list elements for Job 0 (default)	
	PT_vector_T  v_nex;
	int i_star;
	int j_star;
	double objF_nex;	// F(v_nex)
	double objF_grd;	// Value of objective function after one unit movement
	bool optimumIsFound;
};

struct PT_bsf_reduceElem_T_1 {	// Type of reduce-list elements for Job 1
	// Not used
};

struct PT_bsf_reduceElem_T_2 {	// Type of reduce-list elements for Job 2
	// Not used
};

struct PT_bsf_reduceElem_T_3 {	// Type of reduce-list elements for Job 3
	// Not used
};