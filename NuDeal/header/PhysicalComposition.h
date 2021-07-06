#pragma once
#include "Defines.h"

namespace PhysicalComposition {
struct BaseMat_t{
	int niso;
	double temperature;
	vector<int> idiso;
	vector<double> pnum;
};
struct ThermoHydroMat_t {
	int nmix;
	double temperature;
	vector<int> idmix;
	vector<double> pmass;
	// For macro type mateiral
	double k_c, h_t;
};
struct Boundary_t {
	// 0 : vacuum, 1 : reflective w fraction, 2 : white w fraction, 3 : zero-flux
	// 4 : checker board, 5 : rotational
	int BaseType; 
	double lossrate;
	// 0 : Dirichlet (T0), 1 : Neumann (Laplacian(Q)=q0)
	int THType;
	double boundval;
};
class MatComp {
private:
	vector<BaseMat_t> BaseMat;
	vector<ThermoHydroMat_t> THMat;
	vector<Boundary_t> Boundary;
public:
	void PushBaseMacro(bool isTHFeed, double temperature) 
	{ int i = BaseMat.size(); BaseMat.emplace_back(); if (isTHFeed) BaseMat[i].temperature = temperature; }

	void PushBoundary(int BaseType, double lossrate);

	const BaseMat_t& GetBaseMat(int i) { return BaseMat[i]; }

	const Boundary_t& GetBoundary(int i) { return Boundary[-1 - i]; }
};
}