#pragma once
#include "Defines.h"
#include "Array.h"
#include "PhysicalDomain.h"
#include "ShortCharacteristics.h"

namespace SolutionDriver {

inline double normvec(double *vec, int n) {
	double norm = 0.;
	for (int i = 0; i < n; i++) {
		norm += Math::square(vec[i]);
	}
	return sqrt(norm);
}

class SteadyStateDriver{
public:
	template <typename T> using Array = LinPack::Array_t<T>;
	using RayTracingDomain = PhysicalDomain::RayTracingDomain;
	using FluxBoundary = PhysicalDomain::FluxBoundary;
	using FlatSrcDomain = PhysicalDomain::FlatSrcDomain;
	using FlatXSDomain = PhysicalDomain::FlatXSDomain;
	using DriverSCMOC = Transport::DriverSCMOC;
	using Dimension = Geometry::Dimension;

private:
	RayTracingDomain *Rays;
	FluxBoundary *BndFlux;
	FlatSrcDomain *FSR;
	FlatXSDomain *FXR;

	DriverSCMOC *MOCDriver;

	Dimension mode;

	double keff;
	int ng, nangle_oct, scatord;

	int nmaxiter;
	vector<int> niterOut, niterIn;
	bool givenpsi, givenres;
	double crit_k, crit_psi, crit_res;
	Array<double> prevpsi;

	double RelativeKErr(double keff1, double keff0) { return abs(keff1 - keff0); }

	double RelativePsiErr();

	double Residual();

public:
	SteadyStateDriver(Dimension mode, int ng, int nangle_oct, int scatord)
	{	this->mode = mode; this->ng = ng; this->nangle_oct = nangle_oct; this->scatord = scatord;	}

	void Initialize(RayTracingDomain &Rays, FluxBoundary &BndFlux, FlatSrcDomain &FSR, FlatXSDomain &FXR, DriverSCMOC &MOCDriver);

	void SetCriterionK(double crit_k) { this->crit_k = crit_k; }

	void SetCriterionPsi(double crit_psi)
	{ this->crit_psi = crit_psi; givenpsi = true; prevpsi.Create((*this->FXR).GetNblocks()); }

	void setCriterionRes(double crit_rse) { this->crit_res = crit_res; givenres = true; }

	void SetMaxIter(int nmaxiter) { this->nmaxiter = nmaxiter; }

	void RunSS();
};

}