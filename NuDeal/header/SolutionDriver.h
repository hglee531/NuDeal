#pragma once
#include "Defines.h"
#include "SolutionDefines.h"
#include "Array.h"
#include "PhysicalDomain.h"
#include "ShortCharacteristics.h"
#include "XSLibrary.h"

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
	using DriverSN = Transport::DriverSN;
	using Dimension = Geometry::Dimension;

private:
	RayTracingDomain *Rays;
	FluxBoundary *BndFlux;
	FlatSrcDomain *FSR;
	FlatXSDomain *FXR;

	DriverSN *SnDriver;

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

	void Initialize(RayTracingDomain &Rays, FluxBoundary &BndFlux, FlatSrcDomain &FSR, FlatXSDomain &FXR, DriverSN &SnDriver);

	void SetCriterionK(double crit_k) { this->crit_k = crit_k; }

	void SetCriterionPsi(double crit_psi)
	{ this->crit_psi = crit_psi; givenpsi = true; prevpsi.Create((*this->FXR).GetNblocks()); }

	void setCriterionRes(double crit_rse) { this->crit_res = crit_res; givenres = true; }

	void SetMaxIter(int nmaxiter) { this->nmaxiter = nmaxiter; }

	void RunSS();
};

inline void DebugSolutionPin() {
	using namespace Geometry;
	using namespace PhysicalDomain;
	using namespace Transport;

	double c_cir[3] = { 0.63, 0.63, 0.54 };
	UnitSurf Circle(UnitSurf::SurfType::CIRCLE, c_cir, CartPlane::XY);
	
	double c_cir1[3] = { 0.63, 0.63, 0.6 };
	UnitSurf Circle1(UnitSurf::SurfType::CIRCLE, c_cir1, CartPlane::XY);
	double c_cir2[3] = { 0.63, 0.63, 0.63 };
	UnitSurf Circle2(UnitSurf::SurfType::CIRCLE, c_cir2, CartPlane::XY);
	double c_cir3[3] = { 0.63, 0.63, 0.7 };
	UnitSurf Circle3(UnitSurf::SurfType::CIRCLE, c_cir3, CartPlane::XY);

	double c_zpln0[2] = { -1.0, 0.0 }, c_zpln1[2] = { 1.0, 1.26 };
	UnitSurf Zpln0(UnitSurf::SurfType::ZPLN, c_zpln0), Zpln1(UnitSurf::SurfType::ZPLN, c_zpln1);

	double c_xpln0[2] = { -1.0, 0.0 }, c_xpln1[2] = { 1.0, 1.26 };
	double c_ypln0[2] = { -1.0, 0.0 }, c_ypln1[2] = { 1.0, 1.26 };
	UnitSurf Xpln0(UnitSurf::SurfType::XPLN, c_xpln0), Xpln1(UnitSurf::SurfType::XPLN, c_xpln1);
	UnitSurf Ypln0(UnitSurf::SurfType::YPLN, c_ypln0), Ypln1(UnitSurf::SurfType::YPLN, c_ypln1);

	UnitVol Cylinder;
	Cylinder << Circle << Zpln0 << Zpln1;
	Cylinder.Finalize();

	//UnitVol Cylinder1;
	//Cylinder1 << Circle1 << Zpln0 << Zpln1;
	//Cylinder1.Finalize();

	//UnitVol Cylinder2;
	//Cylinder2 << Circle2 << Zpln0 << Zpln1;
	//Cylinder2.Finalize();

	//UnitVol Cylinder3;
	//Cylinder3 << Circle3 << Zpln0 << Zpln1;
	//Cylinder3.Finalize();

	UnitVol Box;
	Box << Zpln0 << Zpln1 << Ypln0 << Ypln1 << Xpln0 << Xpln1;
	Box.Finalize();

	bool isbounded;
	double2 xs, ys, zs;
	isbounded = Cylinder.GetBoundBox(xs, ys, zs);
	cout << "Bound info, (" << isbounded << ")" << endl;
	if (isbounded) cout << "X : [" << xs.x << "," << xs.y << "] ";
	if (isbounded) cout << "Y : [" << ys.x << "," << ys.y << "] ";
	if (isbounded) cout << "Z : [" << zs.x << "," << zs.y << "] " << endl;
	if (isbounded) cout << "Vol = " << Cylinder.GetVolume() << "cm^3" << endl;

	isbounded = Box.GetBoundBox(xs, ys, zs);
	cout << "Bound info, (" << isbounded << ")" << endl;
	if (isbounded) cout << "X : [" << xs.x << "," << xs.y << "] ";
	if (isbounded) cout << "Y : [" << ys.x << "," << ys.y << "] ";
	if (isbounded) cout << "Z : [" << zs.x << "," << zs.y << "] " << endl;
	if (isbounded) cout << "Vol = " << Box.GetVolume() << "cm^3" << endl;

	//Cylinder Volume = 1.15427140641 cm3
	//UnitComp BoxComp(0, Box), CylComp(0, Cylinder);
	UnitComp BoxComp(6, Box), CylComp(0, Cylinder);
	//CylComp.Finalize();
	BoxComp.Finalize();
	
	//UnitComp CylComp1(6, Cylinder1), CylComp2(6, Cylinder2), CylComp3(6, Cylinder3);
	//CylComp1.Finalize();
	//CylComp2.Finalize();
	//CylComp3.Finalize();

	GeometryHandler GeoHandle;
	double origin[3] = { 0.0, 0.0, 0.0 }, L[3] = { 1.26, 1.26, 1.26 };
	GeoHandle.SetOrdinates(origin, L);
	GeoHandle << BoxComp << CylComp;
	//GeoHandle << CylComp1 << CylComp2 << CylComp3;
	//GeoHandle << CylComp2;
	GeoHandle.FinalizeVolumes();
	//GeoHandle.Discretize(Dimension::TwoD, 1.27, 0.1, 1);
	//int3 blockFSR = make_int3(1, 1, 1), blockFXR = make_int3(1, 1, 1);
	//GeoHandle.Discretize(Dimension::TwoD, 0.5, 0.4, 1);
	//int3 blockFSR = make_int3(3, 3, 1), blockFXR = make_int3(3, 3, 1);
	//GeoHandle.Discretize(Dimension::TwoD, 0.01, 0.01, 5);
	//int3 blockFSR = make_int3(5, 5, 1), blockFXR = make_int3(5, 5, 1);
	GeoHandle.Discretize(Dimension::TwoD, 0.1, 0.001, 1);
	int3 blockFSR = make_int3(1, 1, 1), blockFXR = make_int3(1, 1, 1);

	std::cout << "nnode = " << GeoHandle.GetNnode() << ", Lv = " << GeoHandle.GetDivLevel() << std::endl;

	RayTracingDomain RaySP(GeoHandle);
	FluxBoundary BndFlux(Dimension::TwoD, RaySP.GetNnodeZero().x, RaySP.GetNnodeZero().y, RaySP.GetNnodeZero().z);

	FlatSrcDomain FSR(blockFSR, RaySP);
	FlatXSDomain FXR(blockFXR, FSR);

	Library::XSLibrary_t XSLib;
	XSLib.ReadMacro();

	vector<int> quadparam(2);
	quadparam[0] = 16; quadparam[1] = 4;
	//int nangle_oct = CalNangleOct(AngQuadType::UNIFORM, quadparam);
	int nangle_oct = CalNangleOct(AngQuadType::Bi3, quadparam);

	BndFlux.Initialize(XSLib.GetNumGroup(), nangle_oct);
	FSR.Initialize(XSLib.GetNumGroup(), nangle_oct);
	FXR.InitializeMacro(XSLib.GetNumGroup(), XSLib.GetScatOrder(), false);
	FXR.SetMacroXS(XSLib, GeoHandle.GetMatIds());

	//DriverSN SnDriver(XSLib.GetNumGroup(), AngQuadType::UNIFORM, quadparam);
	DriverSN SnDriver(XSLib.GetNumGroup(), AngQuadType::Bi3, quadparam);
	SnDriver.Initialize(RaySP, BndFlux, FSR, FXR);
	//SnDriver.RaySweep();

	SteadyStateDriver SSDriver(Dimension::TwoD, XSLib.GetNumGroup(), nangle_oct, XSLib.GetScatOrder());
	SSDriver.Initialize(RaySP, BndFlux, FSR, FXR, SnDriver);
	SSDriver.SetCriterionPsi(1.e-5);
	SSDriver.SetCriterionK(5.e-8);
	SSDriver.SetMaxIter(10000);

	SSDriver.RunSS();
}

inline void DebugSolutionBox() {
	using namespace Geometry;
	using namespace PhysicalDomain;
	using namespace Transport;

	double c_zpln0[2] = { -1.0, 0.0 }, c_zpln1[2] = { 1.0, 1.26 };
	UnitSurf Zpln0(UnitSurf::SurfType::ZPLN, c_zpln0), Zpln1(UnitSurf::SurfType::ZPLN, c_zpln1);

	double c_xpln0[2] = { -1.0, 0.0 }, c_xpln1[2] = { 1.0, 1.26 };
	double c_ypln0[2] = { -1.0, 0.0 }, c_ypln1[2] = { 1.0, 1.26 };
	UnitSurf Xpln0(UnitSurf::SurfType::XPLN, c_xpln0), Xpln1(UnitSurf::SurfType::XPLN, c_xpln1);
	UnitSurf Ypln0(UnitSurf::SurfType::YPLN, c_ypln0), Ypln1(UnitSurf::SurfType::YPLN, c_ypln1);

	double c_xpln2[2] = { -1.0, 0.0 }, c_xpln3[2] = { 1.0, 2.52 };
	double c_ypln2[2] = { -1.0, 0.0 }, c_ypln3[2] = { 1.0, 2.52 };
	UnitSurf Xpln2(UnitSurf::SurfType::XPLN, c_xpln2), Xpln3(UnitSurf::SurfType::XPLN, c_xpln3);
	UnitSurf Ypln2(UnitSurf::SurfType::YPLN, c_ypln2), Ypln3(UnitSurf::SurfType::YPLN, c_ypln3);

	double c_xpln4[2] = { -1.0, 0.0 }, c_xpln5[2] = { 1.0, 1.2 };
	double c_ypln4[2] = { -1.0, 0.0 }, c_ypln5[2] = { 1.0, 1.2 };
	UnitSurf Xpln4(UnitSurf::SurfType::XPLN, c_xpln4), Xpln5(UnitSurf::SurfType::XPLN, c_xpln5);
	UnitSurf Ypln4(UnitSurf::SurfType::YPLN, c_ypln4), Ypln5(UnitSurf::SurfType::YPLN, c_ypln5);

	UnitVol Box;
	Box << Zpln0 << Zpln1 << Ypln0 << Ypln1 << Xpln0 << Xpln1;
	Box.Finalize();

	UnitVol LBox;
	LBox << Zpln0 << Zpln1 << Ypln2 << Ypln3 << Xpln2 << Xpln3;
	LBox.Finalize();
	
	UnitVol SBox;
	SBox << Zpln0 << Zpln1 << Ypln4 << Ypln5 << Xpln4 << Xpln5;
	SBox.Finalize();
	
	UnitComp BoxComp(0, Box), LBoxComp(6, LBox), SBoxComp(0, SBox);
	BoxComp.Finalize();
	LBoxComp.Finalize();
	SBoxComp.Finalize();

	GeometryHandler GeoHandle;
	double origin[3] = { 0.0, 0.0, 0.0 }, L[3] = { 2.52, 2.52, 1.26 };
	GeoHandle.SetOrdinates(origin, L);
	GeoHandle << BoxComp << LBoxComp << SBoxComp;
	GeoHandle.FinalizeVolumes();
	GeoHandle.Discretize(Dimension::TwoD, 0.1, 0.05, 1);
	int3 blockFSR = make_int3(1, 1, 1), blockFXR = make_int3(1, 1, 1);

	std::cout << "nnode = " << GeoHandle.GetNnode() << ", Lv = " << GeoHandle.GetDivLevel() << std::endl;

	RayTracingDomain RaySP(GeoHandle);
	FluxBoundary BndFlux(Dimension::TwoD, RaySP.GetNnodeZero().x, RaySP.GetNnodeZero().y, RaySP.GetNnodeZero().z);

	FlatSrcDomain FSR(blockFSR, RaySP);
	FlatXSDomain FXR(blockFXR, FSR);

	Library::XSLibrary_t XSLib;
	XSLib.ReadMacro();

	vector<int> quadparam(2);
	quadparam[0] = 16; quadparam[1] = 4;
	//int nangle_oct = CalNangleOct(AngQuadType::UNIFORM, quadparam);
	int nangle_oct = CalNangleOct(AngQuadType::Bi3, quadparam);

	BndFlux.Initialize(XSLib.GetNumGroup(), nangle_oct);
	FSR.Initialize(XSLib.GetNumGroup(), nangle_oct);
	FXR.InitializeMacro(XSLib.GetNumGroup(), XSLib.GetScatOrder(), false);
	FXR.SetMacroXS(XSLib, GeoHandle.GetMatIds());

	//DriverSN SnDriver(XSLib.GetNumGroup(), AngQuadType::UNIFORM, quadparam);
	DriverSN SnDriver(XSLib.GetNumGroup(), AngQuadType::Bi3, quadparam);
	SnDriver.Initialize(RaySP, BndFlux, FSR, FXR);
	//SnDriver.RaySweep();

	SteadyStateDriver SSDriver(Dimension::TwoD, XSLib.GetNumGroup(), nangle_oct, XSLib.GetScatOrder());
	SSDriver.Initialize(RaySP, BndFlux, FSR, FXR, SnDriver);
	SSDriver.SetCriterionPsi(1.e-5);
	SSDriver.SetCriterionK(5.e-8);
	SSDriver.SetMaxIter(10000);

	SSDriver.RunSS();
}
}