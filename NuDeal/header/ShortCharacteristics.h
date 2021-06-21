#pragma once
#include "Defines.h"
#include "Array.h"
#include "PhysicalDomain.h"

namespace Transport {
enum class AngQuadType{
	UNIFORM, // Uniform distribution on (theta, phi) space
	LS, // Level symmetric
	GC, // Gauss-Chebyshev
	Bi3, // Bickley-3 on polar, Uniform on azimuthal
	G_Bi3, // Bickley-3 on polar, Gaussian quadrature on aximuthal
	QR // Quadruple Range
};

class AngularQuadrature {
public:
	AngQuadType quadtype;
	int nangle_oct;
	vector<double> weights;
	vector<double3> omega;
public:
	void CreateSet(AngQuadType type, vector<int> parameters);

	AngQuadType GetType() { return quadtype; }
	int GetNanglesOct() { return nangle_oct; }
	const auto& GetWeights() { return weights; }
	const auto& GetOmega() { return omega; }
};

class AsymptoticExp {
private:
	double tol;
	double2 range;
	vector<double> PtsExp;
	vector<int> PtsXval;
public:
	AsymptoticExp(double tol, double xL, double xR);

	double ExpFast(double xval);

	double ExpSafe(double xval);
};

class DriverSCMOC{
	template <typename T> using Array = LinPack::Array_t<T>;
	using RaySegment = PhysicalDomain::RaySegmentDomain;
	using FlatXS = PhysicalDomain::FlatXSDomain;
	using FlatSrc = PhysicalDomain::FlatSrcDomain;
	using Dimension = Geometry::Dimension;
	using CartAxis = Geometry::CartAxis;
	using NodeInfo_t = Geometry::NodeInfo_t;
private:
	Dimension mode;
	AngularQuadrature QuadSet;
	AsymptoticExp Exponent;

	int ng, scatorder;

	int3 Nxyz, nxyz;
	double3 lxyz0;
	int nnode;
	int divlevel;
	const NodeInfo_t* nodeinfo;
	const int* innodeLv;

	Array<double> trackL; // trackL[nangle_oct]
	double3 wtIF;
	//vector<double> wtIF; // wtIF[3]
	//Array<double> wtIF; // weights from interfaces along x,y,z axes, wtIF[nangle_oct][3]

	double *scalarflux;
	const double *srcF, *srcS, *srcEx;
	const double *xst;

	int nFXR, nFSR;
	vector<int> idFXR;
	const int *idFSR;
	const double *wFSR;

private:
	void AvgInFlux(double *outx, double *outy, double *outz, double *avgin);

	void TrackAnode(int inode, double *angflux, double *src);

public:
	void Initialize(RaySegment& rhs, FlatSrc& FSR, const FlatXS& FXR);

	DriverSCMOC(int ng, AngQuadType quadtype, vector<int> quadparameter);
};
}