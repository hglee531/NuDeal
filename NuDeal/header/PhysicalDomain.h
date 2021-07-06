#pragma once
#include "Defines.h"
#include "UnitGeo.h"
#include "Library.h"
#include "GeoHandle.h"
#include "Array.h"
#include "PhysicalComposition.h"
#include "Math.hpp"

namespace PhysicalDomain {

template <typename T> using Array = LinPack::Array_t<T>;

enum dir6 {
	xleft, xright,
	yleft, yright,
	zleft, zright
};

struct ConnectInfo_t {
	array<int, 6> NeighLv, NeighId, NeighInodes; // 0~6 : xl, xr, yl, yr, zl, zr
	int thisLv, thisId;
};

struct CompileInfo_t {
	int idvol;
	vector<int> inodes;
	vector<double> weightnodes;
	double volsum;
};

class BaseDomain {
public:
	using Dimension = Geometry::Dimension;
	using NodeInfo_t = Geometry::NodeInfo_t;
	using GeometryHandler = Geometry::GeometryHandler;

protected:
	bool isalloc = false;
	Dimension mode;
	double x0, y0, z0;
	double lx0, ly0, lz0;

	int Nx, Ny, Nz;
	int nx, ny, nz;
	int nnode, divlevel;
	vector<int> nnodeLv;
	vector<vector<int>> upperdivmap, lowerdivmap, serialdivmap;
	vector<int> innodeLv;

	vector<NodeInfo_t> nodeInfo;
public:
	BaseDomain() {};

	BaseDomain(const GeometryHandler &rhs) { Create(rhs); }

	void Create(const GeometryHandler &rhs);

	Dimension GetDimension() const { return mode; }

	void GetBaseSizes (int3 &Nxyz, int3 &nxyz, int &nnode, int &divlevel) const {
		Nxyz.x = Nx; Nxyz.y = Ny; Nxyz.z = Nz; nxyz.x = nx; nxyz.y = ny; nxyz.z = nz;
		nnode = this->nnode; divlevel = this->divlevel;
	}

	double3 GetNodeSizes() const { return make_double3(lx0, ly0, lz0); }

	int3 GetNnodeZero() const { return make_int3(Nx, Ny, Nz); }

	const auto& GetNnodeLv() const { return nnodeLv; }
	const auto& GetUpperdivmap() const { return upperdivmap; }
	const auto& GetLowerdivmap() const { return lowerdivmap; }
	const auto& GetSerialdivmap() const { return serialdivmap; }
	const auto& GetBaseNodeInfo() const { return nodeInfo; }
	const auto& GetInNodeLv() const { return innodeLv; }
};

class ConnectedDomain : public BaseDomain {
protected:
	vector<ConnectInfo_t> connectInfo;

private:
	bool FindNeighbor(dir6 srchdir, array<int,3> ixyz, array<int,2> LvnId, array<int,2> &NeighLvnId);

	void RecursiveConnect(int &ActiveLv, int thisidonset);

public:
	ConnectedDomain() {}

	ConnectedDomain(const GeometryHandler &rhs) : BaseDomain(rhs) { CreateConnectInfo(); };

	void CreateConnectInfo();

	void PrintConnectInfo() const;

	const auto& GetConnectInfo() { return connectInfo; }
};

class RayTracingDomain : public BaseDomain {
public:
	template <typename T> using Array = LinPack::Array_t<T>;

private:
	int edgeid[8];
	int3 dxyz[8];
	int3 range0[8][2];
	int3 range[8][2];

	int nloop[8];

protected:
	Array<int> sweepnodes;
	Array<int> sweepIdInBox;
	
	Array<int2> sweeploop; // x : beginning index in sweepnodes, y : storage direction

private:
	void RecursiveTrace(array<queue<int2>, 8> &Qsweeploop, int Lv, int id, int &iternode, int Oct);

	void SetTraceOrder(array<queue<int2>, 8> &Qsweeploop);

public:
	RayTracingDomain() {}

	RayTracingDomain(const GeometryHandler &rhs) : BaseDomain(rhs) { CreateTracingDomain(); }

	void CreateTracingDomain();

	void PrintTracingInfo();

	const auto& GetSweepNodes() const { return sweepnodes; }
	const auto& GetSweepIdInBox() const { return sweepIdInBox; }
	const auto& GetSweepLoop() const { return sweeploop; }
};

/*
	When compiling sub-domains, the different volumes in a block are diffentiated.
	However, when the same volume appears more than once segregated by another volume between it,
	  they should be split.
	For now, it is not capable of spliting those volumes.
*/
class CompiledDomain {
public:
	using Dimension = Geometry::Dimension;
	using NodeInfo_t = Geometry::NodeInfo_t;

protected:
	int Bx, By, Bz;
	int3 block;
	int nblocks;

	vector<int> compilemap, compileid;
	vector<double> compileW;
	vector<CompileInfo_t> compileInfo;

	void Decompile();

public:
	CompiledDomain(int3 block, const BaseDomain &rhs);

	CompiledDomain(int3 block, const CompiledDomain &rhs);

	void GetCompileSizes(int &Bx, int &By, int &Bz, int3 &block, int &nblocks) const 
	{ Bx = this->Bx; By = this->By; Bz = this->Bz; block = this->block; nblocks = this->nblocks; }

	int GetNblocks() const { return nblocks; }
	const auto& GetCompileMap() const { return compilemap; }
	const auto& GetCompileInfo() const { return compileInfo; }	
	const auto& GetDecompileId() const { return compileid; }
	const auto& GetDecompileWeights() const { return compileW; }

	void PrintCompileInfo(string filename = "Compile.out") const;
};

class FluxBoundary {
private:
	int ng, nangle_oct;
	int Nx, Ny, Nz;
	Array<double> xbndflux, ybndflux, zbndflux;
public:
	FluxBoundary(int Nx, int Ny, int Nz) { this->Nx = Nx; this->Ny = Ny; this->Nz = Nz; }

	void Initialize(int ng, int nangle_oct);

	void InitBndFluxUnity();

	void VaccumFlux();

	void UpdateBoundary(int leftcond, int rightcond);

	auto& GetXBndFlux() { return xbndflux; }
	auto& GetYBndFlux() { return ybndflux; }
	auto& GetZBndFlux() { return zbndflux; }
};

class FlatXSDomain : public CompiledDomain {
public:
	using XSLibrary = Library::XSLibrary_t;

private:
	bool isMicro;
	bool isTHfeed;

	int ng, scatorder;

	Array<int> idiso;
	Array<double> pnum, temperature;
	Array<double> xst, xsnf, chi, xskf;
	Array<double> xssm;

public:
	FlatXSDomain(int3 block, const CompiledDomain &rhs) : CompiledDomain(block, rhs) {}

	void InitializeMacro(int ng, int scatorder, bool isTHfeed);

	void Initialize(int ng, int scatorder, int niso);

	void SetMacroXS(const XSLibrary &XS, const vector<int> imat);

	int GetScatOrder() const { return scatorder; }
	const auto& GetTotalXS() const { return xst; }
	const auto& GetNuFisXS() const { return xsnf; }
	const auto& GetFisSpectrum() const { return chi; }
	const auto& GetScatMatrix() const { return xssm; }
};

class FlatSrcDomain : public CompiledDomain {
	// for neutronics
private:
	int ng, nangle_oct, scatorder;
	bool isEx;

	Array<double> flux, psi;
	Array<double> srcS, srcEx;

	Array<double> src;

	void AccumulateSource(const vector<int> &idFXR, const Array<double> &chi);

public:
	FlatSrcDomain(int3 block, const BaseDomain &rhs) : CompiledDomain(block, rhs) {}

	void Initialize(int ng, int nangle_octant, int scatorder = 0, bool isEx = false);

	void InitFluxUnity() { std::fill(flux.begin(), flux.end(), Math::PI*4.0); }

	void InitFluxZero() { std::fill(flux.begin(), flux.end(), 0.0); }

	double UpdateFisSource(const FlatXSDomain &FXR);

	void UpdateScatSource(const FlatXSDomain &FXR);

	void UpdateAngularSource(const FlatXSDomain &FXR);

	void NormalizePsi(double keff);

	auto& GetScalarFlux() { return flux; }
	const auto& GetFisSrc() const { return psi; }
	const auto& GetScatSrc() const { return srcS; }
	const auto& GetExternalSrc() const { return srcEx; }
	const auto& GetTotalSrc() const { return src; }
};

inline void DebugPhysicalDomain() {

		using SurfType = Geometry::UnitSurf::SurfType;
		using CartPlane = Geometry::CartPlane;
		using UnitSurf = Geometry::UnitSurf;
		using UnitVol = Geometry::UnitVol;
		using UnitComp = Geometry::UnitComp;
		using GeometryHandler = Geometry::GeometryHandler;

		double c_cir[3] = { 0.0, 0.0, 0.54 };
		UnitSurf Circle(SurfType::CIRCLE, c_cir, CartPlane::XY);

		double c_zpln0[2] = { -1.0, 0.0 };
		double c_zpln1[2] = { 1.0, 3.0 };
		UnitSurf Zpln0(SurfType::ZPLN, c_zpln0), Zpln1(SurfType::ZPLN, c_zpln1);

		double c_xpln0[2] = { -1.0, -0.63 };
		double c_xpln1[2] = { 1.0, 0.63 };
		double c_ypln0[2] = { -1.0, -0.63 };
		double c_ypln1[2] = { 1.0, 0.63 };
		UnitSurf Xpln0(SurfType::XPLN, c_xpln0), Xpln1(SurfType::XPLN, c_xpln1);
		UnitSurf Ypln0(SurfType::YPLN, c_ypln0), Ypln1(SurfType::YPLN, c_ypln1);

		UnitVol Cylinder;
		Cylinder.Append(Circle);
		Cylinder.Append(Zpln0);
		Cylinder.Append(Zpln1);
		Cylinder.Finalize();

		UnitVol Box;
		Box.Append(Zpln0); Box.Append(Zpln1);
		Box.Append(Ypln0); Box.Append(Ypln1);
		Box.Append(Xpln0); Box.Append(Xpln1);
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

		UnitComp BoxComp(0, Box), CylComp(1, Cylinder);
		CylComp.Finalize();
		GeometryHandler GeoHandle;
		double origin[3] = { -0.63,-0.63,0 }, L[3] = { 1.26, 1.26,3 };
		GeoHandle.SetOrdinates(origin, L);
		GeoHandle.append(BoxComp);
		GeoHandle.append(CylComp);
		GeoHandle.FinalizeVolumes();
		GeoHandle.Discretize(Geometry::Dimension::TwoD, 0.05, 0.005, 6);

		RayTracingDomain RaySP(GeoHandle);
		int3 blockFSR, blockFXR;
		blockFSR.x = 1; blockFSR.y = 1; blockFSR.z = 1;
		blockFXR.x = 2; blockFXR.y = 2; blockFXR.z = 1;
		FlatSrcDomain FSR(blockFSR, RaySP);
		FlatXSDomain FXR(blockFXR, FSR);
	}
};