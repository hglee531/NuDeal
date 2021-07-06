#pragma once
#include "Defines.h"
#include "UnitGeo.h"

/*
  The boundary conditions are not treated.
	It is so complicating works for me and for now.
	If sometime it is applied,
	(1) Boundary surface definition in UnitGeo.h
	(2) B.S search for zero-level nodes with the normal vector calculation
	(3) Weights for reflective-B.C integration in AngularQuadrature class
	(4) Update boundary fluxes in xbndflux, ybndflux, zbndflux in BoundarySurface class
*/

namespace Geometry {

constexpr double mintau_default = 0.01;
constexpr double maxtau_default = 0.0001;
constexpr double eps_dis = 1.e-10;

enum class Dimension {
	ThreeD,
	TwoD,
	OneD
};

struct NodeInfo_t {
	double3 midpt;
	//vector<int> idvols;
	//vector<double> vol;
	int idvol;
	double volcm3;
};

class splittree {
private :
	splittree *subnode = nullptr, *uppernode = nullptr;
	int id, nleaf;
	NodeInfo_t thisinfo;

public:
	splittree() {}

	splittree(int id) { Assign(id); }

	splittree(int id, splittree &uppernode) { Assign(id, uppernode); }

	~splittree() { if (subnode != nullptr) delete[]subnode; }

	void Assign(int id);

	void Assign(int id, splittree &uppernode);

	void Branching(int nleaf);

	void RecordNodeInfo(int idvol, double volcm3, double3 midpt);

	splittree* GetPtrSubnode() { return subnode; }

	int GetId() { return id; }

	int GetNleaf() { return nleaf; }

	int GetUpperId() { return GetId(); }

	NodeInfo_t GetNodeInfo() { return thisinfo; }
};

class GeometryHandler {
private:
	bool issetord, isfinal;
	double x0, y0, z0, Lx, Ly, Lz;

	double mintau, maxtau;
	Dimension mode;

	int nvol, nnode; // nnodeLv : Number of nodes at each level, size : [divlevel]
	int Nx, Ny, Nz; // upperdivmap : Node indices of one-level higher nodes, size [divlevel][nnodeLv[ir]]
	int divlevel;
	double lx0, ly0, lz0;
	vector<int> nnodeLv;
	vector<vector<int>> upperdivmap;
	vector<vector<NodeInfo_t>> info;
	
	vector<UnitVol> Volumes;
	queue<UnitVol> UVbuf;

	vector<int> imat;
	queue<int> matidbuf;

	vector<double3> BdL, BdR;
	vector<double> pow3;

	void init();

	int FindVolId(double3 ptL, double3 ptR, bool lowest);

	void RecursiveSplit(double3 ptL, double3 ptR, int thisLv, splittree &thisnode);

	void LevelCount(int thisLv, splittree &thisnode);

	void RecordDiscInfo(int thisLv, splittree &thisnode);
public:
	GeometryHandler() { init(); }

	GeometryHandler(double origin[3], double L[3]);
	
	void SetOrdinates(double origin[3], double L[3]);

	void append(UnitVol &avol, int matid) { UVbuf.push(avol); matidbuf.push(matid); }

	void append(UnitComp &acomp);

	GeometryHandler& operator<<(UnitComp &acomp) { append(acomp); return *this; }

	void Create(int ncomp, UnitComp *Comps) { for (int i = 0; i < ncomp; i++) append(Comps[i]); }

	void FinalizeVolumes();

	bool Discretize(Dimension Dim, double minlen = mintau_default, double maxlen = maxtau_default , int groundLv = 0);

	void PrintDiscInfo() const;

	Dimension GetDimension() const { return mode; }

	void GetSizeScalars(int3 &N, int &nnode, int &divlevel, double3 &origin, double3 &Width0) const;

	const auto& GetNnodeLv() const { return nnodeLv; }

	const auto& GetUpperdivmap() const { return upperdivmap; }

	const auto& GetNodeinfo() const { return info; }

	const auto& GetUnitVolumes() const { return Volumes; }
	
	const auto& GetMatIds() const { return imat; }
};


inline void DebugGeomHandle() {

	using SurfType = UnitSurf::SurfType;

	double c_cir[3] = { 0.0, 0.0, 0.54 };
	UnitSurf Circle(SurfType::CIRCLE, c_cir, CartPlane::XY);

	double c_zpln0[2] = { -1.0, 0.0 };
	double c_zpln1[2] = { 1.0, 3.0 };
	UnitSurf Zpln0(SurfType::ZPLN, c_zpln0), Zpln1(SurfType::ZPLN, c_zpln1);

	double c_xpln0[2] = { -1.0, -0.63 };
	double c_xpln1[2] = { 1.0, 1.89 };
	double c_ypln0[2] = { -1.0, -0.63 };
	double c_ypln1[2] = { 1.0, 1.89 };
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
	UnitComp CylComp2(1, Cylinder), CylComp3(1, Cylinder), CylComp4(1, Cylinder);
	CylComp.Finalize(); CylComp2.Finalize(); CylComp3.Finalize(); CylComp4.Finalize();
//	CylComp.Relocate(-0.63, 0.63, 0.);	CylComp2.Relocate(-0.63, -0.63, 0.);
//	CylComp3.Relocate(0.63, 0.63, 0.);	CylComp4.Relocate(0.63, -0.63, 0.);
	CylComp2.Relocate(1.26, 0., 0.);
	CylComp3.Relocate(1.26, 1.26, 0.);	CylComp4.Relocate(0., 1.26, 0.);
	GeometryHandler GeoHandle;
	double origin[3] = { -0.63,-0.63,0 }, L[3] = { 2.52, 2.52,3 };
	GeoHandle.SetOrdinates(origin, L);
	GeoHandle.append(BoxComp);
	GeoHandle.append(CylComp);
	GeoHandle.append(CylComp2);
	GeoHandle.append(CylComp3);
	GeoHandle.append(CylComp4);
	GeoHandle.FinalizeVolumes();
	GeoHandle.Discretize(Geometry::Dimension::TwoD, 0.2, 0.03);
}

}