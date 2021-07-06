#pragma once
#include "Defines.h"

namespace Geometry
{

constexpr double eps_geo = 1.e-10;

enum class CartAxis {
	X,
	Y,
	Z
};

enum class CartPlane {
	XY,
	YZ,
	XZ
};

class UnitSurf {

public :
	
	// Aliasing & Enumerator

	using Equation_t = array<double, 10>;

	enum class TransformType {
		RELOC,
		ROTAT,
		SCALE,
		SHAER
	};

	enum class SurfType {
		XPLN,
		YPLN,
		ZPLN,
		PLN,
		CIRCLE,
		ELLIPSE,
		SPHERE,
		ELLIPSOID,
		GENERAL
	};

private:

	// The realm defined with a surface would be,
	// (c_xs*x+c_x)*x + (c_ys*y+c_y)*y + (c_zs*z+c_z)*z - c < 0
	bool is_curve = false;
	string name;
	Equation_t eq = { 0.0, };

	double& c_xs = eq[0];
	double& c_ys = eq[1];
	double& c_zs = eq[2];
	double& c_xy = eq[3];
	double& c_yz = eq[4];
	double& c_xz = eq[5];
	double& c_x = eq[6];
	double& c_y = eq[7];
	double& c_z = eq[8];
	double& c = eq[9];

	static Equation_t CalCoeffs(const Equation_t& in, double xyzc[4][4]);
	static int GetLocalExCurve(int axis, double coeff[6], double sol[2][2]);

public:

	static bool GetTriplePoint(double coeff0[4], double coeff1[4], double coeff2[4], double sol[3]);

private:

	void Transform(double invTM[4][4]);
	void SetCurve();
	bool GetPointVar2(CartAxis axis, double var[2], double &sol);

public:

	UnitSurf() {}
	UnitSurf(const Equation_t& eq) 
	{ Create(eq); }
	UnitSurf(SurfType surftype, double *_coeff, CartPlane cartesianPln = CartPlane::XY)
	{ Create(surftype, _coeff, cartesianPln); }
	UnitSurf(const UnitSurf &rhs) 
	{ Create(rhs);	}

	void Create(SurfType surftype, double *_coeff, CartPlane cartesianPln = CartPlane::XY);
	void Create(const Equation_t& eq, const string& name = "") 
	{ this->eq = eq; SetCurve(); }
	void Create(const UnitSurf& rhs) { eq = rhs.eq; is_curve = rhs.is_curve; }

	const Equation_t& GetEquation() const { return eq; }
	bool IsCurve() { return is_curve; }

	void Relocate(double dx, double dy, double dz);
	void Relocate(double3 d) { Relocate(d.x, d.y, d.z); }
	void Rotate(double cos, double sin, CartAxis Ax);
	void Rotate(double2 c, CartAxis Ax) { Rotate(c.x, c.y, Ax); }

	int GetIntersection(CartPlane cartPlane, const array<double, 2>& val, array<double, 2>& sol);
	bool IsInside(double x, double y, double z, bool includeOn = false);
	int GetLocalExSelf(double so1[6][3]);
	int GetLocalExPln(double CoeffPln[4], double sol[6][3]);

	UnitSurf& operator=(const UnitSurf& rhs) { *this = UnitSurf(rhs); return *this; }
};

class UnitVol {

public:
	
	using Transform = UnitSurf::TransformType;
	using SurfType = UnitSurf::SurfType;

private:

	string name;
	vector<UnitSurf> Surfaces;
	bool isbounded = false;
	double2 xlr, ylr, zlr;
	double vol;

	int OneIntersection(int idx, CartPlane CartPlane, const array<double, 2>& val, array<double, 2>& sol);
	void ResidentFilter(int &codes, int acode, double localEx[6][3], double corners[6][3]);
	bool CalBoundBox();
	double CalVolume(int nx = 300, int ny = 300, int nz = 300);

public:

	UnitVol() {}
	UnitVol(const string& name) { this->name = name; }
	UnitVol(int nsurf, const UnitSurf *_Surfaces) { Create(nsurf, _Surfaces); }
	UnitVol(const UnitSurf &_Surface) { Create(_Surface); }
	UnitVol(const UnitVol& rhs) { Create(rhs); }

	void Create(int _nsurf, const UnitSurf *_Surfaces);
	void Create(const vector<UnitSurf>& Surfaces);
	void Create(const UnitSurf& _Surfaces);
	void Create(const UnitVol& rhs);

	bool IsAlloc() const { return !Surfaces.empty(); }
	const auto& GetSurfaces() const { return Surfaces; }
	void Append(const UnitSurf &asurf) 
	{ Surfaces.emplace_back(asurf); isbounded = false; }
	int GetNumSurfaces() const { return Surfaces.size(); }
	double GetVolume() const { return vol; }
	bool GetBoundBox (double2& xlr, double2& ylr, double2& zlr) const;

	void Relocate(double dx, double dy, double dz);
	void Relocate(double3 d) { Relocate(d.x, d.y, d.z); }
	void Rotate(double cos, double sin, CartAxis Ax);
	void Rotate(double2 c, CartAxis Ax) { Rotate(c.x, c.y, Ax); }
	void Translation(deque<int> move);

	bool IsInside(double x, double y, double z, bool includeOn = false);
	int GetIntersection(CartPlane CartPlane, const array<double, 2>& val, vector<double>& sol);
	void Finalize();

	UnitVol& operator=(const UnitVol& rhs) {
		this->Create(rhs.GetSurfaces());
		rhs.GetBoundBox(xlr, ylr, zlr);
		this->vol = rhs.GetVolume();
		return *this;
	}

	UnitVol& operator<<(const UnitSurf &asurf) { Append(asurf); return *this; }
};

inline void DebugUnitGeo() {
	
	using SurfType = UnitSurf::SurfType;

	double c_cir[3] = { 0.0, 0.0, 0.54 };
	UnitSurf Circle(SurfType::CIRCLE, c_cir, CartPlane::XY);

	double c_zpln0[2] = { -1.0, 0.0 };
	double c_zpln1[2] = { 1.0, 3.0 };
	UnitSurf Zpln0(SurfType::ZPLN, c_zpln0), Zpln1(SurfType::ZPLN, c_zpln1);

	double c_sphere[10] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -9.0 };
	UnitSurf Sphere(SurfType::GENERAL, c_sphere);

	UnitVol Svol(Sphere);
	UnitVol Cylinder;
	Cylinder.Append(Circle);
	Cylinder.Append(Zpln0);
	Cylinder.Append(Zpln1);
	Svol.Finalize();
	Cylinder.Finalize();

	bool isbounded;
	double2 xs, ys, zs;
	isbounded = Svol.GetBoundBox(xs, ys, zs);
	cout << "Bound info, (" << isbounded << ")" << endl;
	if (isbounded) cout << "X : [" << xs.x << "," << xs.y << "] ";
	if (isbounded) cout << "Y : [" << ys.x << "," << ys.y << "] ";
	if (isbounded) cout << "Z : [" << zs.x << "," << zs.y << "] " << endl;
	if (isbounded) cout << "Vol = " << Svol.GetVolume() << "cm^3" << endl;

	isbounded = Cylinder.GetBoundBox(xs, ys, zs);
	cout << "Bound info, (" << isbounded << ")" << endl;
	if (isbounded) cout << "X : [" << xs.x << "," << xs.y << "] ";
	if (isbounded) cout << "Y : [" << ys.x << "," << ys.y << "] ";
	if (isbounded) cout << "Z : [" << zs.x << "," << zs.y << "] " << endl;
	if (isbounded) cout << "Vol = " << Cylinder.GetVolume() << "cm^3" << endl;
}

class UnitComp : public UnitVol{
private:

	int imat;

public:
	
	UnitComp(int imat, const UnitVol &Volume) { Create(imat, Volume); }
	UnitComp(const UnitComp& rhs) { Create(rhs); }

	void Create(int nvol, const UnitVol &Volume);
	void Create(const UnitComp& rhs);

	int GetMatId() const { return imat; }
};

}