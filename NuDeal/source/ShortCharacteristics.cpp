#include "ShortCharacteristics.h"
#include "Math.hpp"

const double w_Bi3[5][5] = {
	{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000 },
	{ 0.218863, 0.781137, 0.000000, 0.000000, 0.000000 },
	{ 0.044355, 0.270249, 0.685397, 0.000000, 0.000000 },
	{ 0.010165, 0.070812, 0.289244, 0.629779, 0.000000 },
	{ 0.002782, 0.019393, 0.089318, 0.298508, 0.590000 }
};
const double sin_Bi3[5][5] = {
	{ 0.707106, 0.000000, 0.000000, 0.000000, 0.000000 },
	{ 0.373917, 0.901601, 0.000000, 0.000000, 0.000000 },
	{ 0.166542, 0.523932, 0.928619, 0.000000, 0.000000 },
	{ 0.079818, 0.266013, 0.597389, 0.941545, 0.000000 },
	{ 0.041631, 0.140115, 0.331811, 0.645183, 0.949724 }
};

namespace Transport {
void AngularQuadrature::CreateSet(AngQuadType type, vector<int> parameters){
	using Math::PI;
	quadtype = type;
	switch (type)
	{
	case Transport::AngQuadType::UNIFORM:
	{
		int nazi = parameters[0], npolar = parameters[1];

		nangle_oct = CalNangleOct(type, parameters);

		double3 thisomega;
		double dtheta = PI * 0.5 / nazi, dphi = PI * 0.5 / npolar;
		for (int ip = 0; ip < npolar; ip++) {
			double phi = dphi * 0.5 * static_cast<double>(1 + 2 * ip);
			thisomega.z = cos(phi);
			for (int ia = 0; ia < nazi; ia++) {
				double theta = dtheta * 0.5 * static_cast<double>(1 + 2 * ia);
				thisomega.x = sin(phi) * cos(theta); thisomega.y = sin(phi) * sin(theta);
				omega.push_back(thisomega);
				weights.push_back(dtheta * (cos(phi - dphi * 0.5) - cos(phi + dphi * 0.5)));
			}
		}

		break;
	}
	case Transport::AngQuadType::LS: 
	{
		break;
	}
	case Transport::AngQuadType::GC:
		break;
	case Transport::AngQuadType::Bi3:
	{
		int nazi = parameters[0], npolar = parameters[1];
		int iBi3 = npolar - 1;

		nangle_oct = CalNangleOct(type, parameters);

		double3 thisomega;
		double dtheta = PI * 0.5 / nazi;
		for (int ip = 0; ip < npolar; ip++) {
			double sinphi = sin_Bi3[iBi3][ip];
			thisomega.z = sqrt(1. - sinphi * sinphi);
			for (int ia = 0; ia < nazi; ia++) {
				double theta = dtheta * 0.5 * static_cast<double>(1 + 2 * ia);
				if (nazi == 1) theta = dtheta * 0.45;
				thisomega.x = sinphi * cos(theta); thisomega.y = sinphi * sin(theta);
				omega.push_back(thisomega);
				weights.push_back(dtheta * w_Bi3[iBi3][ip]);
			}
		}

		break;
	}
	case Transport::AngQuadType::G_Bi3:
		break;
	case Transport::AngQuadType::QR:
		break;
	default:
		break;
	}
}

AsymptoticExp::AsymptoticExp(double tol, double xL, double xR) {
	this->tol = tol; range.x = xL; range.y = xR;
	int npts = (range.y - range.x) / tol + 1;
	double dx = (range.y - range.x) / static_cast<double>(npts);
	for (int i = 0; i < npts + 1; i++) {
		double pointXval = range.x + dx * static_cast<double>(i);
		double pointExp = exp(pointXval);
		PtsXval.push_back(1. - pointXval); PtsExp.push_back(pointExp);
	}
}

double AsymptoticExp::ExpFast(double xval) {
	int i = (xval - range.x) / tol;
	double pointExp = PtsExp[i], pointXval = PtsXval[i];
	return (pointExp * (xval + pointXval));
}

double AsymptoticExp::ExpSafe(double xval) {
	if ((xval - range.x) * (xval - range.y) > 0) return expf(xval);
	int i = (xval - range.x) / tol;
	double pointExp = PtsExp[i], pointXval = PtsXval[i];
	return (pointExp * (xval - pointXval));
}

void DriverSN::Initialize(RaySegment& Rays, Boundary &BndFlux, FlatSrc& FSR, const FlatXS& FXR) {
	constexpr double r3  = 0.33333333333333333333333333;
	constexpr double r6  = 0.16666666666666666666666667;
	constexpr double r12 = 0.08333333333333333333333333;
	
	mode = Rays.GetDimension();
	isY = (mode != Dimension::OneD);
	isZ = (mode == Dimension::ThreeD);

	bandrank = 1;
	if (isY) bandrank++;
	if (isZ) bandrank++;

	Rays.GetBaseSizes(Nxyz, nxyz, nnode, divlevel);
	lxyz0 = Rays.GetNodeSizes();
	trackL.Create(QuadSet.nangle_oct); //wtIF.resize(3);

	int nx = nxyz.x, ny = nxyz.y, nz = nxyz.z;
	int onbot = nx * ny * (nz - 1);
	edgemod[0] = 0; edgemod[1] = nx; edgemod[2] = nx * (ny - 1); edgemod[3] = nx * ny - 1;
	edgemod[4] = onbot; edgemod[5] = onbot + nx; edgemod[6] = onbot + nx * (ny - 1); edgemod[7] = nx * ny*nz - 1;

	bool graytable[8][3]{
	{true, true, true},
	{false, true, true},
	{true, false, true},
	{false, false, true},
	{true, true, false},
	{false, true, false},
	{true, false, false},
	{false, false, false}
	};

	for (int Oct = 0; Oct < 8; Oct++) {
		edgeidx0[Oct].x = (graytable[Oct][0]) ? 0 : Nxyz.x - 1;
		edgeidx0[Oct].y = (graytable[Oct][1]) ? 0 : Nxyz.y - 1;
		edgeidx0[Oct].z = (graytable[Oct][2]) ? 0 : Nxyz.z - 1;
		edgeidx[Oct].x = (graytable[Oct][0]) ? 0 : nxyz.x - 1;
		edgeidx[Oct].y = (graytable[Oct][1]) ? 0 : nxyz.y - 1;
		edgeidx[Oct].z = (graytable[Oct][2]) ? 0 : nxyz.z - 1;
		dxyz[Oct].x = (graytable[Oct][0]) ? 1 : -1;
		dxyz[Oct].y = (graytable[Oct][1]) ? 1 : -1;
		dxyz[Oct].z = (graytable[Oct][2]) ? 1 : -1;
	}

	rnxyz.x = 1. / static_cast<double>(ny * nz);
	rnxyz.y = 1. / static_cast<double>(nx * nz);
	rnxyz.z = 1. / static_cast<double>(nx * ny);

	rnLxyz.push_back(make_double3(12./static_cast<double>(nx*lxyz0.x), 12./static_cast<double>(ny*lxyz0.y), 12./static_cast<double>(nz*lxyz0.z)));
	Lpnxyz.push_back(make_double3(static_cast<double>(lxyz0.x) / static_cast<double>(nx), static_cast<double>(lxyz0.y) / static_cast<double>(ny), static_cast<double>(lxyz0.z) / static_cast<double>(nz)));
	for (int i = 1; i < divlevel; i++) {
		rnLxyz.push_back(rnLxyz[i - 1]);
		rnLxyz[i].x *= static_cast<double>(nx);
		rnLxyz[i].y *= static_cast<double>(ny);
		rnLxyz[i].z *= static_cast<double>(nz);
		Lpnxyz.push_back(Lpnxyz[i - 1]);
		Lpnxyz[i].x /= static_cast<double>(nx);
		Lpnxyz[i].y /= static_cast<double>(ny);
		Lpnxyz[i].z /= static_cast<double>(nz);
	}

	rnsqrxyz.x = 1. / static_cast<double>(nx * nx);
	rnsqrxyz.y = 1. / static_cast<double>(ny * ny);
	rnsqrxyz.z = 1. / static_cast<double>(nz * nz);

	nodeinfo = Rays.GetBaseNodeInfo().data(); innodeLv = Rays.GetInNodeLv().data();
	sweepnodes = &Rays.GetSweepNodes(); sweeploop = &Rays.GetSweepLoop();
	scatorder = FXR.GetScatOrder();

	nloop = (sweeploop->size() / 8) - 1;
	nFXR = FXR.GetNblocks(); nFSR = FSR.GetNblocks();

	vector<double> rpow3, rpow9;

	double val = 1;
	for (int i = 0; i < divlevel; i++) { rpow3.push_back(val); rpow9.push_back(val*val); val /= 3.; }

	double Ax = lxyz0.y * lxyz0.z, Ay = lxyz0.x * lxyz0.z, Az = lxyz0.x * lxyz0.y;
	double Atot = Ax;
	if (mode != Dimension::OneD) Atot += Ay; 
	if (mode == Dimension::ThreeD) Atot += Az;

	wtIF.x = Ax / Atot; 
	wtIF.y = (mode != Dimension::OneD) ? Ay / Atot : 0.0; 
	wtIF.z = (mode == Dimension::ThreeD) ? Az / Atot : 0.0; 

	trackL.Create(2, QuadSet.nangle_oct, divlevel); wtVol.resize(QuadSet.nangle_oct);
	wtAband.Create(3, 3, 3, QuadSet.nangle_oct); std::fill(wtAband.begin(), wtAband.end(), 0.0);
	wtVband.Create(3, 3, 3, QuadSet.nangle_oct); std::fill(wtVband.begin(), wtVband.end(), 0.0);
	trackLband.Create(QuadSet.nangle_oct, divlevel);

	// So far, only for 2D
	wtbandP[0][0].Create(1, 2, QuadSet.nangle_oct); std::fill(wtbandP[0][0].begin(), wtbandP[0][0].end(), 0.0);
	wtbandP[0][1].Create(1, 2, QuadSet.nangle_oct); std::fill(wtbandP[0][1].begin(), wtbandP[0][1].end(), 0.0);
	wtbandP[1][0].Create(1, 2, QuadSet.nangle_oct); std::fill(wtbandP[1][0].begin(), wtbandP[1][0].end(), 0.0);
	wtbandP[1][1].Create(1, 2, QuadSet.nangle_oct); std::fill(wtbandP[1][1].begin(), wtbandP[1][1].end(), 0.0);
	wtbandT[0][0].Create(1, 2, QuadSet.nangle_oct); std::fill(wtbandT[0][0].begin(), wtbandT[0][0].end(), 0.0);
	wtbandT[0][1].Create(2, 2, QuadSet.nangle_oct); std::fill(wtbandT[0][1].begin(), wtbandT[0][1].end(), 0.0);
	wtbandT[1][0].Create(2, 2, QuadSet.nangle_oct); std::fill(wtbandT[1][0].begin(), wtbandT[1][0].end(), 0.0);
	wtbandT[1][1].Create(3, 2, QuadSet.nangle_oct); std::fill(wtbandT[1][1].begin(), wtbandT[1][1].end(), 0.0);

	areaIF = make_double3(Ax / lxyz0.z, Ay / lxyz0.z, Az);
	LsqrIF = make_double3(Math::square(lxyz0.x)*r12, Math::square(lxyz0.y)*r12, Math::square(lxyz0.z)*r12);

	VBand.resize(QuadSet.nangle_oct);

	for (int iang = 0; iang < QuadSet.nangle_oct; iang++) {
		double3 omega = QuadSet.omega[iang];
		double3 Ldeproj; // de-projected length
		Ldeproj.x = lxyz0.x / omega.x; Ldeproj.y = lxyz0.y / omega.y; Ldeproj.z = lxyz0.z / omega.z;

		double Lmin = Ldeproj.x;
		if (mode != Dimension::OneD) Lmin = (Lmin > Ldeproj.y) ? Ldeproj.y : Lmin;
		if (mode == Dimension::ThreeD) Lmin = (Lmin > Ldeproj.z) ? Ldeproj.z : Lmin;

		double3 Lproj = make_double3(0.0, 0.0, 0.0), Lrem = lxyz0;
		Lproj.x = Lmin * omega.x;	Lrem.x = lxyz0.x - Lproj.x;
		if (mode != Dimension::OneD) { Lproj.y = Lmin * omega.y; Lrem.y = lxyz0.y - Lproj.y; }
		if (mode == Dimension::ThreeD) { Lproj.z = Lmin * omega.z; Lrem.z = lxyz0.z - Lproj.z; }

		// trackL generate
		double A1 = 0., A2 = 0., A3 = 0.;
		A1 = Lrem.y*Lrem.z;
		if (mode != Dimension::OneD) {
			A1 += Lrem.x*Lrem.z;
			A2 += Lproj.x * Lrem.z + Lproj.y * Lrem.z;
		}
		if (mode == Dimension::ThreeD) {
			A1 += Lrem.x*Lrem.y;
			A2 += Lproj.x * Lrem.y + Lproj.y * Lrem.x + Lproj.z * (Lrem.x + Lrem.y);
			A3 = 2.0 * (Lproj.y * Lproj.z + Lproj.x * Lproj.y + Lproj.x * Lproj.z);
		}
		double Atot_VpL = A1 + A2 * 0.5 + A3 * r6;
		//trackL(iang, 0) = Lmin / Atot * (1 + A2 / 2.0 + A3 / 3.0);
		trackL(0, iang, 0) = Lmin * (A1 + A2 * 0.5 + A3 * r6) / Atot; // Length for flux attenuation
		trackL(1, iang, 0) = Lmin * (A1 + A2 * r3 + A3 * r12) / Atot_VpL; // Length for flux average

		// trackLband generate
		trackLband(iang, 0) = Lmin;

		// wtband generate
		double3 Arec = make_double3(Lrem.y*Lrem.z, Lrem.x*Lrem.z, Lrem.x*Lrem.y);
		double3 Apar[3] = {
			make_double3(0.0, Lproj.x * Lrem.z, Lproj.x * Lrem.y),
			make_double3(Lproj.y * Lrem.z, 0.0, Lproj.y * Lrem.x),
			make_double3(Lproj.z * Lrem.y, Lproj.z * Lrem.x, 0.0)	};
		double3 Atri[3] = {
			make_double3(0.0, Lproj.x * Lproj.z, Lproj.x * Lproj.y),
			make_double3(Lproj.y * Lproj.z, 0.0, Lproj.x * Lproj.y),
			make_double3(Lproj.y * Lproj.z, Lproj.x * Lproj.z, 0.0)	};

		//wtAband
		// to rectangle
		wtAband(0, 0, 0, iang) = Arec.x / Ax; wtAband(0, 1, 1, iang) = Arec.y / Ay; wtAband(0, 2, 2, iang) = Arec.z / Az;
		// to parallelogram
		wtAband(1, 1, 0, iang) = Apar[0].y / Ay; wtAband(1, 2, 0, iang) = Apar[0].z / Az;
		wtAband(1, 0, 1, iang) = Apar[1].x / Ax; wtAband(1, 2, 1, iang) = Apar[1].z / Az;
		wtAband(1, 0, 2, iang) = Apar[2].x / Ax; wtAband(1, 1, 2, iang) = Apar[2].y / Ay;
		// to triangle
		wtAband(2, 1, 0, iang) = Atri[0].y / Ay; wtAband(2, 2, 0, iang) = Atri[0].z / Az;
		wtAband(2, 0, 1, iang) = Atri[1].x / Ax; wtAband(2, 2, 1, iang) = Atri[1].z / Az;
		wtAband(2, 0, 2, iang) = Atri[2].x / Ax; wtAband(2, 1, 2, iang) = Atri[2].y / Ay;
		
		// wtVband
		double V = Atot_VpL;
		// to rectangle
		wtVband(0, 0, 0, iang) = Arec.x / V; wtVband(0, 1, 1, iang) = Arec.y / V; wtVband(0, 2, 2, iang) = Arec.z / V;
		// to parallelogram
		wtVband(1, 1, 0, iang) = Apar[0].y / V; wtVband(1, 2, 0, iang) = Apar[0].z / V;
		wtVband(1, 0, 1, iang) = Apar[1].x / V; wtVband(1, 2, 1, iang) = Apar[1].z / V;
		wtVband(1, 0, 2, iang) = Apar[2].x / V; wtVband(1, 1, 2, iang) = Apar[2].y / V;
		// to triangle
		wtVband(2, 1, 0, iang) = Atri[0].y / V; wtVband(2, 2, 0, iang) = Atri[0].z / V;
		wtVband(2, 0, 1, iang) = Atri[1].x / V; wtVband(2, 2, 1, iang) = Atri[1].z / V;
		wtVband(2, 0, 2, iang) = Atri[2].x / V; wtVband(2, 1, 2, iang) = Atri[2].y / V;
		
		wtVol[iang].x = (Arec.x + 0.5 * (Apar[1].x + Apar[2].x) + r6 * (Atri[1].x + Atri[2].x) ) / V;	
		wtVol[iang].y = (Arec.y + 0.5 * (Apar[0].y + Apar[2].y) + r6 * (Atri[0].y + Atri[2].y) ) / V;	
		wtVol[iang].z = (Arec.z + 0.5 * (Apar[0].z + Apar[1].z) + r6 * (Atri[0].z + Atri[1].z) ) / V;

		double3 del = make_double3(0.5*(Lrem.x - Lproj.x), 0.5*(Lrem.y - Lproj.y), 0.0);

		wtbandP[0][0](0, 0, iang) = Lrem.y;
		wtbandP[0][0](0, 1, iang) = Lrem.x;
		wtbandP[0][1](0, 0, iang) = 0.5*Lrem.y*Lproj.y/LsqrIF.y;
		wtbandP[0][1](0, 1, iang) = 0.5*Lrem.x*Lproj.x/LsqrIF.x;
		wtbandP[1][0](0, 0, iang) = -0.5*Lrem.y*Lproj.y;
		wtbandP[1][0](0, 1, iang) = -0.5*Lrem.x*Lproj.x;
		wtbandP[1][1](0, 0, iang) = Lrem.y*r12*(Math::square(Lrem.y) - 3.*Math::square(Lproj.y))/LsqrIF.y;
		wtbandP[1][1](0, 1, iang) = Lrem.x*r12*(Math::square(Lrem.x) - 3.*Math::square(Lproj.x))/LsqrIF.x;

		wtbandT[0][0](0, 0, iang) = Lproj.x; 
		wtbandT[0][0](0, 1, iang) = Lproj.y;
		wtbandT[0][1](0, 0, iang) = -Lproj.x*del.x/LsqrIF.x;
		wtbandT[0][1](0, 1, iang) = -Lproj.y*del.y/LsqrIF.y;
		wtbandT[0][1](1, 0, iang) = -Lproj.x*Lproj.x/LsqrIF.x;
		wtbandT[0][1](1, 1, iang) = -Lproj.y*Lproj.y/LsqrIF.y;
		wtbandT[1][0](0, 0, iang) = Lproj.x*del.y;
		wtbandT[1][0](0, 1, iang) = Lproj.y*del.x; 
		wtbandT[1][0](1, 0, iang) = Lproj.x*Lproj.y;
		wtbandT[1][0](1, 1, iang) = Lproj.y*Lproj.x;
		wtbandT[1][1](0, 0, iang) = -Lproj.x*del.x*del.y/LsqrIF.x;  
		wtbandT[1][1](0, 1, iang) = -Lproj.y*del.x*del.y/LsqrIF.y;
		wtbandT[1][1](1, 0, iang) = -Lproj.x*(del.x*Lproj.y + del.y*Lproj.x)/LsqrIF.x;  
		wtbandT[1][1](1, 1, iang) = -Lproj.y*(del.x*Lproj.y + del.y*Lproj.x)/LsqrIF.y;
		wtbandT[1][1](2, 0, iang) = -Lproj.x*2.*Lproj.x*Lproj.y/LsqrIF.x; 
		wtbandT[1][1](2, 1, iang) = -Lproj.y*2.*Lproj.x*Lproj.y/LsqrIF.y;
		
		VBand[iang] = (Lrem.x + Lrem.y) + 0.5*(Lproj.x + Lproj.y);

		for (int lv = 1; lv < divlevel; lv++) {
			trackL(0, iang, lv) = rpow3[lv] * trackL(0, iang, 0);
			trackL(1, iang, lv) = rpow3[lv] * trackL(1, iang, 0);
			trackLband(iang, lv) = rpow3[lv] * trackLband(iang, 0);
		}
	}

	scalarflux = &FSR.GetScalarFlux();
	src = &FSR.GetTotalSrc();
	xst = &FXR.GetTotalXS();
	xbndflux = &BndFlux.GetXBndFlux();
	ybndflux = &BndFlux.GetYBndFlux();
	zbndflux = &BndFlux.GetZBndFlux();

	const int *idFSR2FXR = FXR.GetDecompileId().data();
	idFSR = FSR.GetDecompileId().data(); wFSR = FSR.GetDecompileWeights().data();

	for (int i = 0; i < nnode; i++) idFXR.push_back(idFSR2FXR[idFSR[i]]);
}

DriverSN::DriverSN(int ng, AngQuadType quadtype, vector<int> quadparameter)
	: Exponent(0.0001, -40.0, 0.0) 
{
	this->ng = ng;
	QuadSet.CreateSet(quadtype, quadparameter);
	ntotGnAng = ng * QuadSet.nangle_oct;
}

void DriverSN::RaySweep() {
	for (int Oct = 0; Oct < 8; Oct++) {
		int inode;
		int thisLv = 0, thisId, iFSR, iFXR;
		double wtFSR;

		// Workspace for saving outgoing angular fluxes
		vector<double> avgfluxIF(ng), avgfluxVol(ng);
		Array<realphi> fluxInx, fluxIny, fluxInz;
		Array<realphi> fluxUtDX, fluxUtDY, fluxUtDZ;

		fluxInx.Create(ng, QuadSet.nangle_oct, divlevel);
		if (isY) fluxIny.Create(ng, QuadSet.nangle_oct, nxyz.x, divlevel);
		if (isZ) fluxInz.Create(ng, QuadSet.nangle_oct, nxyz.x*nxyz.y, divlevel);
		fluxUtDX.Create(ng, QuadSet.nangle_oct, divlevel);
		if (isY) fluxUtDY.Create(ng, QuadSet.nangle_oct, divlevel);
		if (isZ) fluxUtDZ.Create(ng, QuadSet.nangle_oct, divlevel);
		
		vector<int3> idstore(divlevel);
		idstore[0].x = edgeidx0[Oct].x;
		idstore[0].y = edgeidx0[Oct].y;
		idstore[0].z = edgeidx0[Oct].z;

		int3 ids = make_int3(0, 0, 0);

		for (int il = 0; il < nloop; il++) {
			int traceL = (*sweeploop)(il, Oct).x;
			int traceR = (*sweeploop)(il + 1, Oct).x, storedir = (*sweeploop)(il+1, Oct).y;

			realphi *flux_thisx;
			realphi *flux_thisy = NULL;
			realphi *flux_thisz = NULL;

			if (thisLv == 0) {
				flux_thisx = &(*xbndflux)(0, 0, idstore[thisLv].y + idstore[thisLv].z * Nxyz.y, Oct);
				if (isY) flux_thisy = &(*ybndflux)(0, 0, idstore[thisLv].x + idstore[thisLv].z * Nxyz.x, Oct);
				if (isZ) flux_thisz = &(*zbndflux)(0, 0, idstore[thisLv].x + idstore[thisLv].y * Nxyz.x, Oct);
			}
			else {
				flux_thisx = &fluxInx(0, 0, thisLv);
				if (isY) flux_thisy = &fluxIny(0, 0, idstore[thisLv].x, thisLv);
				if (isZ) flux_thisz = &fluxInz(0, 0, idstore[thisLv].x + idstore[thisLv].y * 3, thisLv);
			}

			realphi *angflux = flux_thisx;

			for (int itrace = traceL; itrace < traceR; itrace++) {
				int inode = (*sweepnodes)(itrace, Oct);
				int iFXR = idFXR[inode], iFSR = idFSR[inode];
				double wtFSR = wFSR[inode];
				
				for (int iangle = 0; iangle < QuadSet.nangle_oct; iangle++) {
					int stride = iangle * ng;
//#define NODALMEAN
#ifdef NODALMEAN
					AvgInFlux(iangle, flux_thisx + stride, flux_thisy + stride, flux_thisz + stride, avgfluxIF, avgfluxVol);
					TrackAnode(iangle, thisLv, iFXR, iFSR, wtFSR, angflux + stride, avgfluxIF, avgfluxVol, &(*src)(0, iangle, iFSR, Oct));
#else
					TrackBands(iangle, thisLv, iFXR, iFSR, wtFSR, 
						flux_thisx + stride, flux_thisy + stride, flux_thisz + stride, &(*src)(0, iangle, iFSR, Oct));
#endif
				}

#ifdef NODALMEAN
				UpdateFluxIn(angflux, flux_thisy, flux_thisz);
#endif

				idstore[thisLv].x += dxyz[Oct].x;
				//ids.x++;

				if (isY) flux_thisy += ntotGnAng * dxyz[Oct].x;
				if (isZ) flux_thisz += ntotGnAng * dxyz[Oct].x;
			}

			// Update outgoing fluxes of the upper level
			realphi *angfluxUp, *angfluxLow;
			if (storedir > 0) {
				ids.y++;

				if (thisLv > 0) {
					if (thisLv == 1) angfluxUp = &(*xbndflux)(0, 0, idstore[0].y + idstore[0].z * Nxyz.y, Oct);
					else angfluxUp = &fluxInx(0, 0, thisLv - 1);

					angfluxLow = &fluxInx(0, 0, thisLv);

					UpperOutFlux_X(thisLv, Oct, ids, angfluxUp, angfluxLow);

					angfluxUp = &fluxUtDX(0, 0, thisLv);
					LowerOutFlux_X(thisLv-1, Oct, ids, angfluxLow, angfluxUp);
				}
				idstore[thisLv].x = (thisLv > 0) ? edgeidx[Oct].x : edgeidx0[Oct].x;
				idstore[thisLv].y += dxyz[Oct].y;

				ids.x = 0;
			}
			if (storedir > 1) {
				if (thisLv > 0 && isY) {
					if (thisLv == 1) angfluxUp = &(*ybndflux)(0, 0, idstore[0].x + idstore[0].z * Nxyz.x, Oct);
					else angfluxUp = &fluxIny(0, 0, idstore[thisLv - 1].x, thisLv - 1);

					angfluxLow = &fluxIny(0, 0, 0, thisLv);

					UpperOutFlux_Y(thisLv, Oct, ids, angfluxUp, angfluxLow);

					angfluxUp = &fluxUtDY(0, 0, thisLv);
					LowerOutFlux_Y(thisLv-1, Oct, ids, angfluxLow, angfluxUp);
				}
				idstore[thisLv].y = (thisLv > 0) ? edgeidx[Oct].y : edgeidx0[Oct].y;
				idstore[thisLv].z += dxyz[Oct].z;

				ids.y = 0; ids.z++;
			}
			if (storedir > 2) {
				if (thisLv > 0) {
					if (isZ) {
						if (thisLv == 1) angfluxUp = &(*zbndflux)(0, 0, idstore[0].x + idstore[0].y * Nxyz.x, Oct);
						else angfluxUp = &fluxInz(0, 0, idstore[thisLv - 1].x + idstore[thisLv - 1].y * nxyz.x, thisLv - 1);

						angfluxLow = &fluxInz(0, 0, 0, thisLv);

						UpperOutFlux_Z(thisLv, Oct, ids, angfluxUp, angfluxLow);

						angfluxUp = &fluxUtDZ(0, 0, thisLv);
						LowerOutFlux_Z(thisLv-1, Oct, ids, angfluxLow, angfluxUp);
					}
					idstore[thisLv].z = (thisLv > 0) ? edgeidx[Oct].z : edgeidx0[Oct].z;
					ids.z = 0;
					thisLv--;
					idstore[thisLv].x += dxyz[Oct].x;
					ids = idstore[thisLv];
					ids.x = edgeidx[Oct].x + dxyz[Oct].x * ids.x;
					ids.y = edgeidx[Oct].y + dxyz[Oct].y * ids.y;
					ids.z = edgeidx[Oct].z + dxyz[Oct].z * ids.z;
				}
			}
			
			// Update incoming fluxes of the lower level
			if (storedir == 0) {
				int nextLv = thisLv + 1;
				realphi *angfluxLowX, *angfluxLowY = NULL, *angfluxLowZ = NULL;
				realphi *angfluxUpX, *angfluxUpY, *angfluxUpZ;

				angfluxLowX = &fluxUtDX(0, 0, nextLv);
				if (isY) angfluxLowY = &fluxUtDY(0, 0, nextLv);
				if (isZ) angfluxLowZ = &fluxUtDZ(0, 0, nextLv);

				if (thisLv == 0) {
					angfluxUpX = &(*xbndflux)(0, 0, idstore[thisLv].y + idstore[thisLv].z * Nxyz.y, Oct);
					if (isY) angfluxUpY = &(*ybndflux)(0, 0, idstore[thisLv].x + idstore[thisLv].z * Nxyz.x, Oct);
					if (isZ) angfluxUpZ = &(*zbndflux)(0, 0, idstore[thisLv].x + idstore[thisLv].y * Nxyz.x, Oct);
				}
				else {
					angfluxUpX = &fluxInx(0, 0, thisLv);
					if (isY) angfluxUpY = &fluxIny(0, 0, idstore[thisLv].x, thisLv);
					if (isZ) angfluxUpZ = &fluxInz(0, 0, idstore[thisLv].x + idstore[thisLv].y * 3, thisLv);
				}

				SaveUpperFlux(angfluxLowX, angfluxUpX);
				if (isY) SaveUpperFlux(angfluxLowY, angfluxUpY);
				if (isZ) SaveUpperFlux(angfluxLowZ, angfluxUpZ);

				angfluxUpX = angfluxLowX;	angfluxUpY = angfluxLowY;	angfluxUpZ = angfluxLowZ;
				angfluxLowX = &fluxInx(0, 0, nextLv);
				if (isY) angfluxLowY = &fluxIny(0, 0, 0, nextLv);
				if (isZ) angfluxLowZ = &fluxInz(0, 0, 0, nextLv);

				ids = make_int3(0, 0, 0);

				LowerOutFlux_X(thisLv, Oct, ids, angfluxLowX, angfluxUpX);
				if (isY) LowerOutFlux_Y(thisLv, Oct, ids, angfluxLowY, angfluxUpY);
				if (isZ) LowerOutFlux_Z(thisLv, Oct, ids, angfluxLowZ, angfluxUpZ);

				thisLv = nextLv;
				idstore[thisLv].x = edgeidx[Oct].x;
				idstore[thisLv].y = edgeidx[Oct].y;
				idstore[thisLv].z = edgeidx[Oct].z;
			}
		}

		Oct = Oct;
	}
}

//void DriverSN::RaySweepRigorous() {
//	for (int Oct = 0; Oct < 8; Oct++) {
//		int inode;
//		int thisLv = 0, thisId, iFSR, iFXR;
//		double wtFSR;
//
//		// Workspace for saving outgoing angular fluxes
//		vector<double> avgfluxIF(ng), avgfluxVol(ng);
//		vector<Array<realphi>> fluxInx(divlevel), fluxIny(divlevel), fluxInz(divlevel);
//
//		vector<int2> stridexy(divlevel);
//
//		int mx = 1, my = 1, mz = 1;
//		for (int lv = 0; lv < divlevel; lv++) {
//			stridexy[lv].x = Nxyz.x*mx; stridexy[lv].y = Nxyz.y*my;
//
//			fluxInx[lv].Create(ng, QuadSet.nangle_oct, Nxyz.y*my*Nxyz.z*mz);
//			if (isY) fluxIny[lv].Create(ng, QuadSet.nangle_oct, Nxyz.x*my*Nxyz.z*mz);
//			if (isZ) fluxInz[lv].Create(ng, QuadSet.nangle_oct, Nxyz.x*mx*Nxyz.y*my);
//
//			std::fill(fluxInx[lv].begin(), fluxInx[lv].end(), 0.0);
//			if (isY) std::fill(fluxIny[lv].begin(), fluxIny[lv].end(), 0.0);
//			if (isZ) std::fill(fluxInz[lv].begin(), fluxInz[lv].end(), 0.0);
//
//			mx *= 3;
//			if (isY) my *= 3;
//			if (isZ) mz *= 3;
//		}
//
//		std::copy(&(*xbndflux)(0, 0, 0, Oct), &(*xbndflux)(0, 0, 0, Oct + 1), fluxInx[0].begin());
//		if (isY) std::copy(&(*ybndflux)(0, 0, 0, Oct), &(*ybndflux)(0, 0, 0, Oct + 1), fluxIny[0].begin());
//		if (isZ) std::copy(&(*zbndflux)(0, 0, 0, Oct), &(*zbndflux)(0, 0, 0, Oct + 1), fluxInz[0].begin());
//
//		vector<int3> idstore(divlevel);
//		idstore[0].x = edgeidx0[Oct].x;
//		idstore[0].y = edgeidx0[Oct].y;
//		idstore[0].z = edgeidx0[Oct].z;
//		for (int lv = 1; lv < divlevel; lv++) {
//			idstore[lv].x = idstore[lv - 1].x * nxyz.x + edgeidx[Oct].x;
//			idstore[lv].y = idstore[lv - 1].y * nxyz.y + edgeidx[Oct].y;
//			idstore[lv].z = idstore[lv - 1].z * nxyz.z + edgeidx[Oct].z;
//		}
//		vector<int3> idinbox(divlevel);
//		vector<int3> upids(divlevel);
//
//		for (int il = 0; il < nloop; il++) {
//			int traceL = (*sweeploop)(il, Oct).x;
//			int traceR = (*sweeploop)(il + 1, Oct).x, storedir = (*sweeploop)(il + 1, Oct).y;
//
//			realphi *flux_thisx;
//			realphi *flux_thisy = NULL;
//			realphi *flux_thisz = NULL;
//
//			int stridex = stridexy[thisLv].x;
//			int stridey = stridexy[thisLv].y;
//
//
//			flux_thisx = &fluxInx[thisLv](0, 0, idstore[thisLv].y + idstore[thisLv].z * stridey);
//			if (isY) flux_thisy = &fluxIny[thisLv](0, 0, idstore[thisLv].x + idstore[thisLv].z * stridex);
//			if (isZ) flux_thisz = &fluxInz[thisLv](0, 0, idstore[thisLv].x + idstore[thisLv].y * stridex);
//
//			realphi *angflux = flux_thisx;
//
//			for (int itrace = traceL; itrace < traceR; itrace++) {
//				int inode = (*sweepnodes)(itrace, Oct);
//				int iFXR = idFXR[inode], iFSR = idFSR[inode];
//				double wtFSR = wFSR[inode];
//
//				for (int iangle = 0; iangle < QuadSet.nangle_oct; iangle++) {
//					int stride = iangle * ng;
//				//#define NODALMEAN
//#ifdef NODALMEAN
//					AvgInFlux(iangle, flux_thisx + stride, flux_thisy + stride, flux_thisz + stride, avgfluxIF, avgfluxVol);
//					TrackAnode(iangle, thisLv, iFXR, iFSR, wtFSR, angflux + stride, avgfluxIF, avgfluxVol, &(*src)(0, iangle, iFSR, Oct));
//#else
//					TrackBands(iangle, thisLv, iFXR, iFSR, wtFSR,
//						flux_thisx + stride, flux_thisy + stride, flux_thisz + stride, &(*src)(0, iangle, iFSR, Oct));
//#endif
//				}
//
//#ifdef NODALMEAN
//				UpdateFluxIn(angflux, flux_thisy, flux_thisz);
//#endif
//
//				int m = dxyz[Oct].x;
//				for (int lv = thisLv; lv < divlevel; lv++) {
//					idstore[lv].x += m;
//					m *= 3;
//				}
//
//				if (isY) flux_thisy += ntotGnAng*dxyz[Oct].x;
//				if (isZ) flux_thisz += ntotGnAng*dxyz[Oct].x;
//			}
//
//			// Update outgoing fluxes of the upper level
//			realphi *angfluxUp, *angfluxLow;
//#ifndef LinearFlux
//			if (storedir > 0) {
//				if (thisLv > 0) {
//					int strideyup = stridexy[thisLv-1].y, strideylow = stridexy[thisLv].y;
//					int iyup = idstore[thisLv - 1].y, iylow = idstore[thisLv].y;
//					int izup = idstore[thisLv - 1].z, izlow = idstore[thisLv].z;
//
//					angfluxUp = &fluxInx[thisLv - 1](0, 0, iyup + izup * strideyup);
//					angfluxLow = &fluxInx[thisLv](0, 0, iylow + izlow * strideylow);
//
//					for (int i = 0; i < ntotGnAng; i++) {
//						angfluxUp[i] += rnxyz.x * angfluxLow[i];
//					}
//				}
//
//				if (thisLv == 0) {
//					idstore[thisLv].x = edgeidx0[Oct].x;
//					for (int lv = thisLv + 1; lv < divlevel; lv++) idstore[lv].x = idstore[lv - 1].x * 3 + edgeidx[Oct].x;
//				}
//				else {
//					for (int lv = thisLv; lv < divlevel; lv++) idstore[lv].x = idstore[lv - 1].x * 3 + edgeidx[Oct].x;
//				}
//
//				int m = dxyz[Oct].y;
//				for (int lv = thisLv; lv < divlevel; lv++) {
//					idstore[lv].y += m;
//					m *= 3;
//				}
//			}
//			if (storedir > 1) {
//				if (thisLv > 0 && isY) {
//					int stridexup = stridexy[thisLv - 1].x, stridexlow = stridexy[thisLv].x;
//					int ixup = idstore[thisLv - 1].x, ixlow = idstore[thisLv].x;
//					int izup = idstore[thisLv - 1].z, izlow = idstore[thisLv].z;
//
//					angfluxUp = &fluxIny[thisLv - 1](0, 0, ixup + izup * stridexup);
//					angfluxLow = &fluxIny[thisLv](0, 0, ixlow + izlow * stridexlow);
//
//					for (int j = 0; j < 3; j++) {
//						for (int i = 0; i < ntotGnAng; i++) {
//							angfluxUp[i] += rnxyz.y * angfluxLow[i + j * dxyz[Oct].x * ntotGnAng];
//						}
//					}
//				}
//
//				if (thisLv == 0) {
//					idstore[thisLv].y = edgeidx0[Oct].y;
//					for (int lv = thisLv + 1; lv < divlevel; lv++) idstore[lv].y = idstore[lv - 1].y * 3 + edgeidx[Oct].y;
//				}
//				else {
//					for (int lv = thisLv; lv < divlevel; lv++) idstore[lv].y = idstore[lv - 1].y * 3 + edgeidx[Oct].y;
//				}
//
//				int m = dxyz[Oct].z;
//				for (int lv = thisLv; lv < divlevel; lv++) {
//					idstore[lv].z += m;
//					m *= 3;
//				}
//			}
//			if (storedir > 2) {
//				if (thisLv == 0) {
//					idstore[thisLv].z = edgeidx0[Oct].z;
//					for (int lv = thisLv + 1; lv < divlevel; lv++) idstore[lv].z = idstore[lv - 1].z * 3 + edgeidx[Oct].z;
//				}
//				else {
//					for (int lv = thisLv; lv < divlevel; lv++) idstore[lv].z = idstore[lv - 1].z * 3 + edgeidx[Oct].z;
//				}
//				if (thisLv > 0) {
//					if (isZ) {
//					int stridexup = stridexy[thisLv - 1].x, stridexlow = stridexy[thisLv].x;
//						int ixup = idstore[thisLv - 1].x, ixlow = idstore[thisLv].x;
//						int iyup = idstore[thisLv - 1].y, iylow = idstore[thisLv].y;
//
//						angfluxUp = &fluxInz[thisLv - 1](0, 0, ixup + iyup * stridexup);
//						angfluxLow = &fluxInz[thisLv](0, 0, ixlow + iylow * stridexlow);
//
//						for (int k = 0; k < 3; k++) {
//							for (int j = 0; j < 3; j++) {
//								for (int i = 0; i < ntotGnAng; i++) {
//									angfluxUp[i] += rnxyz.y * angfluxLow[i + (k * dxyz[Oct].y * stridexlow + j * dxyz[Oct].x) * ntotGnAng];
//								}
//							}
//						}
//					}
//					thisLv--;
//					int m = dxyz[Oct].x;
//					for (int lv = thisLv; lv < divlevel; lv++) {
//						idstore[lv].x += m;
//						m *= 3;
//					}
//				}
//			}
//#endif
//			// Update incoming fluxes of the lower level
//			if (storedir == 0) {
//				int nextLv = thisLv + 1;
//
//				int strideyup = stridexy[thisLv].y, strideylow = stridexy[nextLv].y;
//				int stridexup = stridexy[thisLv].x, stridexlow = stridexy[nextLv].x;
//				int ixup = idstore[thisLv].x, ixlow = idstore[nextLv].x;
//				int iyup = idstore[thisLv].y, iylow = idstore[nextLv].y;
//				int izup = idstore[thisLv].z, izlow = idstore[nextLv].z;
//
//				
//				angfluxUp = &fluxInx[thisLv](0, 0, iyup + izup * strideyup);
//				angfluxLow = &fluxInx[nextLv](0, 0, iylow + izlow * strideylow);
//
//
//				for (int k = 0; k < nxyz.z; k++) {
//					for (int j = 0; j < nxyz.y; j++) {
//						for (int i = 0; i < ntotGnAng; i++) {
//							angfluxLow[i + (k * dxyz[Oct].z* strideylow + j * dxyz[Oct].y) * ntotGnAng] = angfluxUp[i];
//						}
//					}
//				}
//				for (int i = 0; i < ntotGnAng; i++) angfluxUp[i] = zerophi;
//
//				if (isY) {
//					angfluxUp = &fluxIny[thisLv](0, 0, ixup + izup * stridexup);
//					angfluxLow = &fluxIny[nextLv](0, 0, ixlow + izlow * stridexlow);
//
//					for (int k = 0; k < nxyz.z; k++) {
//						for (int j = 0; j < nxyz.x; j++) {
//							for (int i = 0; i < ntotGnAng; i++) {
//								angfluxLow[i + (k * dxyz[Oct].z* stridexlow + j * dxyz[Oct].x) * ntotGnAng] = angfluxUp[i];
//							}
//						}
//					}
//					for (int i = 0; i < ntotGnAng; i++) angfluxUp[i] = zerophi;
//				}
//				if (isZ) {
//					angfluxUp = &fluxInz[thisLv](0, 0, ixup + iyup * stridexup);
//					angfluxLow = &fluxInz[nextLv](0, 0, ixlow + iylow * stridexlow);
//
//					for (int k = 0; k < nxyz.y; k++) {
//						for (int j = 0; j < nxyz.x; j++) {
//							for (int i = 0; i < ntotGnAng; i++) {
//								angfluxLow[i + (k * dxyz[Oct].y * stridexlow + j * dxyz[Oct].x) * ntotGnAng] = angfluxUp[i];
//							}
//						}
//					}
//					for (int i = 0; i < ntotGnAng; i++) angfluxUp[i] = zerophi;
//				}
//
//				thisLv = nextLv;
//			}
//
//			for (int lv = 1; lv < divlevel; lv++) {
//				upids[lv].x = idstore[lv].x / 3;
//				upids[lv].y = idstore[lv].y / 3;
//				upids[lv].z = idstore[lv].z / 3;
//			}
//			for (int lv = 1; lv < divlevel; lv++) {
//				idinbox[lv].x = idstore[lv].x % 3;
//				idinbox[lv].y = idstore[lv].y % 3;
//				idinbox[lv].z = idstore[lv].z % 3;
//			}
//			idinbox[0] = idstore[0];
//
//			
//		std::copy(fluxInx[0].begin(), fluxInx[0].end(), &(*xbndflux)(0, 0, 0, Oct));
//		if (isY) std::copy(fluxIny[0].begin(), fluxIny[0].end(), &(*ybndflux)(0, 0, 0, Oct));
//		if (isZ) std::copy(fluxInz[0].begin(), fluxInz[0].end(), &(*zbndflux)(0, 0, 0, Oct));
//		}
//	}
//}
}