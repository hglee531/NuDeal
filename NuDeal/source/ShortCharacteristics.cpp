#include "ShortCharacteristics.h"
#include "Math.hpp"

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
		break;
	case Transport::AngQuadType::GC:
		break;
	case Transport::AngQuadType::Bi3:
		break;
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

void DriverSCMOC::Initialize(RaySegment& Rays, Boundary &BndFlux, FlatSrc& FSR, const FlatXS& FXR) {
	constexpr double r3  = 0.33333333333333333333333333;
	constexpr double r6  = 0.16666666666666666666666667;
	constexpr double r12 = 0.08333333333333333333333333;
	mode = Rays.GetDimension();
	Rays.GetBaseSizes(Nxyz, nxyz, nnode, divlevel);
	lxyz0 = Rays.GetNodeSizes();
	trackL.Create(QuadSet.nangle_oct); //wtIF.resize(3);

	int nx = nxyz.x, ny = nxyz.y, nz = nxyz.z;
	int onbot = nx * ny * (nz - 1);
	edgemod[0] = 0; edgemod[1] = nx; edgemod[2] = nx * (ny - 1); edgemod[3] = nx * ny - 1;
	edgemod[4] = onbot; edgemod[5] = onbot + nx; edgemod[6] = onbot + nx * (ny - 1); edgemod[7] = nx * ny*nz - 1;

	nodeinfo = Rays.GetBaseNodeInfo().data(); innodeLv = Rays.GetInNodeLv().data();
	sweepnodes = &Rays.GetSweepNodes(); sweeploop = &Rays.GetSweepLoop();
	scatorder = FXR.GetScatOrder();

	nloop = (sweeploop->size() / 8) - 1;
	nFXR = FXR.GetNblocks(); nFSR = FSR.GetNblocks();

	vector<double> rpow3, rpow9;

	double val = 1;
	for (int i = 0; i < divlevel; i++) { rpow3.push_back(val); rpow9.push_back(val*val); val /= 3.; }

	double Ax = lxyz0.y * lxyz0.z, Ay = lxyz0.x * lxyz0.z, Az = lxyz0.x * lxyz0.y;
	double Atot = Ax + Ay + Az;

	//wtIF.resize(3); wtIF[0] = Ax / Atot; wtIF[1] = Ay / Atot; wtIF[2] = Az / Atot;
	wtIF.x = Ax / Atot; wtIF.y = Ay / Atot; wtIF.z = Az / Atot;

	trackL.Create(2, QuadSet.nangle_oct, divlevel);
	wtAband.Create(3, 3, 3, QuadSet.nangle_oct); std::fill(wtAband.begin(), wtAband.end(), 0.0);
	wtVband.Create(3, 3, 3, QuadSet.nangle_oct); std::fill(wtVband.begin(), wtVband.end(), 0.0);
	trackLband.Create(QuadSet.nangle_oct, divlevel);

	for (int iang = 0; iang < QuadSet.nangle_oct; iang++) {
		double3 omega = QuadSet.omega[iang];
		double3 Ldeproj; // de-projected length
		Ldeproj.x = lxyz0.x / omega.x; Ldeproj.y = lxyz0.y / omega.y; Ldeproj.z = lxyz0.z / omega.z;

		double Lmin = min({ Ldeproj.x, Ldeproj.y, Ldeproj.z });
		double3 Lproj, Lrem;
		Lproj.x = Lmin * omega.x; Lproj.y = Lmin * omega.y; Lproj.z = Lmin * omega.z;
		Lrem.x = lxyz0.x - Lproj.x; Lrem.y = lxyz0.y - Lproj.y; Lrem.z = lxyz0.z - Lproj.z;

		// trackL generate
		double A1 = Lrem.x*Lrem.y + Lrem.y*Lrem.z + Lrem.x*Lrem.z;
		double A2 = Lproj.x * (Lrem.y + Lrem.z) + Lproj.y * (Lrem.x + Lrem.z) + Lproj.z * (Lrem.x + Lrem.y);
		double A3 = 2.0 * (Lproj.y * Lproj.z + Lproj.x * Lproj.y + Lproj.x * Lproj.z);
		double Atot_VpL = A1 + A2 * 0.5 + A3 * r6;
		//trackL(iang, 0) = Lmin / Atot * (1 + A2 / 2.0 + A3 / 3.0);
		trackL(0, iang, 0) = Lmin * (A1 + A2 * 0.5 + A3 * r6) / Atot; // Length for flux attenuation
		trackL(1, iang, 0) = Lmin * (A1 + A2 * r3 + A3 * r12) / Atot; // Length for flux average

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
		double V = Arec.x + Arec.y + Arec.z;
		V += (Apar[0].y+Apar[0].z + Apar[1].x+Apar[1].z + Apar[2].x+Apar[2].y) * 0.5;
		V += (Atri[0].y+Atri[0].z + Atri[1].x+Atri[1].z + Atri[2].x+Atri[2].y) * r6;
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
		
		for (int lv = 1; lv < divlevel; lv++) {
			trackL(0, iang, lv) = rpow3[lv] * trackL(0, iang, 0);
			trackL(1, iang, lv) = rpow3[lv] * trackL(1, iang, 0);
			trackLband(iang, lv) = rpow3[lv] * trackLband(iang, 0);
		}

		double tau = Lmin;
		double exptau = exp(-tau);

		double lossrate[3], srcloss[3];
		lossrate[0] = exptau; lossrate[1] = (1. - lossrate[0]) / tau; lossrate[2] = (1. - lossrate[1]) / tau;
		srcloss[0] = 1. - exptau; srcloss[1] = 1. - srcloss[0] / tau; srcloss[2] = 0.5 - srcloss[1] / tau;

		double accrate[3], srcacc[3];
		accrate[0] = (1. - lossrate[0]) / tau; accrate[1] = (1. - accrate[0]) / tau; accrate[2] = (0.5 - accrate[1]) / tau;
		srcacc[0] = 1. - srcloss[0] / tau; srcacc[1] = 0.5 - srcacc[0] / tau; srcacc[2] = r6 - srcacc[1] / tau;

		double outflux[2] = { 0. }, avgflux[2] = { 0. };
		for (int outdir = 0; outdir < 3; outdir++) {
			for (int indir = 0; indir < 3; indir++) {
				double totalloss = 0., totalacc = 0.;
				for (int band = 0; band < 3; band++) {
					totalloss += lossrate[band] * wtAband(band, outdir, indir, iang);
					totalacc += accrate[band] * wtVband(band, outdir, indir, iang);
				}
				outflux[0] += totalloss;
				avgflux[0] += totalacc;

				totalloss = totalacc = 0.;
				for (int band = 0; band < 3; band++) {
					totalloss += srcloss[band] * wtAband(band, outdir, indir, iang);
					totalacc += srcacc[band] * wtVband(band, outdir, indir, iang);
				}
				outflux[1] += totalloss;
				avgflux[1] += totalacc;
			}
		}

		double asmflux[2], asmavg[2];
		double asyL[2] = { trackL(0,iang,0), trackL(1,iang,0) };
		double rates[2] = {	exp(-asyL[0]), (1. - exp(-asyL[1])) / asyL[1] };
		asmflux[0] = rates[0]; asmflux[1] = 1. - rates[0];
		asmavg[0] = rates[1]; asmavg[1] = 1. - rates[1];

		asyL[0] = 0.;
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

DriverSCMOC::DriverSCMOC(int ng, AngQuadType quadtype, vector<int> quadparameter) 
	: Exponent(0.0001, -40.0, 0.0) 
{
	this->ng = ng;
	QuadSet.CreateSet(quadtype, quadparameter);
	ntotGnAng = ng * QuadSet.nangle_oct;
}

void DriverSCMOC::RaySweep() {
	for (int Oct = 0; Oct < 8; Oct++) {
		int inode;
		int thisLv = 0, thisId, iFSR, iFXR;
		double wtFSR;

		// Workspace for saving outgoing angular fluxes
		Array<double> fluxInx, fluxIny, fluxInz;

		fluxInx.Create(ng, QuadSet.nangle_oct, divlevel);
		fluxIny.Create(ng, QuadSet.nangle_oct, nxyz.x, divlevel);
		fluxInz.Create(ng, QuadSet.nangle_oct, nxyz.x*nxyz.y, divlevel);
		
		vector<int3> idstore(divlevel);

		for (int il = 0; il < nloop; il++) {
			int traceL = (*sweeploop)(il, Oct).x;
			int traceR = (*sweeploop)(il + 1, Oct).x, storedir = (*sweeploop)(il+1, Oct).y;

			double *flux_thisx;
			double *flux_thisy;
			double *flux_thisz;

			if (thisLv == 0) {
				flux_thisx = &(*xbndflux)(0, 0, idstore[thisLv].y + idstore[thisLv].z * Nxyz.y, Oct);
				flux_thisy = &(*ybndflux)(0, 0, idstore[thisLv].x + idstore[thisLv].z * Nxyz.x, Oct);
				flux_thisz = &(*zbndflux)(0, 0, idstore[thisLv].x + idstore[thisLv].y * Nxyz.x, Oct);
			}
			else {
				flux_thisx = &fluxInx(0, 0, thisLv);
				flux_thisy = &fluxIny(0, 0, idstore[thisLv].x, thisLv);
				flux_thisz = &fluxInz(0, 0, idstore[thisLv].x + idstore[thisLv].y * 3, thisLv);
			}

			double *angflux = flux_thisx;

			for (int itrace = traceL; itrace < traceR; itrace++) {
				int inode = (*sweepnodes)(itrace, Oct);
				int iFXR = idFXR[inode], iFSR = idFSR[inode];
				double wtFSR = wFSR[inode];
				
				for (int iangle = 0; iangle < QuadSet.nangle_oct; iangle++) {
					int stride = iangle * ng;
//#define NODALMEAN
#ifdef NODALMEAN
					AvgInFlux(flux_thisx + stride, flux_thisy + stride, flux_thisz + stride, angflux + stride);
					TrackAnode(iangle, thisLv, iFXR, iFSR, wtFSR, angflux + stride, &(*src)(0, iangle, iFSR, Oct));
#else
					TrackBands(iangle, thisLv, iFXR, iFSR, wtFSR, 
						flux_thisx + stride, flux_thisy + stride, flux_thisz + stride, &(*src)(0, iangle, iFSR, Oct));
#endif
				}

#ifdef NODALMEAN
				//UpdateFluxIn(angflux, flux_thisy, flux_thisz);
#endif

				idstore[thisLv].x++;

				flux_thisy += ntotGnAng;
				flux_thisz += ntotGnAng;
			}

			// Update outgoing fluxes of the upper level
			double *angfluxUp, *angfluxLow;
			if (storedir > 0) {
				if (thisLv > 0) {
					if (thisLv == 1) angfluxUp = &(*xbndflux)(0, 0, idstore[0].y + idstore[0].z * Nxyz.y, Oct);
					else angfluxUp = &fluxInx(0, 0, thisLv - 1);

					angfluxLow = &fluxInx(0, 0, thisLv);

					UpperOutFlux_X(angfluxUp, angfluxLow);
				}
				idstore[thisLv].x = 0;
				idstore[thisLv].y++;
			}
			if (storedir > 1) {
				if (thisLv > 0) {
					if (thisLv == 1) angfluxUp = &(*ybndflux)(0, 0, idstore[0].x + idstore[0].z * Nxyz.x, Oct);
					else angfluxUp = &fluxIny(0, 0, idstore[thisLv - 1].x, thisLv - 1);

					angfluxLow = &fluxIny(0, 0, 0, thisLv);

					UpperOutFlux_Y(angfluxUp, angfluxLow);
				}
				idstore[thisLv].y = 0;
				idstore[thisLv].z++;
			}
			if (storedir > 2) {
				if (thisLv > 0) {
					if (thisLv == 1) angfluxUp = &(*zbndflux)(0, 0, idstore[0].x + idstore[0].y * Nxyz.x, Oct);
					else angfluxUp = &fluxInz(0, 0, idstore[thisLv - 1].x + idstore[thisLv - 1].y * nxyz.x, thisLv - 1);

					angfluxLow = &fluxInz(0, 0, 0, thisLv);

					UpperOutFlux_Z(angfluxUp, angfluxLow);

					idstore[thisLv].z = 0;
					thisLv--;
					idstore[thisLv].x++;
				}
			}
			
			// Update incoming fluxes of the lower level
			if (storedir == 0) {
				int nextLv = thisLv + 1;
				int id0, id1;
				double *angfluxLowX, *angfluxLowY, *angfluxLowZ;
				double *angfluxUpX, *angfluxUpY, *angfluxUpZ;

				angfluxLowX = &fluxInx(0, 0, nextLv);
				angfluxLowY = &fluxIny(0, 0, 0, nextLv);
				angfluxLowZ = &fluxInz(0, 0, 0, nextLv);

				if (thisLv == 0) {
					angfluxUpX = &(*xbndflux)(0, 0, idstore[thisLv].y + idstore[thisLv].z * Nxyz.y, Oct);
					angfluxUpY = &(*ybndflux)(0, 0, idstore[thisLv].x + idstore[thisLv].z * Nxyz.x, Oct);
					angfluxUpZ = &(*zbndflux)(0, 0, idstore[thisLv].x + idstore[thisLv].y * Nxyz.x, Oct);
				}
				else {
					angfluxUpX = &fluxInx(0, 0, thisLv);
					angfluxUpY = &fluxIny(0, 0, idstore[thisLv].x, thisLv);
					angfluxUpZ = &fluxInz(0, 0, idstore[thisLv].x + idstore[thisLv].y * 3, thisLv);
				}

				LowerOutFlux_X(angfluxLowX, angfluxUpX);
				LowerOutFlux_Y(angfluxLowY, angfluxUpY);
				LowerOutFlux_Z(angfluxLowZ, angfluxUpZ);

				thisLv = nextLv;
			}
		}
	}
}
}