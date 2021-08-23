#pragma once
#include "Defines.h"
#include "SolutionDefines.h"
#include "Array.h"
#include "PhysicalDomain.h"

constexpr double r9 = 1 / 9.0;

namespace Transport {
enum class AngQuadType{
	UNIFORM, // Uniform distribution on (theta, phi) space
	LS, // Level symmetric
	GC, // Gauss-Chebyshev
	Bi3, // Bickley-3 on polar, Uniform on azimuthal
	G_Bi3, // Bickley-3 on polar, Gaussian quadrature on aximuthal
	QR // Quadruple Range
};

inline int CalNangleOct(AngQuadType quadtype, vector<int> parameters) {
	int nangle_oct;
	switch (quadtype)
	{
	case Transport::AngQuadType::UNIFORM:
	{
		int nazi = parameters[0], npolar = parameters[1];
		nangle_oct = nazi * npolar;
		break;
	}
	case Transport::AngQuadType::LS:
		break;
	case Transport::AngQuadType::GC:
		break;
	case Transport::AngQuadType::Bi3:
	{
		int nazi = parameters[0], npolar = parameters[1];
		nangle_oct = nazi * npolar;
		break;
	}
	case Transport::AngQuadType::G_Bi3:
		break;
	case Transport::AngQuadType::QR:
		break;
	default:
		break;
	}

	return nangle_oct;
}

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
	vector<double> PtsXval;
public:
	AsymptoticExp(double tol, double xL, double xR);

	double ExpFast(double xval);

	double ExpSafe(double xval);
};

inline double Avg_Step(double flux0, double flux1, double tau, double srcflux) {
	// Flat source assumption (step characteristics)
	// Effective at high absorbing regions
	return srcflux - (flux1 - flux0) / tau;
}

inline double Avg_DD(double flux0, double flux1) {
	return 0.5 * (flux0 + flux1);
}

class DriverSN{
	template <typename T> using Array = LinPack::Array_t<T>;
	using RaySegment = PhysicalDomain::RayTracingDomain;
	using Boundary = PhysicalDomain::FluxBoundary;
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
	int ntotGnAng;
	array<int, 8> edgemod;
	array<int3, 8> edgeidx0, edgeidx;
	array<int3, 8> dxyz;

	int3 Nxyz, nxyz;
	double3 lxyz0;
	int nnode;
	int divlevel;
	const NodeInfo_t* nodeinfo;
	const int* innodeLv;

	int nloop;
	const Array<int> *sweepnodes;
	const Array<int2> *sweeploop;

	Array<double> trackL; // trackL[nangle_oct]
	double3 wtIF;
	vector<double3> wtVol;
	Array<double> wtAband, wtVband; // (3,3,3) per angle. for v,u = x,y,z, (v,u-rectangle), (v,u-parallelogram), (v,u-triangle)
	Array<double> trackLband;

	Array<double> wtbandP[2][2], wtbandT[2][2]; // (exp order, nangle)
	double3 areaIF, LsqrIF;
	vector<double> VBand;
	vector<double3> rnLxyz, Lpnxyz;
	double3 rnsqrxyz;
	
	int bandrank;
	double3 rnxyz;
	bool isY, isZ;

	Array<double> *scalarflux;
	const Array<double> *src;
	const Array<double> *xst;
	Array<realphi> *xbndflux, *ybndflux, *zbndflux;

	int nFXR, nFSR;
	vector<int> idFXR;
	const int *idFSR;
	const double *wFSR;

private:
	inline void AvgInFlux(int iangle, double *outx, double *outy, double *outz, vector<double> &avgfluxIF,  vector<double> &avgfluxVol) {
		for (int ig = 0; ig < ng; ig++) {
			avgfluxIF[ig] = wtIF.x * outx[ig]; avgfluxVol[ig] = wtVol[iangle].x * outx[ig];
			if (isY) { avgfluxIF[ig] += wtIF.y * outy[ig]; avgfluxVol[ig] += wtVol[iangle].y * outy[ig]; }
			if (isZ) { avgfluxIF[ig] += wtIF.z * outz[ig]; avgfluxVol[ig] += wtVol[iangle].z * outz[ig]; }
		}
	}

	inline void TrackAnode(int iangle, int thisLv, int iFXR, int iFSR, double wtFSR, double *angflux, vector<double> &avgfluxIF, vector<double> &avgfluxVol, const double *src) {
		for (int ig = 0; ig < ng; ig++) {
			double sigT = (*xst)(ig, iFXR);
			double angflux0[2] = { avgfluxIF[ig], avgfluxVol[ig] }, srcflux = src[ig] / sigT;

			double tau[2] = { 
				trackL(0, iangle, thisLv) * sigT,
				trackL(1, iangle, thisLv) * sigT };
			//double rates[2] = {
			//	exp(-tau[0]),
			//	(1. - exp(-tau[1])) / tau[1],
			//};
			double rates[2] = {
				Exponent.ExpFast(-tau[0]),
				(1. - Exponent.ExpFast(-tau[1])) / tau[1],
			};

			angflux[ig] = angflux0[0] * rates[0] + srcflux * (1. - rates[0]);
			(*scalarflux)(ig, iFSR) += (angflux0[1] * rates[1] + srcflux * (1. - rates[1])) * wtFSR * QuadSet.weights[iangle];

			//angflux[ig] = srcflux + (angflux0 - srcflux) * exp(-tau);// Exponent.ExpFast(-tau);
			//(*scalarflux)(ig, iFSR) += Avg_Step(angflux0, angflux[ig], tau, srcflux) * wtFSR * QuadSet.weights[iangle];

			/*if (scatorder >= 1) {
			scalarflux[ig + (iFXR + nFXR) * ng] +=
			scalarflux[ig + (iFXR + 2 * nFXR) * ng] +=
			scalarflux[ig + (iFXR + 3 * nFXR) * ng] +=
			} */
		}
	}

	void TrackBands(int iangle, int thisLv, int iFXR, int iFSR, double wtFSR, double *xflux, double *yflux, double *zflux, const double *src) {
		constexpr double r6 = 1. / 6.;
		double Len = trackLband(iangle, thisLv);
		for (int ig = 0; ig < ng; ig++) {
			double sigT = (*xst)(ig, iFXR);
			double tau = Len * sigT;
			//double srcflux = src[ig] / sigT;
			double srcflux = src[ig] * Len;
			double flux0[3] = { xflux[ig], (isY ? yflux[ig] : 0.0), (isZ ? zflux[ig] : 0.0) };

			double exptau = exp(-tau);
			//double exptau = Exponent.ExpFast(-tau);

			double reducedexp[5];
			reducedexp[0] = exptau;
			reducedexp[1] = (1. - exptau) / tau;
			reducedexp[2] = (1. - reducedexp[1]) / tau;
			reducedexp[3] = (0.5 - reducedexp[2]) / tau;
			reducedexp[4] = (r6 - reducedexp[3]) / tau;

			double flux1[3] = { 0.0, 0.0, 0.0 };

			for (int outdir = 0; outdir < bandrank; outdir++) {
				for (int indir = 0; indir < bandrank; indir++) {
					double totalloss = 0., totalacc = 0.;
					for (int band = 0; band < bandrank; band++) {
						totalloss += reducedexp[band] * wtAband(band, outdir, indir, iangle);
						totalacc  += reducedexp[band+1]  * wtVband(band, outdir, indir, iangle);
					}
					flux1[outdir] += flux0[indir] * totalloss;
					(*scalarflux)(ig, iFSR) += flux0[indir] * totalacc * wtFSR * QuadSet.weights[iangle];

					totalloss = totalacc = 0.;
					for (int band = 0; band < bandrank; band++) {
						totalloss += reducedexp[band+1] * wtAband(band, outdir, indir, iangle);
						totalacc  += reducedexp[band+2]  * wtVband(band, outdir, indir, iangle);
					}
					flux1[outdir] += srcflux * totalloss;
					(*scalarflux)(ig, iFSR) += srcflux * totalacc * wtFSR * QuadSet.weights[iangle];
				}
				//(*scalarflux)(ig, iFSR) += r6*(flux1[outdir] + flux0[outdir])*wtFSR*QuadSet.weights[iangle];
			}

			//double meanflux = flux1[0] * wtIF.x + flux1[1] * wtIF.y + flux1[2] * wtIF.y;
			//xflux[ig] = meanflux;
			//if (isY) yflux[ig] = meanflux;
			//if (isZ) zflux[ig] = meanflux;

			xflux[ig] = flux1[0]; 
			if (isY) yflux[ig] = flux1[1];
			if (isZ) zflux[ig] = flux1[2];
		}
	}

	void TrackBands(int iangle, int thisLv, int iFXR, int iFSR, double wtFSR, double2 *xflux, double2 *yflux, double2 *zflux, const double *src) {
		constexpr double r6 = 1. / 6., r24 = 1. / 24;
		double Len = trackLband(iangle, thisLv);
		double *areaIF_arr = &areaIF.x;
		for (int ig = 0; ig < ng; ig++) {
			double sigT = (*xst)(ig, iFXR);
			double tau = Len * sigT;
			//double srcflux = src[ig] / sigT;
			double srcflux = src[ig] * Len;

			double exptau = exp(-tau);
			//double exptau = Exponent.ExpFast(-tau);

			double reducedexp[6];
			reducedexp[0] = exptau;
			reducedexp[1] = (1. - exptau) / tau;
			reducedexp[2] = (1. - reducedexp[1]) / tau;
			reducedexp[3] = (0.5 - reducedexp[2]) / tau;
			reducedexp[4] = (r6 - reducedexp[3]) / tau;
			reducedexp[5] = (r24 - reducedexp[4]) / tau;

			double flux0[3][2], flux1[3][2] = { {0.0}, };
			flux0[0][0] = xflux[ig].x;  flux0[0][1] = xflux[ig].y;
			if (isY) { flux0[1][0] = yflux[ig].x; flux0[1][1] = yflux[ig].y; }

			// Only for 2D problems
			double scalar = 0.;
			for (int m_in = 0; m_in < 2; m_in++) {
				for (int m_out = 0; m_out < 2; m_out++) {
					for (int axis = 0; axis < 2; axis++) {
						int inaxi = axis;
						double wtband = wtbandP[m_in][m_out](0, inaxi, iangle);
						flux1[axis][m_out] += flux0[inaxi][m_in] * reducedexp[0] * wtband;
						if (m_in == 0) flux1[axis][m_out] += srcflux * reducedexp[1] * wtband;
						if (m_out == 0) {
							scalar += flux0[inaxi][m_in] * reducedexp[1] * wtband;
							if (m_in == 0) scalar += srcflux * reducedexp[2] * wtband;
						}

						inaxi = (axis + 1) % 2;
						for (int ordexp = 0; ordexp < 1 + m_in + m_out; ordexp++) {
							wtband = wtbandT[m_in][m_out](ordexp, inaxi, iangle);
							flux1[axis][m_out] += flux0[inaxi][m_in] * reducedexp[1 + ordexp] * wtband;
							if (m_in == 0) flux1[axis][m_out] += srcflux * reducedexp[2 + ordexp] * wtband;
							if (m_out == 0) {
								scalar += flux0[inaxi][m_in] * reducedexp[2 + ordexp] * wtband;
								if (m_in == 0) scalar += srcflux * reducedexp[3 + ordexp] * wtband;
							}
						}
					}
				}
			}

			for (int axis = 0; axis < 2; axis++) {
				for (int m = 0; m < 2; m++) {
					flux1[axis][m] /= areaIF_arr[axis];
				}
			}
			scalar /= VBand[iangle];

			xflux[ig] = make_double2(flux1[0][0], flux1[0][1]);
			if (isY) yflux[ig] = make_double2(flux1[1][0], flux1[1][1]);

			(*scalarflux)(ig, iFSR) += scalar * wtFSR * QuadSet.weights[iangle];
		}
	}

	inline void UpdateFluxIn(double *angflux, double *fluxiny, double *fluxinz) {
		if (isY) std::copy(angflux, angflux + ntotGnAng, fluxiny);
		if (isZ) std::copy(angflux, angflux + ntotGnAng, fluxinz);
	}

	inline void UpperOutFlux_X(int thisLv, int Oct, int3 ids, double *angfluxUp, double *angfluxLow) {
		for (int i = 0; i < ntotGnAng; i++) {
			angfluxUp[i] += rnxyz.x * angfluxLow[i];
		}
	}
	
	inline void UpperOutFlux_Y(int thisLv, int Oct, int3 ids, double *angfluxUp, double *angfluxLow) {
		for (int j = 0; j < nxyz.x; j++) {
			for (int i = 0; i < ntotGnAng; i++) {
				angfluxUp[i] += rnxyz.y * angfluxLow[i + j * ntotGnAng];
			}
		}
	}

	inline void UpperOutFlux_Z(int thisLv, int Oct, int3 ids, double *angfluxUp, double *angfluxLow) {
		for (int j = 0; j < nxyz.x * nxyz.y; j++) {
			for (int i = 0; i < ntotGnAng; i++) {
				angfluxUp[i] += rnxyz.z * angfluxLow[i + j * ntotGnAng];
			}
		}
	}

	inline void SaveUpperFlux(double *angfluxLow, double *angfluxUp) {
		std::copy(angfluxUp, angfluxUp + ntotGnAng, angfluxLow);
		std::fill(angfluxUp, angfluxUp + ntotGnAng, _ZERO);
	}

	inline void LowerOutFlux_X(int thisLv, int Oct, int3 ids, double *angfluxLow, double *angfluxUp) {
		std::copy(angfluxUp, angfluxUp + ntotGnAng, angfluxLow);
	}

	inline void LowerOutFlux_Y(int thisLv, int Oct, int3 ids, double *angfluxLow, double *angfluxUp) {
		for (int i = 0; i < nxyz.x; i++)
			std::copy(angfluxUp, angfluxUp + ntotGnAng, angfluxLow + i*ntotGnAng);
	}

	inline void LowerOutFlux_Z(int thisLv, int Oct, int3 ids, double *angfluxLow, double *angfluxUp) {
		for (int i = 0; i < nxyz.x*nxyz.y; i++)
			std::copy(angfluxUp, angfluxUp + ntotGnAng, angfluxLow + i * ntotGnAng);
	}

	inline void UpperOutFlux_X(int thisLv, int Oct, int3 ids, double2 *angfluxUp, double2 *angfluxLow) {
		double2 wt1;
		wt1.x = static_cast<double>(2 * ids.y - 1 - nxyz.y)*0.5*rnLxyz[thisLv - 1].y;
		wt1.y = rnsqrxyz.y;
		for (int i = 0; i < ntotGnAng; i++) {
			angfluxUp[i].x += rnxyz.x * angfluxLow[i].x;
			angfluxUp[i].y += rnxyz.x * (angfluxLow[i].x * wt1.x + angfluxLow[i].y * wt1.y);
		}
	}

	inline void UpperOutFlux_Y(int thisLv, int Oct, int3 ids, double2 *angfluxUp, double2 *angfluxLow) {
		double2 wt1;
		wt1.x = static_cast<double>(2 * edgeidx[Oct].x + 1 - nxyz.x)*0.5*rnLxyz[thisLv - 1].x;
		//wt1.x = static_cast<double>(1 - nxyz.x)*0.5*rnLxyz[thisLv - 1].x;
		wt1.y = rnsqrxyz.x;
		for (int j = 0; j < nxyz.x; j++) {
			for (int i = 0; i < ntotGnAng; i++) {
				angfluxUp[i].x += rnxyz.y * angfluxLow[i + j * ntotGnAng].x;
				angfluxUp[i].y += rnxyz.y * (angfluxLow[i + j * ntotGnAng].x * wt1.x + angfluxLow[i + j * ntotGnAng].y * wt1.y);
			}
			wt1.x += static_cast<double>(dxyz[Oct].x) * rnLxyz[thisLv - 1].x;
			//wt1.x += rnLxyz[thisLv - 1].x;
		}
	}

	inline void UpperOutFlux_Z(int thisLv, int Oct, int3 ids, double2 *angfluxUp, double2 *angfluxLow) {
		return;
	}

	inline void SaveUpperFlux(double2 *angfluxLow, double2 *angfluxUp) {
		std::copy(angfluxUp, angfluxUp + ntotGnAng, angfluxLow);
		std::fill(angfluxUp, angfluxUp + ntotGnAng, _ZERO2);
	}

	inline void LowerOutFlux_X(int thisLv, int Oct, int3 ids, double2 *angfluxLow, double2 *angfluxUp) {
		std::copy(angfluxUp, angfluxUp + ntotGnAng, angfluxLow);
		double wt1 = 0.5*Lpnxyz[thisLv].y*static_cast<double>(2 * ids.y + 1 - nxyz.y);
		for (int i = 0; i < ntotGnAng; i++) {
			angfluxLow[i].x += angfluxUp[i].y*wt1;
		}
	}

	inline void LowerOutFlux_Y(int thisLv, int Oct, int3 ids, double2 *angfluxLow, double2 *angfluxUp) {
		for (int i = 0; i < nxyz.x; i++)
			std::copy(angfluxUp, angfluxUp + ntotGnAng, angfluxLow + i * ntotGnAng);
		double wt1 = 0.5*Lpnxyz[thisLv].x*static_cast<double>(2 * edgeidx[Oct].x + 1 - nxyz.x);
		//double wt1 = 0.5*Lpnxyz[thisLv].x*static_cast<double>(1 - nxyz.x);
		for (int j = 0; j < nxyz.x; j++) {
			for (int i = 0; i < ntotGnAng; i++) {
				angfluxLow[i + j * ntotGnAng].x += angfluxUp[i].y*wt1;
			}
			wt1 += static_cast<double>(dxyz[Oct].x) * Lpnxyz[thisLv].x;
			//wt1 += Lpnxyz[thisLv].x;
		}
	}

	inline void LowerOutFlux_Z(int thisLv, int Oct, int3 ids, double2 *angfluxLow, double2 *angfluxUp) {
		return;
	}

public:
	void Initialize(RaySegment &Rays, Boundary &BndFlux, FlatSrc& FSR, const FlatXS& FXR);

	DriverSN(int ng, AngQuadType quadtype, vector<int> quadparameter);

	void RaySweep();

	void RaySweepRigorous();

	void RaySweepLinearFlux();
};
}