#pragma once
#include "Defines.h"
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
		break;
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

class DriverSCMOC{
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
	Array<double> wtAband, wtVband; // (3,3,3) per angle. for v,u = x,y,z, (v,u-rectangle), (v,u-parallelogram), (v,u-triangle)
	Array<double> trackLband;
	//vector<double> wtIF; // wtIF[3]
	//Array<double> wtIF; // weights from interfaces along x,y,z axes, wtIF[nangle_oct][3]

	Array<double> *scalarflux;
	const Array<double> *src;
	const Array<double> *xst;
	Array<double> *xbndflux, *ybndflux, *zbndflux;

	int nFXR, nFSR;
	vector<int> idFXR;
	const int *idFSR;
	const double *wFSR;

private:
	inline void AvgInFlux(double *outx, double *outy, double *outz, double *avgin) {
		for (int ig = 0; ig < ng; ig++) {
			avgin[ig] = wtIF.x * outx[ig] + wtIF.y * outy[ig] + wtIF.z * outz[ig];
		}
	}

	inline void TrackAnode(int iangle, int thisLv, int iFXR, int iFSR, double wtFSR, double *angflux, const double *src) {
		for (int ig = 0; ig < ng; ig++) {
			double sigT = (*xst)(ig, iFXR);
			double angflux0 = angflux[ig], srcflux = src[ig] / sigT;

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

			angflux[ig] = angflux0 * rates[0] + srcflux * (1. - rates[0]);
			(*scalarflux)(ig, iFSR) += (angflux0 * rates[1] + srcflux * (1. - rates[1])) * wtFSR * QuadSet.weights[iangle];

			//angflux[ig] = srcflux + (angflux0 - srcflux) * exp(-tau);// Exponent.ExpFast(-tau);
			//(*scalarflux)(ig, iFSR) += Avg_Step(angflux0, angflux[ig], tau, srcflux) * wtFSR * QuadSet.weights[iangle];

			/*if (scatorder >= 1) {
			scalarflux[ig + (iFXR + nFXR) * ng] +=
			scalarflux[ig + (iFXR + 2 * nFXR) * ng] +=
			scalarflux[ig + (iFXR + 3 * nFXR) * ng] +=
			} */
		}
	}

	inline void TrackBands(int iangle, int thisLv, int iFXR, int iFSR, double wtFSR, double *xflux, double *yflux, double *zflux, const double *src) {
		constexpr double r6 = 0.1666666666666666666666666666666667;
		for (int ig = 0; ig < ng; ig++) {
			double sigT = (*xst)(ig, iFXR);
			double tau = trackLband(iangle, thisLv) * sigT;
			double srcflux = src[ig] / sigT;
			double flux0[3] = { xflux[ig], yflux[ig], zflux[ig] };

			//double exptau = exp(-tau);
			double exptau = Exponent.ExpFast(-tau);
			exptau = exp(-tau);
			
			double lossrate[3], srcloss[3];
			lossrate[0] = exptau; lossrate[1] = (1. - lossrate[0]) / tau; lossrate[2] = (1. - lossrate[1]) / tau;
			srcloss[0] = 1. - exptau; srcloss[1] = 1. - srcloss[0] / tau; srcloss[2] = 0.5 - srcloss[1] / tau;

			double accrate[3], srcacc[3];
			accrate[0] = (1. - lossrate[0]) / tau; accrate[1] = (1. - accrate[0]) / tau; accrate[2] = (0.5 - accrate[1]) / tau;
			srcacc[0] = 1. - srcloss[0] / tau; srcacc[1] = 0.5 - srcacc[0] / tau; srcacc[2] = r6 - srcacc[1] / tau;

			double flux1[3] = { 0.0, 0.0, 0.0 };

			for (int outdir = 0; outdir < 3; outdir++) {
				for (int indir = 0; indir < 3; indir++) {
					double totalloss = 0., totalacc = 0.;
					for (int band = 0; band < 3; band++) {
						totalloss += lossrate[band] * wtAband(band, outdir, indir, iangle);
						totalacc  += accrate[band]  * wtVband(band, outdir, indir, iangle);
					}
					flux1[outdir] += flux0[indir] * totalloss;
					(*scalarflux)(ig, iFSR) += flux0[indir] * totalacc * wtFSR * QuadSet.weights[iangle];

					totalloss = totalacc = 0.;
					for (int band = 0; band < 3; band++) {
						totalloss += srcloss[band] * wtAband(band, outdir, indir, iangle);
						totalacc  += srcacc[band]  * wtVband(band, outdir, indir, iangle);
					}
					flux1[outdir] += srcflux * totalloss;
					(*scalarflux)(ig, iFSR) += srcflux * totalacc * wtFSR * QuadSet.weights[iangle];
				}
				//(*scalarflux)(ig, iFSR) += r6*(flux1[outdir] + flux0[outdir])*wtFSR*QuadSet.weights[iangle];
			}

			xflux[ig] = flux1[0]; yflux[ig] = flux1[1]; zflux[ig] = flux1[2];
		}
	}

	inline void UpdateFluxIn(double *angflux, double *fluxiny, double *fluxinz) {
		std::copy(angflux, angflux + ntotGnAng, fluxiny);
		std::copy(angflux, angflux + ntotGnAng, fluxinz);
	}

	inline void UpperOutFlux_X(double *angfluxUp, double *angfluxLow) {
		for (int i = 0; i < ntotGnAng; i++) {
			angfluxUp[i] += r9 * angfluxLow[i];
		}
	}
	
	inline void UpperOutFlux_Y(double *angfluxUp, double *angfluxLow) {
		for (int j = 0; j < nxyz.x; j++) {
			for (int i = 0; i < ntotGnAng; i++) {
				angfluxUp[i] += r9 * angfluxLow[i + j * ntotGnAng];
			}
		}
	}

	inline void UpperOutFlux_Z(double *angfluxUp, double *angfluxLow) {
		for (int j = 0; j < nxyz.x * nxyz.y; j++) {
			for (int i = 0; i < ntotGnAng; i++) {
				angfluxUp[i] += r9 * angfluxLow[i + j * ntotGnAng];
			}
		}
	}

	inline void LowerOutFlux_X(double *angfluxLow, double *angfluxUp) {
		std::copy(angfluxUp, angfluxUp + ntotGnAng, angfluxLow);
	}

	inline void LowerOutFlux_Y(double *angfluxLow, double *angfluxUp) {
		for (int i = 0; i < nxyz.x; i++) 
			std::copy(angfluxUp + i * ntotGnAng, angfluxUp + (i + 1)*ntotGnAng, angfluxLow);

	}

	inline void LowerOutFlux_Z(double *angfluxLow, double *angfluxUp) {
		for (int i = 0; i < nxyz.x*nxyz.y; i++) 
			std::copy(angfluxUp + i * ntotGnAng, angfluxUp + (i + 1)*ntotGnAng, angfluxLow);
	}

public:
	void Initialize(RaySegment &Rays, Boundary &BndFlux, FlatSrc& FSR, const FlatXS& FXR);

	DriverSCMOC(int ng, AngQuadType quadtype, vector<int> quadparameter);

	void RaySweep();
};
}