#include "Defines.h"
#include "NuDEAL.h"

#include "PhysicalDomain.h"
#include "XSLibrary.h"
#include "ShortCharacteristics.h"
#include "SolutionDriver.h"

using namespace Geometry;
using namespace PhysicalDomain;
using namespace Transport;
using namespace SolutionDriver;

int main(int argc, char *argv[])
{
	NuDEAL::Master_t Master;

	Master.Initialize(argc, argv);

	double c_cir[3] = { 0.0, 0.0, 0.54 };
	UnitSurf Circle(UnitSurf::SurfType::CIRCLE, c_cir, CartPlane::XY);

	double c_zpln0[2] = { -1.0, 0.0 }, c_zpln1[2] = { 1.0, 1. };
	UnitSurf Zpln0(UnitSurf::SurfType::ZPLN, c_zpln0), Zpln1(UnitSurf::SurfType::ZPLN, c_zpln1);

	double c_xpln0[2] = { -1.0, 0.0 }, c_xpln1[2] = { 1.0, 1. };
	double c_ypln0[2] = { -1.0, 0.0 }, c_ypln1[2] = { 1.0, 1. };
	UnitSurf Xpln0(UnitSurf::SurfType::XPLN, c_xpln0), Xpln1(UnitSurf::SurfType::XPLN, c_xpln1);
	UnitSurf Ypln0(UnitSurf::SurfType::YPLN, c_ypln0), Ypln1(UnitSurf::SurfType::YPLN, c_ypln1);

	UnitVol Cylinder;
	Cylinder << Circle << Zpln0 << Zpln1;
	//Cylinder.Finalize();

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

	UnitComp BoxComp(0, Box), CylComp(0, Cylinder);
	//CylComp.Finalize(); 
	BoxComp.Finalize();

	GeometryHandler GeoHandle;
	double origin[3] = { 0.0, 0.0, 0.0 }, L[3] = { 1., 1., 1. };
	GeoHandle.SetOrdinates(origin, L);
	GeoHandle << BoxComp;// << CylComp;
	GeoHandle.FinalizeVolumes();
	GeoHandle.Discretize(Dimension::TwoD, 1.1, 0.001, 1);

	RayTracingDomain RaySP(GeoHandle);
	FluxBoundary BndFlux(RaySP.GetNnodeZero().x, RaySP.GetNnodeZero().y, RaySP.GetNnodeZero().z);
	
	int3 blockFSR = make_int3(1, 1, 1), blockFXR = make_int3(1, 1, 1);
	FlatSrcDomain FSR(blockFSR, RaySP);
	FlatXSDomain FXR(blockFXR, FSR);

	Library::XSLibrary_t XSLib;
	XSLib.ReadMacro();

	vector<int> quadparam(2);
	quadparam[0] = 1; quadparam[1] = 1;
	int nangle_oct = CalNangleOct(AngQuadType::UNIFORM, quadparam);
	
	BndFlux.Initialize(XSLib.GetNumGroup(), nangle_oct);
	FSR.Initialize(XSLib.GetNumGroup(), nangle_oct);
	FXR.InitializeMacro(XSLib.GetNumGroup(), XSLib.GetScatOrder(), false);
	FXR.SetMacroXS(XSLib, GeoHandle.GetMatIds());

	DriverSCMOC MOCDriver(XSLib.GetNumGroup(), AngQuadType::UNIFORM, quadparam);
	MOCDriver.Initialize(RaySP, BndFlux, FSR, FXR);
	//MOCDriver.RaySweep();

	SteadyStateDriver SSDriver(Dimension::TwoD, XSLib.GetNumGroup(), nangle_oct, XSLib.GetScatOrder());
	SSDriver.Initialize(RaySP, BndFlux, FSR, FXR, MOCDriver);
	SSDriver.SetCriterionPsi(1.e-5);
	SSDriver.SetMaxIter(1000);

	SSDriver.RunSS();

	Master.Finalize();

	return EXIT_SUCCESS;
}