#include "Defines.h"
#include "UnitGeo.h"
#include "GeoHandle.h"
#include "Input.h"
#include "HardCodeParam.h"
#include "XS.h"
#include "Array.hpp"
#include "PhysicsProp.h"

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);

	//LinPack::Array_t<double> a;
	//LinPack::Array_t<double> b;
	//LinPack::Array_t<int> c;

	//vector<double> doubleVec(100);

	//IO::InputManager_t Parser;
	//std::string file = std::string(argv[1]);
	//Parser.ReadInput(file);

	//Geometry::DebugUnitGeo();
	//Geometry::DebugGeomHandle();
	//PhysicalDomain::DebugPhysicalDomain();

	//XS::XSLib C5G7Data(false, ng, nXSset, 0);
	//C5G7Data.UploadXSData(typeXS, XSSet, XSSM);

	MPI_Finalize();

	return EXIT_SUCCESS;
}