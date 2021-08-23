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

	//Geometry::DebugGeomHandle();
	SolutionDriver::DebugSolutionPin();
	//SolutionDriver::DebugSolutionBox();

	Master.Finalize();

	return EXIT_SUCCESS;
}