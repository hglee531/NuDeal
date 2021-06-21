#include "NuDEAL.h"
#include "MPI.hpp"
#include "XSLibrary.h"

namespace NuDEAL
{

void Master_t::Initialize(int argc, char *argv[])
{
	MPI::Init(&argc, &argv);

	XS = new Library::XSLibrary_t;


	XS->ReadMacro();
}


void Master_t::Finalize()
{


	MPI::Finalize();
}

}