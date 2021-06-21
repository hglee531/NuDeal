#pragma once
#include "Defines.h"
#include "Library.h"

namespace NuDEAL
{

class Master_t final
{
private:

	Library::XSLibrary_t *XS = nullptr;

public:

	void Initialize(int argc, char *argv[]);


	void Finalize();

};

}