#include "Defines.h"
#include "NuDEAL.h"

int main(int argc, char *argv[])
{
	NuDEAL::Master_t Master;

	Master.Initialize(argc, argv);



	Master.Finalize();

	return EXIT_SUCCESS;
}