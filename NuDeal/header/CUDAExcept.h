#pragma once
#include "cuda.h"
#include "cublas_v2.h"
#include "cusparse.h"
#include "curand.h"

// CUDA Error Check Macros

#define cudaCheckError(error) \
if (error != cudaError::cudaSuccess) { \
	std::cout << cudaGetErrorString(error) << std::endl; \
	std::cout << "In file " << __FILE__ \
		<< " line " << __LINE__ \
		<< " func " << __func__ << std::endl; \
	::exit(EXIT_FAILURE); \
}

#define cudaCheckErrorRAND(error) \
if (error != curandStatus_t::CURAND_STATUS_SUCCESS) { \
	std::cout << "cuRAND error code : " << error << std::endl; \
	std::cout << "In file " << __FILE__ \
		<< " line " << __LINE__ \
		<< " func " << __func__ << std::endl; \
	::exit(EXIT_FAILURE); \
}

#define cudaCheckErrorBLAS(error) \
if (error != cublasStatus_t::CUBLAS_STATUS_SUCCESS) { \
	std::cout << "cuBLAS error code : " << error << std::endl; \
	std::cout << "In file " << __FILE__ \
		<< " line " << __LINE__ \
		<< " func " << __func__ << std::endl; \
	::exit(EXIT_FAILURE); \
}

#define cudaCheckErrorSPARSE(error) \
if (error != cusparseStatus_t::CUSPARSE_STATUS_SUCCESS) { \
	std::cout << "cuSPARSE error code : " << error << std::endl; \
	std::cout << "In file " << __FILE__ \
		<< " line " << __LINE__ \
		<< " func " << __func__ << std::endl; \
	::exit(EXIT_FAILURE); \
}