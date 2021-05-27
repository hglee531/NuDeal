#pragma once
#include "Array.h"
#include "CUDAExcept.h"

namespace LinPack
{

template <typename T>
void Array_t<T>::SetDimension(size_type nx, size_type ny, size_type nz, size_type nw)
{
	this->nx = nx;
	this->ny = ny;
	this->nz = nz;
	this->nw = nw;
	this->nxy = nx * ny;
	this->nxyz = nx * ny * nz;
	this->n = nx * ny * nz * nw;
}

template <typename T>
void Array_t<T>::Alias(const Array_t<T>& rhs)
{
	SetDimension(rhs.nx, rhs.ny, rhs.nz, rhs.nw);
	ptr_host = rhs.ptr_host;
	ptr_device = rhs.ptr_device;
}

template <typename T>
void Array_t<T>::AliasHost(const Array_t<T>& rhs)
{
	SetDimension(rhs.nx, rhs.ny, rhs.nz, rhs.nw);
	ptr_host = rhs.ptr_host;
}

template <typename T>
void Array_t<T>::AliasHost(pointer ptr, size_type nx, size_type ny, size_type nz, size_type nw)
{
	SetDimension(nx, ny, nz, nw);
	ptr_host = ptr;
}

template <typename T>
void Array_t<T>::AliasDevice(const Array_t<T>& rhs)
{
	SetDimension(rhs.nx, rhs.ny, rhs.nz, rhs.nw);
	ptr_device = rhs.ptr_device;
}

template <typename T>
void Array_t<T>::AliasDevice(pointer ptr, size_type nx, size_type ny, size_type nz, size_type nw)
{
	SetDimension(nx, ny, nz, nw);
	ptr_device = ptr;
}

template <typename T>
void Array_t<T>::ResizeHost(size_type nx, size_type ny, size_type nz, size_type nw)
{
	ClearHost();
	SetDimension(nx, ny, nz, nw);
	ptr_host = new T[n];
}

template <typename T>
void Array_t<T>::ResizeHost(const_pointer ptr, size_type nx, size_type ny, size_type nz, size_type nw)
{
	ResizeHost(nx, ny, nz, nw);
	std::copy(std::execution::par, ptr, ptr + n, begin());
}

template <typename T>
void Array_t<T>::ResizeDevice(size_type nx, size_type ny, size_type nz, size_type nw)
{
	ClearDevice();
	SetDimension(nx, ny, nz, nw);
	cudaCheckError( cudaMalloc(&ptr_device, n * sizeof(T)) ); 
}

template <typename T>
void Array_t<T>::ResizeDevice(const_pointer ptr, size_type nx, size_type ny, size_type nz, size_type nw)
{
	ResizeDevice(nx, ny, nz, nw);
	cudaCheckError( cudaMemcpy(ptr_device, ptr, n * sizeof(T), cudaMemcpyHostToDevice) );
}

template <typename T>
void Array_t<T>::ClearHost()
{
	if (IsHostAlloc()) delete[] ptr_host;
	ptr_host = static_cast<pointer>(NULL);
}

template <typename T>
void Array_t<T>::ClearDevice()
{
	if (IsDeviceAlloc()) cudaFree(ptr_device);
	ptr_device = static_cast<pointer>(NULL);
}

template <typename T>
inline Array_t<T>& Array_t<T>::operator=(const_reference val)
{
	std::fill(std::execution::par, begin(), end(), val);
	return *this;
}

template <typename T>
inline Array_t<T>& Array_t<T>::operator=(const_pointer ptr)
{
	std::copy(std::execution::par, ptr, ptr + size(), begin());
	return *this;
}

template <typename T> template <typename U>
inline Array_t<T>& Array_t<T>::operator=(const Array_t<U>& rhs)
{
	if (!IsHostAlloc() && !IsHostAlias()) ResizeHost(rhs.nx, rhs.ny, rhs.nz, rhs.nw);
	#pragma omp parallel for schedule(guided)
	for (index_type i = 0; i < size(); ++i) ptr_host[i] = static_cast<T>(rhs.ptr_host[i]);
	return *this;
}

template <typename T>
inline Array_t<T>& Array_t<T>::operator=(const Array_t<T>& rhs)
{
	if (!IsHostAlloc() && !IsHostAlias()) ResizeHost(rhs.nx, rhs.ny, rhs.nz, rhs.nw);
	std::copy(std::execution::par, rhs.begin(), rhs.end(), begin());
	return *this;
}

template <typename T>
inline Array_t<T>& Array_t<T>::operator=(Array_t<T>&& rhs)
{
	AliasHost(rhs);
	rhs.ptr_host = NULL;
	return *this;
}

template <typename T>
inline Array_t<T>& Array_t<T>::operator+=(const_reference val)
{
	#pragma omp parallel for schedule(guided)
	for (index_type i = 0; i < size(); ++i) ptr_host[i] += val;
	return *this;
}

template <typename T> template <typename U>
inline Array_t<T>& Array_t<T>::operator+=(const Array_t<U>& rhs)
{
	#pragma omp parallel for schedule(guided)
	for (index_type i = 0; i < size(); ++i) ptr_host[i] += rhs.ptr_host[i];
	return *this;
}

template <typename T>
inline Array_t<T>& Array_t<T>::operator+=(const Array_t<T>& rhs)
{
	#pragma omp parallel for schedule(guided)
	for (index_type i = 0; i < size(); ++i) ptr_host[i] += rhs.ptr_host[i];
	return *this;
}

template <typename T>
inline Array_t<T>& Array_t<T>::operator-=(const_reference val)
{
	#pragma omp parallel for schedule(guided)
	for (index_type i = 0; i < size(); ++i) ptr_host[i] -= val;
	return *this;
}

template <typename T> template <typename U>
inline Array_t<T>& Array_t<T>::operator-=(const Array_t<U>& rhs)
{
	#pragma omp parallel for schedule(guided)
	for (index_type i = 0; i < size(); ++i) ptr_host[i] -= rhs.ptr_host[i];
	return *this;
}

template <typename T>
inline Array_t<T>& Array_t<T>::operator-=(const Array_t<T>& rhs)
{
	#pragma omp parallel for schedule(guided)
	for (index_type i = 0; i < size(); ++i) ptr_host[i] -= rhs.ptr_host[i];
	return *this;
}

template <typename T>
inline Array_t<T>& Array_t<T>::operator*=(const_reference val)
{
	#pragma omp parallel for schedule(guided)
	for (index_type i = 0; i < size(); ++i) ptr_host[i] *= val;
	return *this;
}

template <typename T> template <typename U>
inline Array_t<T>& Array_t<T>::operator*=(const Array_t<U>& rhs)
{
	#pragma omp parallel for schedule(guided)
	for (index_type i = 0; i < size(); ++i) ptr_host[i] *= rhs.ptr_host[i];
	return *this;
}

template <typename T>
inline Array_t<T>& Array_t<T>::operator*=(const Array_t<T>& rhs)
{
	#pragma omp parallel for schedule(guided)
	for (index_type i = 0; i < size(); ++i) ptr_host[i] *= rhs.ptr_host[i];
	return *this;
}

template <typename T>
inline Array_t<T>& Array_t<T>::operator/=(const_reference val)
{
	#pragma omp parallel for schedule(guided)
	for (index_type i = 0; i < size(); ++i) ptr_host[i] /= val;
	return *this;
}

template <typename T> template <typename U>
inline Array_t<T>& Array_t<T>::operator/=(const Array_t<U>& rhs)
{
	#pragma omp parallel for schedule(guided)
	for (index_type i = 0; i < size(); ++i) ptr_host[i] /= rhs.ptr_host[i];
	return *this;
}

template <typename T>
inline Array_t<T>& Array_t<T>::operator/=(const Array_t<T>& rhs)
{
	#pragma omp parallel for schedule(guided)
	for (index_type i = 0; i < size(); ++i) ptr_host[i] /= rhs.ptr_host[i];
	return *this;
}

}