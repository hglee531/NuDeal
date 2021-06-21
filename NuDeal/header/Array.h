#pragma once
#include "Defines.h"
#include "CUDAExcept.h"
#include "CUDATypeTraits.h"

#define _HDINLINE_ __inline__ __host__ __device__

namespace LinPack
{

namespace _Detail // Namesapce Detail
{

template <typename T, unsigned N>
class _Static_array
{
public:
	T _Elems[N] = { 0, };

	_Static_array() = default;

	_HDINLINE_ T *begin() noexcept { return _Elems; }
	_HDINLINE_ const T* begin() const noexcept { return _Elems; }
	_HDINLINE_ T *end() noexcept { return _Elems + N; }
	_HDINLINE_ const T* end() const noexcept { return _Elems + N; }
	_HDINLINE_ T& front() noexcept { return *_Elems; }
	_HDINLINE_ const T& front() const noexcept { return *_Elems; }
	_HDINLINE_ T& back() noexcept { return *(_Elems + N - 1); }
	_HDINLINE_ const T& back() const noexcept { return *(_Elems + N - 1); }
	_HDINLINE_ T& operator[] (size_t idx) noexcept { return _Elems[idx]; }
	_HDINLINE_ const T& operator[] (size_t idx) const noexcept { return _Elems[idx]; }
};

template <unsigned N, typename... Ts>
using _call_if_t = typename std::enable_if_t<
	std::conjunction_v<std::is_integral<Ts>...> && sizeof...(Ts) <= N && N>;

template <typename T>
using _not_portable_t = typename std::enable_if_t<!CUDATypeTraits::is_portable_v<T>>;

template <typename T>
using _portable_t = typename std::enable_if_t<CUDATypeTraits::is_portable_v<T>>;

template <typename T>
struct _host_memory
{
	static constexpr T* malloc(size_t n) 
	{ 
		return n == 0 ? static_cast<T*>(nullptr) : new T[n]; 
	}
	static void free(T *ptr) 
	{ 
		if (!ptr) return;
		delete[] ptr; 
	}
	static void memset(T *begin, T *end, const T& val = T{})
	{
		if (!(end - begin)) return;
		std::fill(std::execution::par, begin, end, val);
	}
	static void memcpy(const T* first, const T* last, T* dest)
	{
		if (!(last - first)) return;
		std::copy(std::execution::par, first, last, dest);
	}
};

template <typename T>
struct _device_memory
{
	static constexpr T* malloc(size_t n) 
	{ 
		T *ptr = nullptr;
		if (n == 0) return ptr;
		cudaCheckError( cudaMalloc(&ptr, n * sizeof(T)) );
		return ptr;
	}
	static constexpr void free(T *ptr) 
	{ 
		if (!ptr) return;
		cudaCheckError( cudaFree(ptr) ); 
	}
	static constexpr void memset(T *begin, T *end)
	{
		if (!(end - begin)) return;
		cudaCheckError( cudaMemset(begin, 0x0, (end - begin) * sizeof(T)) );
	}
	static constexpr void memset(T *begin, T *end, const T& val)
	{
		size_t sz = end - begin;
		if (!sz) return;
		T *ptr = new T[sz]; std::fill(std::execution::par, ptr, ptr + sz, val);
		cudaCheckError( cudaMemcpy(begin, ptr, sz * sizeof(T), cudaMemcpyHostToDevice) );
		delete[] ptr;
	}
	static constexpr void memcpy(const T* first, const T* last, T* dest)
	{
		if (!(last - first)) return;
		cudaCheckError( 
			cudaMemcpy(dest, first, (last - first) * sizeof(T), cudaMemcpyDeviceToDevice) 
		);
	}
};

template <typename T> struct _memory_traits : false_type {};
template <typename T> struct _memory_traits<_host_memory<T>> : true_type{};
template <typename T> struct _memory_traits<_device_memory<T>> : true_type{};

template <typename T, typename _Memory,
  typename = typename std::enable_if_t<_memory_traits<_Memory>::value>>
class _Vector_t
{
public:

	using memory_type = _Memory;
	using value_type = T;
	using pointer = T*;
	using const_pointer = const T*;
	using reference = T&;
	using const_reference = const T&;
	using size_type = std::size_t;
	using difference_type = std::ptrdiff_t;
	using iterator = T*;
	using const_iterator = const T*;

protected:
	
	pointer _Myfirst = nullptr, _Mylast = nullptr;

public:

	_Vector_t() noexcept = default;
	explicit _Vector_t(const size_type count) :
		_Myfirst{ memory_type::malloc(count) },
		_Mylast{ _Myfirst + count }
	{
		memory_type::memset(_Myfirst, _Mylast);
	}
	explicit _Vector_t(const size_type count, const T& val) :
		_Vector_t{ count }
	{
		memory_type::memset(_Myfirst, _Mylast, val);
	}
	
public:

	_Vector_t(const _Vector_t& rhs) :
		_Myfirst{ memory_type::malloc(rhs.size()) },
		_Mylast{ _Myfirst + rhs.size() }
	{
		memory_type::memcpy(rhs._Myfirst, rhs._Mylast, _Myfirst);
	}
	_Vector_t(_Vector_t&& rhs) noexcept :
		_Myfirst{ std::move(rhs._Myfirst) },
		_Mylast{ std::move(rhs._Mylast) }
	{
		rhs._Myfirst = rhs._Mylast = nullptr;
	}

protected:
	void _assign_range(const_pointer first, const_pointer last)
	{
		size_type sz = last - first;
		size_type old_sz = _Mylast - _Myfirst;
		if (sz != old_sz) {
			memory_type::free(_Myfirst);
			_Myfirst = memory_type::malloc(sz);
			_Mylast = _Myfirst + sz;
		}
		memory_type::memcpy(first, last, _Myfirst);
	}

	void _swap_all(_Vector_t& rhs)
	{
		std::swap(_Myfirst, rhs._Myfirst);
		std::swap(_Mylast, rhs._Mylast);
	}

public:

	_Vector_t& operator=(const _Vector_t& rhs)
	{
		_assign_range(rhs._Myfirst, rhs._Mylast);
		return *this;
	}
	_Vector_t& operator=(_Vector_t&& rhs)
	{
		if (this != std::addressof(rhs)) {
			_swap_all(rhs);
			rhs.clear();
		}
		return *this;
	}

public:
	
	void clear() noexcept
	{
		memory_type::free(_Myfirst);
		_Myfirst = _Mylast = nullptr;
	}

public:

	~_Vector_t() { clear(); }

public:

	void resize(size_type count)
	{
		const size_type old_size = _Mylast - _Myfirst;
		if (count == old_size) return;
		pointer old_first = _Myfirst;
		_Myfirst = memory_type::malloc(count);
		_Mylast = _Myfirst + count;
		memory_type::memcpy(old_first, old_first + min(count, old_size), _Myfirst);
		memory_type::free(old_first);
	}

	void swap(_Vector_t& rhs) noexcept
	{
		_swap_all(rhs);
	}

public:

	_HDINLINE_ constexpr bool empty() const noexcept { return _Myfirst == _Mylast; }
	_HDINLINE_ constexpr size_type size() const noexcept { return _Mylast - _Myfirst; }
	_HDINLINE_ pointer data() noexcept { return _Myfirst; }
	_HDINLINE_ const_pointer data() const noexcept { return _Myfirst; }
	_HDINLINE_ iterator begin() noexcept { return _Myfirst; }
	_HDINLINE_ const_iterator begin() const noexcept { return _Myfirst; }
	_HDINLINE_ iterator end() noexcept { return _Mylast; }
	_HDINLINE_ const_iterator end() const noexcept { return _Mylast; }
	_HDINLINE_ reference front() noexcept { return *_Myfirst; }
	_HDINLINE_ const_reference front() const noexcept { return *_Myfirst; }
	_HDINLINE_ reference back() noexcept { return *(_Mylast - 1); }
	_HDINLINE_ const_reference back() const noexcept { return *(_Mylast - 1); }
	_HDINLINE_ reference operator[](size_type i) noexcept { return *(_Myfirst + i); }
	_HDINLINE_ const_reference operator[](size_type i) const noexcept { return *(_Myfirst + i); }
};

} // Namespace Detail

template <typename, unsigned = 4, typename = void> class Array_t;

// Host Array
template <typename T, unsigned N>
class Array_t<T, N, _Detail::_not_portable_t<T>> : 
	public _Detail::_Vector_t<T, _Detail::_host_memory<T>>
{
private:

	using Host = _Detail::_Vector_t<T, _Detail::_host_memory<T>>;

public:

	using pointer = typename Host::pointer;
	using const_pointer = typename Host::const_pointer;
	using reference = typename Host::reference;
	using const_reference = typename Host::const_reference;
	using size_type = typename Host::size_type;
	using difference_type = typename Host::difference_type;
	using iterator = typename Host::iterator;
	using const_iterator = typename Host::const_iterator;
	using index_type = long long;

private:
		
	using Stride_t = _Detail::_Static_array<size_t, N>;
	template <typename... Ts> using call_if_t = _Detail::_call_if_t<N, Ts...>;

private:

	Stride_t stride = { 0, };

public: // Constructor && Destructor

	Array_t() noexcept = default;
	
	template <typename... Ts, typename = call_if_t<size_type, Ts...>>
	explicit Array_t(size_type first, Ts... pack) :
		Host{ (first * ... * pack) },
		stride{ first, first *= pack... }
	{}

	template <typename... Ts, typename = call_if_t<size_type, Ts...>>
	explicit Array_t(const_pointer ptr, size_type first, Ts... pack) :
		Array_t{ first, pack... }
	{
		Host::memory_type::memcpy(ptr, ptr + Host::size(), Host::begin());
	}

	Array_t(const Array_t& rhs) :
		Host{ rhs },
		stride{ rhs.stride }
	{}

	Array_t(Array_t&& rhs) :
		Host{ std::move(rhs) },
		stride{ std::move(rhs.stride) }
	{}

public:

	Array_t& operator=(const Array_t& rhs)
	{
		Host::operator=(rhs);
		stride = rhs.stride;
		return *this;
	}

	Array_t& operator=(Array_t&& rhs)
	{
		Host::operator=(std::move(rhs));
		stride = std::move(stride);
		return *this;
	}

public:

	void Destroy() { Host::clear(); }

	~Array_t() { stride = { 0, }; }
	
public:

	template <typename... Ts, typename = call_if_t<size_type, Ts...>>
	void Create(size_type first, Ts... pack)
	{
		stride = { first, first *= pack... };
		*(static_cast<Host*>(this)) = Host(stride[sizeof...(Ts)]);
	}

public:

	void Copy(const_pointer first, const_pointer last = nullptr)
	{
		Host::memory_type::memcpy(first, last ? last : first + Host::size(), Host::begin());
	}

public:

	void Fill(const_reference val)
	{
		Host::memory_type::memset(Host::begin(), Host::end(), val);
	}

public:

	template <typename... Ts, typename = call_if_t<size_type, Ts...>>
	void Reshape(size_type first, Ts... pack)
	{
#ifdef _DEBUG
		size_type _n = (first * ... * pack);
		if (!Host::empty()) assert(_n == Host::size());
#endif
		stride = { first, first *= pack...};
	}

public:

	template <unsigned D> constexpr size_type Rank() const noexcept
	{
		static_assert(D <= N && D);
		if constexpr (D > 1) return stride[D - 1] / stride[D - 2];
		else return stride[0];
	}

	template <unsigned D> constexpr size_type Stride() const noexcept
	{
		static_assert(D <= N && D);
		if constexpr (D > 1) return stride[D - 2];
		else return 1;
	}

	template <typename... Ts, typename = call_if_t<index_type, Ts...>>
	inline reference operator() (index_type idx, Ts... seq) noexcept
	{
		if constexpr (sizeof...(Ts) > 0) {
			int i = 0; (..., (idx += seq * stride[i++]));
		}
		return Host::operator[](idx);
	}

	template <typename... Ts, typename = call_if_t<index_type, Ts...>>
	inline const_reference operator() (index_type idx, Ts... seq) const noexcept
	{
		if constexpr (sizeof...(Ts) > 0) {
			int i = 0; (..., (idx += seq * stride[i++]));
		}
		return Host::operator[](idx);
	}

}; // Host Array

template <typename T, unsigned N>
class Array_t<T, N, _Detail::_portable_t<T>> : 
	public _Detail::_Vector_t<T, _Detail::_host_memory<T>>,
	private _Detail::_Vector_t<T, _Detail::_device_memory<T>>
{
private:

	using Host = _Detail::_Vector_t<T, _Detail::_host_memory<T>>;
	using Device = _Detail::_Vector_t<T, _Detail::_device_memory<T>>;

public:

	using pointer = typename Host::pointer;
	using const_pointer = typename Host::const_pointer;
	using reference = typename Host::reference;
	using const_reference = typename Host::const_reference;
	using size_type = typename Host::size_type;
	using difference_type = typename Host::difference_type;
	using iterator = typename Host::iterator;
	using const_iterator = typename Host::const_iterator;
	using index_type = long long;

private:

	using Stride_t = _Detail::_Static_array<size_t, N>;
	template <typename... Ts> using call_if_t = _Detail::_call_if_t<N, Ts...>;

private:

	Stride_t stride = { 0, };

public:

	Array_t() noexcept = default;

	template <typename... Ts, typename = call_if_t<size_type, Ts...>>
	explicit Array_t(size_type first, Ts... pack) :
		Host{ (first * ... * pack) },
		stride{ first, first *= pack... }
	{}

	template <typename... Ts, typename = call_if_t<size_type, Ts...>>
	explicit Array_t(const_pointer ptr, size_type first, Ts... pack) :
		Array_t{ first, pack... }
	{
		Host::memory_type::memcpy(ptr, ptr + Host::size(), Host::begin());
	}

	Array_t(const Array_t& rhs) :
		Host{ rhs },
		Device{ rhs },
		stride{ rhs.stride }
	{}

	Array_t(Array_t&& rhs) :
		Host{ std::move(rhs) },
		Device{ std::move(rhs) },
		stride{ std::move(rhs.stride) }
	{}

public:
	
	Array_t& operator=(const Array_t& rhs)
	{
		Host::operator=(rhs);
		Device::operator=(rhs);
		stride = rhs.stride;
		return *this;
	}

	Array_t& operator=(Array_t&& rhs)
	{
		Host::operator=(std::move(rhs));
		Device::operator=(std::move(rhs));
		stride = std::move(rhs.stride);
		return *this;
	}

public:

	void Destroy() { Host::clear(); }
	void DestroyDevice() { Device::clear(); }

	~Array_t() { stride = { 0, }; }

public:

	template <typename... Ts, typename = call_if_t<size_type, Ts...>>
	void Create(size_type first, Ts... pack)
	{
		stride = { first, first *= pack... };
		*(static_cast<Host*>(this)) = Host(stride[sizeof...(Ts)]);
	}

	template <typename... Ts, typename = call_if_t<size_type, Ts...>>
	void CreateDevice(size_type first, Ts... pack)
	{
		stride = { first, first *= pack... };
		*(static_cast<Device*>(this)) = Device(stride[sizeof...(Ts)]);
	}

public:

	void CopyHtoH(const_pointer first, const_pointer last = nullptr)
	{
		Host::memory_type::memcpy(first, last ? last : first + Host::size(), Host::begin());
	}
	void CopyDtoD(const_pointer first, const_pointer last = nullptr)
	{
		cudaCheckError(
			cudaMemcpy(Device::data(), first, last ? last - first : Device::size() * sizeof(T),
				cudaMemcpyDeviceToDevice)
		);
	}
	void CopyHtoD(const_pointer first, const_pointer last = nullptr)
	{
		cudaCheckError(
			cudaMemcpy(Device::data(), first, last ? last - first : Device::size() * sizeof(T),
				cudaMemcpyHostToDevice)
		);
	}
	void CopyDtoH(const_pointer first, const_pointer last = nullptr)
	{
		cudaCheckError(
			cudaMemcpy(Host::data(), first, last ? last - first : Host::size() * sizeof(T),
				cudaMemcpyDeviceToHost)
		);
	}
	void CopyHtoD()
	{
		Device::resize(Host::size());
		cudaCheckError(
			cudaMemcpy(Device::data(), Host::data(), Device::size() * sizeof(T),
				cudaMemcpyHostToDevice)
		);
	}
	void CopyDtoH()
	{
		Host::resize(Device::size());
		cudaCheckError(
			cudaMemcpy(Host::data(), Device::data(), Host::size() * sizeof(T),
				cudaMemcpyDeviceToHost)
		);
	}

public:

	void Fill(const_reference val)
	{
		Host::memory_type::memset(Host::begin(), Host::end(), val);
	}
	void FillDevice(const_reference val)
	{
		Device::memory_type::memset(Device::data(), Device::end(), val);
	}

public:

	template <typename... Ts, typename = call_if_t<size_type, Ts...>>
	void Reshape(size_type first, Ts... pack)
	{
#ifdef _DEBUG
		size_type _n = (first * ... * pack);
		if (!Host::empty()) assert(_n == Host::size());
		if (!Device::empty()) assert(_n == Device::size());
#endif
		stride = { first, first *= pack... };
	}

public:

	_HDINLINE_ constexpr bool device_empty() const noexcept { return Device::empty(); }
	_HDINLINE_ constexpr size_type device_size() const noexcept { return Device::size(); }
	_HDINLINE_ pointer device_data() noexcept { return Device::data(); }
	_HDINLINE_ const_pointer device_data() const noexcept { return Device::data(); }

	_HDINLINE_ constexpr bool empty() const noexcept { return Host::empty(); }
	_HDINLINE_ constexpr size_type size() const noexcept { return Host::size(); }
	_HDINLINE_ pointer data() noexcept { return Host::data(); }
	_HDINLINE_ const_pointer data() const noexcept { return Host::data(); }
	_HDINLINE_ iterator begin() noexcept { return Host::begin(); }
	_HDINLINE_ const_iterator begin() const noexcept { return Host::begin(); }
	_HDINLINE_ iterator end() noexcept { return Host::end(); }
	_HDINLINE_ const_iterator end() const noexcept { return Host::end(); }
	_HDINLINE_ reference front() noexcept { return Host::front(); }
	_HDINLINE_ const_reference front() const noexcept { return Host::front(); }
	_HDINLINE_ reference back() noexcept { return Host::back(); }
	_HDINLINE_ const_reference back() const noexcept { return Host::back(); }

	template <unsigned D> 
	_HDINLINE_ constexpr size_type Rank() const noexcept
	{
		static_assert(D <= N && D);
		if constexpr (D > 1) return stride[D - 1] / stride[D - 2];
		else return stride[0];
	}

	template <unsigned D> 
	_HDINLINE_ constexpr size_type Stride() const noexcept
	{
		static_assert(D <= N && D);
		if constexpr (D > 1) return stride[D - 2];
		else return 1;
	}

	template <typename... Ts, typename = call_if_t<index_type, Ts...>>
	_HDINLINE_ reference operator() (index_type idx, Ts... seq) noexcept
	{
		if constexpr (sizeof...(Ts) > 0) {
			int i = 0; (..., (idx += seq * stride[i++]));
		}
#ifdef __CUDA_ARCH__
		return Device::operator[](idx);
#else
		return Host::operator[](idx);
#endif
	}

	template <typename... Ts, typename = call_if_t<index_type, Ts...>>
	_HDINLINE_ const_reference operator() (index_type idx, Ts... seq) const noexcept
	{
		if constexpr (sizeof...(Ts) > 0) {
			int i = 0; (..., (idx += seq * stride[i++]));
		}
#ifdef __CUDA_ARCH__
		return Device::operator[](idx);
#else
		return Host::operator[](idx);
#endif
	}

}; // Device Array

}

#undef _HDINLINE_