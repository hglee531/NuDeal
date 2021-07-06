#pragma once
#include "Defines.h"

namespace Math
{

static constexpr double PI = 3.14159265358979323846264338327950288419716939937510582097494459230781;

template <typename T> decltype(auto) square(T a)
{
	return a * a;
}

template <typename T> decltype(auto) cubic(T a)
{
	return a * a * a;
}

}