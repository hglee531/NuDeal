#pragma once
#define LinearFlux

constexpr double _ZERO = 0.0;
constexpr double2 _ZERO2{ 0.0, 0.0 };

constexpr double _UNITY = 1.0;
constexpr double2 _UNITY2{ 1.0, 0.0 };

#ifdef LinearFlux
using realphi = double2;
constexpr realphi zerophi = _ZERO2;
constexpr realphi unityphi = _UNITY2;
#else
using realphi = double;
constexpr realphi zerophi = _ZERO;
constexpr realphi unityphi = _UNITY;
#endif