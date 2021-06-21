#pragma once
#include "Defines.h"
#include "Library.h"
#include "Array.h"

namespace Library
{

class XSLibrary_t
{
public:
	
	template <typename T> using Array_t = LinPack::Array_t<T>;

	enum class Format
	{
		Micro,
		Macro,
	};

	struct MacroXS_t
	{
		enum class Type
		{
			TR,
			ABS,
			NUFIS,
			KAPPAFIS,
			CHI,
			SCAT,
			COUNT,
		};

		static constexpr int num_types = static_cast<int>(Type::COUNT);
		Array_t<double> xs[num_types];
		Array_t<double>& tr = xs[static_cast<int>(Type::TR)];
		Array_t<double>& abs = xs[static_cast<int>(Type::ABS)];
		Array_t<double>& nufis = xs[static_cast<int>(Type::NUFIS)];
		Array_t<double>& kappafis = xs[static_cast<int>(Type::KAPPAFIS)];
		Array_t<double>& chi = xs[static_cast<int>(Type::CHI)];
		Array_t<double>& scat = xs[static_cast<int>(Type::SCAT)];
	};

	struct MicroXS_t
	{
		enum class Type
		{
			TOT,
			ABS,
			FIS,
			NU,
			KAPPA,
			CHI,
			CAP,
			N2N,
			N3N,
			SCAT,
			COUNT
		};

		static constexpr int num_types = static_cast<int>(Type::COUNT);

		Array_t<double> xs[num_types];
		Array_t<double>& total = xs[static_cast<int>(Type::TOT)];
		Array_t<double>& abs = xs[static_cast<int>(Type::ABS)];
		Array_t<double>& fis = xs[static_cast<int>(Type::FIS)];
		Array_t<double>& nu = xs[static_cast<int>(Type::NU)];
		Array_t<double>& kappa = xs[static_cast<int>(Type::KAPPA)];
		Array_t<double>& chi = xs[static_cast<int>(Type::CHI)];
		Array_t<double>& cap = xs[static_cast<int>(Type::CAP)];
		Array_t<double>& n2n = xs[static_cast<int>(Type::N2N)];
		Array_t<double>& n3n = xs[static_cast<int>(Type::N3N)];
		Array_t<double>& scat = xs[static_cast<int>(Type::SCAT)];
	};


	using MacroType = MacroXS_t::Type;
	using MicroType = MicroXS_t::Type;

private:

	Format format = Format::Micro;
	int num_nuclides = 0, num_groups = 0;
	int num_temps = 0, scat_order = 0;

	std::vector<string> material_labels;
	std::vector<MacroXS_t> MacroXS;
	std::vector<MicroXS_t> MicroXS;

public:

	XSLibrary_t() noexcept = default;

	void ReadMacro(const string& file = "./lib/c5g7.xslib");
	void ReadMicro();

public:

	const auto& GetMicroXS() const { return MicroXS; }

	bool IsMacro() const { return format == Format::Macro; }
	bool IsMicro() const { return format == Format::Micro; }
};

}