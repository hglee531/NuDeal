#include "XSLibrary.h"
#include "IOUtil.h"
#include "IOExcept.h"

namespace Library
{

struct MacroInput_t
{
	enum class type {
		SCAT_ORDER, NG, MATERIAL, COUNT
	};
	inline static const vector<string> cards = {
		"SCAT_ORDER", "NG", "MATERIAL"
	};
};

void XSLibrary_t::ReadMacro(const string& file)
{
	using Util = IO::Util_t;
	using Except = IO::Exception_t;
	using Cards = MacroInput_t::type;

	ifstream fin(file);

	if (fin.fail()) Except::Abort(Except::Code::FILE_NOT_FOUND, file);

	this->format = Format::Macro;

	for (int c = 0; c < static_cast<int>(Cards::COUNT); ++c) {

		int count = Util::FindKeyword(fin, MacroInput_t::cards[c]);
		if (!count) Except::Abort(Except::Code::INVALID_INPUT_CARD);

		switch (static_cast<Cards>(c)) {
		case Cards::SCAT_ORDER:
			fin >> this->scat_order;
			break;

		case Cards::NG:
			fin >> this->num_groups;
			break;

		case Cards::MATERIAL:
			fin.seekg(-MacroInput_t::cards[c].size(), ios::cur);

			material_labels.resize(count);
			MacroXS.resize(count);

			for (int mat = 0; mat < count; ++mat) {
				auto line = Util::GetLine(fin);
				auto fields = Util::SplitFields(line, Util::SC::Blank);
				material_labels[mat] = fields[1];

				for (int t = 0; t < MacroXS_t::num_types - 1; ++t) 
					MacroXS[mat].xs[t].Create(num_groups);
				MacroXS[mat].scat.Create(num_groups, num_groups, scat_order + 1);
				
				for (int g = 0; g < num_groups; ++g)
					for (int t = 0; t < MacroXS_t::num_types - 1; ++t)
						fin >> MacroXS[mat].xs[t](g);

				for (int s = 0; s < scat_order + 1; ++s) 
					for (int g = 0; g < num_groups; ++g) 
						for (int gg = 0; gg < num_groups; ++gg) 
							fin >> MacroXS[mat].scat(gg, g, s);

			}

			break;
		}

	}


	fin.close();
}

}