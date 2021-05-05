#include "Input.h"
#include "IOUtil.h"
#include "Exception.h"
#include "IOUtil.h"
#include "Lexer.h"

namespace IO
{

using Parse = Util_t;
using SC = Parse::SC;

void InputManager_t::HashTree_t::ProcessMacro(const string& contents)
{
	string::size_type macro_beg, pos = 0;

	while ((macro_beg = contents.find(SC::HASHTAG, pos)) != string::npos) {
		auto macro_end = Parse::FindEndOfMacro(contents, macro_beg);
		auto macro = contents.substr(macro_beg, macro_end - macro_beg);
		auto v = Parse::ExtractMacro(macro);
		this->macro[v[0]][v[1]] = v[2];
	}
}

void InputManager_t::HashTree_t::Make(const string& file, 
	string::size_type Beg, string::size_type End)
{
	string::size_type beg = Beg, end;

	while ((end = Parse::FindEndPoint(file, beg)) < End) {
		auto name = file.substr(Beg, beg - Beg);

		beg = file.find_first_not_of(SC::LBRACE, beg);
		auto contents = file.substr(beg, end - beg);

		HashTree_t Tnew, &T = this->children[Parse::Trim(name)];
		Tnew.parent = this;
		Tnew.CountLine(name, contents);
		T = std::move(Tnew);

		if (Parse::IsClosed(contents)) 
			T.Make(file, beg, end);
		else 
			T.contents = contents;

		beg = end + 1; Beg = beg;
	}
}

void InputManager_t::HashTree_t::CountLine(const string& name, const string& contents)
{
	num_lines = Parse::LineCount(name) + Parse::LineCount(contents);
	line_info = Parse::LineCount(name, name.find_first_not_of(SC::LF)) + 1;
	line_info += parent->line_info;
	for (const auto& t : parent->children) line_info += t.second.num_lines;
}

InputManager_t::Blocks InputManager_t::GetBlockID(string line) const
{
	int pos_end = line.find(SC::BLANK, 0);
	string block = line.substr(0, pos_end);

	block = Parse::Uppercase(block);
	for (int i = 0; i < BlockNames.size(); ++i)
		if (!block.compare(BlockNames[i]))
			return static_cast<Blocks>(i);
	return Blocks::INVALID;
}

template <typename T>
T InputManager_t::GetCardID(Blocks block, string line) const
{
	static constexpr T INVALID = static_cast<T>(-1);

	int pos_beg = line.find_first_not_of(SC::BLANK);
	if (pos_beg == string::npos) return INVALID;

	int pos_end = line.find(SC::BLANK, pos_beg);
	string card = line.substr(pos_beg, pos_end - pos_beg);

	card = Parse::Uppercase(card);
	int b = static_cast<int>(block);
	for (int i = 0; i < CardNames[b].size(); ++i) 
		if (!card.compare(CardNames[b][i])) 
			return static_cast<T>(i);
	return INVALID;
}

void InputManager_t::ExtractInput(istream& fin)
{
	stringstream strstream;

	strstream << fin.rdbuf();
	
	original = strstream.str();

	std::replace(original.begin(), original.end(), SC::TAB, SC::BLANK);
	std::replace(original.begin(), original.end(), SC::CR, SC::BLANK);
}

void InputManager_t::InspectSyntax(const string& contents)
{
	using Except = Exception_t;
	using Code = Except::Code;

	// Check parenthesis, braces and brackets.
	try {
		Parse::AreBracketsMatched(contents);
	}
	catch (const exception& e) {
		Except::Abort(Code::MISMATCHED_BRAKETS, e.what());
	}

	// Check macro definition
	try {
		Parse::IsMacroValid(contents);
	}
	catch (const exception& e) {
		Except::Abort(Code::INVALID_MACRO, e.what());
	}

	// Check variables : block, card and unit geometry
	try {
		Parse::IsVariableCorrect(contents);
	}
	catch (const exception& e) {
		Except::Abort(Code::INVALID_VARIABLE, e.what());
	}
}

// BLOCK

void InputManager_t::ParseGeometryBlock(HashTree_t& Tree)
{
	using Except = Exception_t;
	using Code = Except::Code;
	using Cards = GeometryCards;

	for (auto& T : Tree.children) {
		
		auto& object = T.second;
		string card = Parse::Trim(T.first);
		Cards ID = GetCardID<Cards>(Blocks::GEOMETRY, card);

		switch (ID)
		{
		case Cards::UNITVOLUME :
			ParseUnitVolumeCard(object); break;
		case Cards::UNITCOMP :
			ParseUnitCompCard(object); break;
		case Cards::INVALID :
			Except::Abort(Code::INVALID_INPUT_CARD, card, object.GetLineInfo());
		}

	}

}

void InputManager_t::ParseMaterialBlock(HashTree_t& Tree)
{

}

void InputManager_t::ParseOptionBlock(HashTree_t& Tree)
{

}

/// GEOMETRY CARDS

void InputManager_t::ParseUnitVolumeCard(HashTree_t& Tree)
{
	using Except = Exception_t;
	using Code = Except::Code;

	string ORIGIN = "ORIGIN";

	for (auto& T : Tree.children) {
		auto name = Parse::Trim(T.first);
		auto& object = T.second;
		auto contents = Parse::EraseSpace(object.contents);
		auto v = Parse::SplitFields(contents, string(1, SC::SEMICOLON));
		
		UnitVolume_t U;

		for (auto i = v.begin(); i != v.end(); ) {
			auto& s = *i;
			if (s.find(ORIGIN) != string::npos) {
				auto delimiter = string(1, s[ORIGIN.size()]);
				auto u = Parse::SplitFields(s, delimiter);

				if (u.size() != 2) 
					Except::Abort(Code::INVALID_ORIGIN_DATA, object.contents, object.GetLineInfo());
				U.origin = u.back();

				i = v.erase(i);
			}
			else ++i;
		}

		v = Parse::SplitFields(v.front(), SC::DAMPERSAND);

		U.equations = v;

		unitVolumes[name] = U;
	}

}

void InputManager_t::ParseUnitCompCard(HashTree_t& Tree)
{
	using Except = Exception_t;
	using Code = Except::Code;

	for (auto& T : Tree.children) {
		auto& object = T.second;
		auto prefix = Parse::EraseSpace(T.first);
		auto v = Parse::SplitFields(prefix, string(1, SC::COLON));
		if (v.size() != 2) 
			Except::Abort(Code::BACKGROUND_MISSED, object.contents, object.GetLineInfo());
		auto name = v.front();
		auto background = v.back();
		auto contents = Parse::EraseSpace(object.contents);
		v = Parse::SplitFields(contents, string(1, SC::SEMICOLON));

		UnitComp_t U;

		U.background = background;

		for (const auto& s : v) {
			auto u = Parse::SplitFields(s, SC::RDBRACKET);
			U.unitvols.push_back(u.front());
			U.displace.push_back(vector<string>());
			u.erase(u.begin());
			for (const auto& k : u) {
				auto lpos = k.find(SC::LPAREN) + 1;
				auto rpos = k.find(SC::RPAREN);
				if (lpos == string::npos && rpos == string::npos)
					U.displace.back().push_back(k);
				else if (lpos != string::npos && rpos != string::npos)
					U.displace.back().push_back(k.substr(lpos, rpos - lpos));
				else
					Except::Abort(Code::MISMATCHED_BRAKETS, object.contents, object.GetLineInfo());
			}
		}

		unitComps[name] = U;
	}
}

void InputManager_t::Preprocess()
{
	InspectSyntax(original);

	modified = Parse::ReplaceMacro(original);

	Parse::ReplaceComments(modified);

	HashTree.ProcessMacro(modified);

	modified = Parse::EraseSpace(modified, string(1, SC::BLANK));

	HashTree.Make(modified);
}

void InputManager_t::ReadInput(string file)
{
	using Except = Exception_t;
	using Code = Except::Code;

	Lexer = new Lexer_t;

	ifstream fin(file);
	
	if (fin.fail()) Except::Abort(Code::FILE_NOT_FOUND, file);

	this->file = file;

	ExtractInput(fin);

	fin.close();

	Preprocess();
	
	for (auto& T : HashTree.children) {
		auto& contents = T.second;
		string block = Parse::Trim(T.first);
		Blocks ID = GetBlockID(block);
		switch (ID)
		{
		case Blocks::GEOMETRY :
			ParseGeometryBlock(contents); break;
		case Blocks::MATERIAL :
			ParseMaterialBlock(contents); break;
		case Blocks::OPTION :
			ParseOptionBlock(contents); break;
		case Blocks::INVALID :
			Except::Abort(Code::INVALID_INPUT_BLOCK, block, contents.GetLineInfo());
		}
	}


}

}