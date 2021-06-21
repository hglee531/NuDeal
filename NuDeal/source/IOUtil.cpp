#include "IOUtil.h"
#include "IOExcept.h"

namespace IO
{

using Except = Exception_t;
using Code = Except::Code;

string Util_t::Uppercase(const string& line)
{
	string l = line;
	std::transform(l.begin(), l.end(), l.begin(), ::toupper);
	return l;
}

int Util_t::Integer(string field)
{
	int val;

	try {
		val = stoi(field);
	}
	catch (invalid_argument) {
		Except::Abort(Code::INVALID_INTEGER, field);
	}

	if (val != stod(field)) Except::Abort(Code::INVALID_INTEGER, field);

	return val;
}

double Util_t::Float(string field)
{
	double val;

	try {
		val = stod(field);
	}
	catch (std::invalid_argument) {
		Except::Abort(Code::INVALID_FLOATING_POINT, field);
	}

	return val;
}

bool Util_t::Logical(string field)
{
	Uppercase(field);
	if (!field.compare("T")) return true;
	else if (!field.compare("TRUE")) return true;
	else if (!field.compare("F")) return false;
	else if (!field.compare("FALSE")) return false;

	Except::Abort(Code::INVALID_LOGICAL, field);
}

string Util_t::Trim(const string& field, const string& delimiter)
{
	string s = field;
	auto pos = s.find_first_not_of(delimiter);
	s.erase(0, pos);
	pos = s.find_last_not_of(delimiter) + 1;
	s.erase(pos);

	return s;
}

int Util_t::LineCount(const string& line, size_type count)
{
	auto beg = line.begin();
	auto end = line.end();

	if (count != string::npos) end = beg + count;

	return static_cast<int>(std::count(beg, end, SC::LF));
}

string Util_t::EraseSpace(const string& line, const string& delimiter)
{
	string s = line;
	for (const auto& i : delimiter)
		s.erase(std::remove(s.begin(), s.end(), i), s.end());
	return s;
}

Util_t::size_type Util_t::FindEndPoint(const string& contents, size_type& pos)
{
	pos = contents.find_first_of(SC::LeftBrace, pos);

	size_type end, beg = pos;

	while ((end = contents.find_first_of(SC::RightBrace, beg)) != string::npos) {
		if (IsClosed(contents.substr(pos, end - pos + 1))) break;
		beg = end + 1;
	}
	return end;
}

vector<string> Util_t::SplitFields(const string& line, const string& delimiter)
{
	vector<string> splitted;

	size_type beg, pos = 0;

	while ((beg = line.find_first_not_of(delimiter, pos)) != string::npos) {
		pos = line.find_first_of(delimiter, beg + 1);
		splitted.push_back(line.substr(beg, pos - beg));
	}

	return splitted;
}

vector<string> Util_t::SplitFields(const string& line, char delimiter)
{
	vector<string> splitted;

	size_type beg, pos = 0;

	while ((beg = line.find_first_not_of(delimiter, pos)) != string::npos) {
		pos = line.find_first_of(delimiter, beg + 1);
		splitted.push_back(line.substr(beg, pos - beg));
	}

	return splitted;
}

double3 Util_t::GetCoordinate(const string& field)
{
	auto b = field.find_first_of(SC::LeftParen);
	auto e = field.find_last_of(SC::RightParen);

	if (b == string::npos) b = -1;
	if (e == string::npos) e = field.size();

	auto origin = SplitFields(field.substr(b + 1, e - b - 1), SC::Comma);
	if (origin.size() != 3) Except::Abort(Code::INVALID_COORDINATE, field);
	return make_double3(Float(origin[0]), Float(origin[1]), Float(origin[2]));
}

bool Util_t::IsClosed(const string& s) 
{
	auto lcount = std::count(s.begin(), s.end(), SC::LeftBrace);
	auto rcount = std::count(s.begin(), s.end(), SC::RightBrace);

	return (lcount == rcount) && (lcount > 0);
}

string Util_t::GetLine(istream& fin)
{
	string line;

	do {

		std::getline(fin, line);
		std::replace(line.begin(), line.end(), SC::Tab, SC::Blank);
		std::replace(line.begin(), line.end(), SC::CR, SC::Blank);

		line = Trim(line);

	} while (line.empty());

	return line;
}

int Util_t::FindKeyword(istream& fin, const string& keyword)
{
	int count = 0;
	auto pos = fin.tellg();
	fin.clear(); fin.seekg(ios::beg);
	while (!fin.eof()) {
		string field; fin >> field;
		if (Uppercase(field) == keyword) {
			if (!count) pos = fin.tellg();
			++count;
		}
	}
	fin.clear(); fin.seekg(pos);
	return count;
}

}