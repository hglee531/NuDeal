#pragma once
#include "Defines.h"
#include "IODeclare.h"

namespace IO
{
// Special Characters
struct SpecialCharacters
{
	static constexpr char LeftParen = '(';
	static constexpr char RightParen = ')';
	static constexpr char LeftBracket = '[';
	static constexpr char RightBracket = ']';
	static constexpr char LeftBrace = '{';
	static constexpr char RightBrace = '}';
	static constexpr char LeftAngle = '<';
	static constexpr char RightAngle = '>';
	static constexpr char Equal = '=';
	static constexpr char Plus = '+';
	static constexpr char Minus = '-';
	static constexpr char Asterisk = '*';
	static constexpr char Slash = '/';
	static constexpr char Caret = '^';
	static constexpr char Hash = '#';
	static constexpr char Dot = '.';
	static constexpr char Comma = ',';
	static constexpr char SemiColon = ';';
	static constexpr char Colon = ':';
	static constexpr char BackSlash = '\\';
	static constexpr char Blank = ' ';
	static constexpr char Tab = '\t';
	static constexpr char CR = '\r';
	static constexpr char LF = '\n';
};

class Util_t
{
public :

	using SC = SpecialCharacters;
	using size_type = string::size_type;

	// String Parsing Utility
	/// Type Conversion
	static string Uppercase(const string& line);
	static int Integer(string field);
	static double Float(string field);
	static bool Logical(string field);
	/// String Manipulation
	static string Trim(const string& field, const string& delimiter = "\n ");
	static int LineCount(const string& line, size_type count = string::npos);
	static string EraseSpace(const string& line, const string& delimiter = "\n ");
	static size_type FindEndPoint(const string& contents, size_type& pos);
	static vector<string> SplitFields(const string& line, const string& delimiter);
	static vector<string> SplitFields(const string& line, char delimiter);
	static double3 GetCoordinate(const string& field);
	/// Filestream Manipulation
	static string GetLine(istream& fin);
	static int FindKeyword(istream& fin, const string& keyword);
	
	static bool IsClosed(const string& s);
};

}