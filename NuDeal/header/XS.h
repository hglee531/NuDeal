#pragma once
#include "Defines.h"
#include "Array.h"

namespace XS{
enum class XSType {
	TOT,
	ABS,
	FIS,
	SCAT,
	Nu,
	Chi,
	Kappa,
	CAP,
	N2N,
	N3N
};

struct XSData {
	template<typename T> using Array = LinPack::Array_t<T>;
	array<bool,10> typeXS;
	Array<double> XSval[10];
	Array<double> XSSM[3];
	//double *ScatKernel, *FisKernel;
};

class XSLib {
private:
	bool isMicro;
	int niso, nset, ng;
	int ntemp;
	int scatorder;

	vector<XSData> XS_set;
	
	void Init() { niso = nset = ng = 0; ntemp = 1; }

	void Append(XSData datum) { XS_set.push_back(datum); }

	void Resize() { XS_set.resize(nset); }
public:
	XSLib() { Init(); }

	XSLib(bool isMicro, int ng, int nset, int scatorder = 0);

	XSLib(ifstream libfile);

	const auto& GetXSSet() const { return XS_set; }

	int GetNiso() { return niso; }

	int GetNg() { return ng; }

	int GetNtemp() { return ntemp; }

	int GetScatOrder() { return scatorder; }

	void UploadXSData(array<bool,10> typeXS, vector<vector<vector<double>>> XS, vector<vector<double>> XSSM);
};
}