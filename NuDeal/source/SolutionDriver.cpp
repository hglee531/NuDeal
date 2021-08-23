#include "SolutionDriver.h"

namespace SolutionDriver {
	double SteadyStateDriver::RelativePsiErr() {
		using CompileInfo_t = PhysicalDomain::CompileInfo_t;

		const Array<double> &psi = FSR->GetFisSrc();
		const vector<CompileInfo_t> &compileInfo = FSR->GetCompileInfo();

		double num = 0., denom = 0.;

		for (int i = 0; i < FSR->GetNblocks(); i++) {
			double vol = compileInfo[i].volsum;
			double pow = psi(i)*vol;
			double dpow = pow - prevpsi(i)*vol;
			num += dpow*dpow;
			denom += pow*pow;
		}

		return num / denom;
	}

	void SteadyStateDriver::Initialize(RayTracingDomain &Rays, FluxBoundary &BndFlux, FlatSrcDomain &FSR, FlatXSDomain &FXR, DriverSN &SnDriver) {
		this->Rays = &Rays; this->BndFlux = &BndFlux;
		this->FSR = &FSR; this->FXR = &FXR;
		this->SnDriver = &SnDriver;

		keff = 1.0; nmaxiter = 1000;
		crit_k = 5.e-7; crit_psi = 1.e-5; crit_res = 1.e-6;
		givenpsi = false; givenres = false;		
	}

	void SteadyStateDriver::RunSS() {
		double keff0;

		FSR->InitFluxUnity();
		FSR->UpdateFisSource(*FXR);
		std::copy(FSR->GetFisSrc().begin(), FSR->GetFisSrc().end(), prevpsi.begin());
		FSR->NormalizePsi(keff);
		FSR->UpdateAngularSource(*FXR);


		bool convk, convpsi, convres;
		convk = false;
		convpsi = (false || !givenpsi);
		convres = (false || !givenres);

		BndFlux->InitBndFluxUnity();

		//vector<double> xnormal(BndFlux->GetXBndFlux().size());
		//vector<double> ynormal(BndFlux->GetYBndFlux().size());
		//vector<double> znormal(BndFlux->GetZBndFlux().size());

		for (int iout = 0; iout < nmaxiter; iout++) {
			keff0 = keff;

			//for (int inn = 0; inn < 100; inn++) {
				FSR->InitFluxZero();
				SnDriver->RaySweep();
				//SnDriver->RaySweepRigorous();
			//}
			//std::copy(BndFlux->GetXBndFlux().begin(), BndFlux->GetXBndFlux().end(), xnormal.begin());
			//std::copy(BndFlux->GetYBndFlux().begin(), BndFlux->GetYBndFlux().end(), ynormal.begin());
			//std::copy(BndFlux->GetZBndFlux().begin(), BndFlux->GetZBndFlux().end(), znormal.begin());

			//double fluxnorm = normvec(FSR->GetScalarFlux().data(), FSR->GetScalarFlux().size());
			//double xnorm = normvec(BndFlux->GetXBndFlux().data(), BndFlux->GetXBndFlux().size());
			//double ynorm = normvec(BndFlux->GetYBndFlux().data(), BndFlux->GetYBndFlux().size());
			//double znorm = normvec(BndFlux->GetZBndFlux().data(), BndFlux->GetZBndFlux().size());

			//for (auto iter = xnormal.begin(); iter != xnormal.end(); iter++) (*iter) /= xnorm;
			//for (auto iter = ynormal.begin(); iter != ynormal.end(); iter++) (*iter) /= ynorm;
			//for (auto iter = znormal.begin(); iter != znormal.end(); iter++) (*iter) /= znorm;

			keff = FSR->UpdateFisSource(*FXR);
			FSR->NormalizePsi(keff);
			FSR->UpdateAngularSource(*FXR);

			std::copy(FSR->GetFisSrc().begin(), FSR->GetFisSrc().end(), prevpsi.begin());

			cout << iout << ' ' << keff << ' ' << RelativeKErr(keff, keff0) << endl;

			convk = RelativeKErr(keff, keff0) < crit_k;
			if (givenpsi) convpsi = RelativePsiErr() < crit_psi;

			if (convk && convpsi && convres) break;

			BndFlux->UpdateBoundary(111, 111);
		}
	}
}