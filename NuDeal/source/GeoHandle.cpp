#include "GeoHandle.h"

namespace Geometry
{
void splittree::Assign(int id) {
	subnode = nullptr; uppernode = nullptr;
	this->id = id; nleaf = 0;
}

void splittree::Assign(int id, splittree &uppernode) {
	subnode = nullptr; this->uppernode = &uppernode;
	this->id = id; nleaf = 0;
}

void splittree::Branching(int nleaf) {
	subnode = new splittree[nleaf];
	this->nleaf = nleaf;
	for (int i = 0; i < nleaf; i++) {
		subnode[i].Assign(i, *this);
	}
}

void splittree::RecordNodeInfo(int idvol, double volcm3, double3 midpt) {
	thisinfo.idvol = idvol;
	thisinfo.volcm3 = volcm3;
	thisinfo.midpt = midpt;
}

inline bool InLocalBox(CartAxis axis, double3 ptL, double3 ptR, double3 aninter) {
	switch (axis) {
	case CartAxis::X:
		if (aninter.x > ptL.x + eps_geo && aninter.x < ptR.x - eps_geo) return true;
		break;
	case CartAxis::Y:
		if (aninter.y > ptL.y + eps_geo && aninter.y < ptR.y - eps_geo) return true;
		break;
	case CartAxis::Z:
		if (aninter.z > ptL.z + eps_geo && aninter.z < ptR.z - eps_geo) return true;
		break;
	}
	return false;
}

void GeometryHandler::init() {
	nvol = nnode = divlevel = 0; isfinal = false; issetord = false;
	x0 = y0 = z0 = Lx = Ly = Lz = 0.0;
}

int GeometryHandler::FindVolId(double3 ptL, double3 ptR, bool lowest) {
	int id = -1;
	
	// Store volumes whose bound boxes cross the finite box
	vector<int> VolInBound;
	for (int i = 0; i < nvol; i++) {
		int inboundcount = 0;
		double OL, OR;
		double Ll, Lr, Rl, Rr;
		double Ldet, Rdet;
		OL = BdL[i].x; OR = BdR[i].x;
		Ll = ptL.x - OL; Lr = ptR.x - OL;
		Rl = ptL.x - OR; Rr = ptR.x - OR;
		if (Ll*Lr < 0 && Rl*Rr < 0) inboundcount++;
		if (Lr*Rr > 0 && Ll*Rl > 0 && Ll*Lr > 0) continue;
		if (mode != Dimension::OneD) {
			OL = BdL[i].y; OR = BdR[i].y;
			Ll = OL - ptL.y; Lr = OL - ptR.y;
			Rl = OR - ptL.y; Rr = OR - ptR.y;
			if (Ll*Lr < 0 && Rl*Rr < 0) inboundcount++;
			if (Lr*Rr > 0 && Ll*Rl > 0 && Ll*Lr > 0) continue;
			if (mode == Dimension::ThreeD) {
				OL = BdL[i].z; OR = BdR[i].z;
				Ll = OL - ptL.z; Lr = OL - ptR.z;
				Rl = OR - ptL.z; Rr = OR - ptR.z;
				if (Ll*Lr < 0 && Rl*Rr < 0) inboundcount++;
				if (Lr*Rr > 0 && Ll*Rl > 0 && Ll*Lr > 0) continue;
				if (inboundcount == 3) return id;
			}
			else {
				if (inboundcount == 2) return id;
			}
		}
		else {
			if (inboundcount == 1) return id;
		}
		VolInBound.push_back(i);
	}

	int nbnd = VolInBound.size();
	if (!lowest) {
		for (int i = 0; i < nbnd; i++) {
			vector<double> interpts;
			int idvol = VolInBound[i];
			array<double, 2> val;

			// Along x axis
			double3 aninter;
			aninter.y = ptL.y; aninter.z = ptL.z;
			val[0] = ptL.y; val[1] = ptL.z;
			int ninter = Volumes[idvol].GetIntersection(CartPlane::YZ, val, interpts);
			for (int j = 0; j < ninter; j++) {
				aninter.x = interpts[j];
				if (InLocalBox(CartAxis::X, ptL, ptR, aninter)) return id;
			}
			if (mode != Dimension::OneD) {
				aninter.y = ptR.y; val[0] = ptR.y;
				ninter = Volumes[idvol].GetIntersection(CartPlane::YZ, val, interpts);
				for (int j = 0; j < ninter; j++) {
					aninter.x = interpts[j];
					if (InLocalBox(CartAxis::X, ptL, ptR, aninter)) return id;
				}
			}
			if (mode == Dimension::ThreeD) {
				aninter.y = ptL.y; aninter.z = ptR.z;
				val[0] = ptL.y; val[1] = ptR.z;
				int ninter = Volumes[idvol].GetIntersection(CartPlane::YZ, val, interpts);
				for (int j = 0; j < ninter; j++) {
					aninter.x = interpts[j];
					if (InLocalBox(CartAxis::X, ptL, ptR, aninter)) return id;
				}
				aninter.y = ptR.y; val[0] = ptR.y;
				ninter = Volumes[idvol].GetIntersection(CartPlane::YZ, val, interpts);
				for (int j = 0; j < ninter; j++) {
					aninter.x = interpts[j];
					if (InLocalBox(CartAxis::X, ptL, ptR, aninter)) return id;
				}
			}


			// Along y axis
			if (mode != Dimension::OneD) {
				aninter.x = ptL.x; aninter.z = ptL.z;
				val[0] = ptL.x; val[1] = ptL.z;
				ninter = Volumes[idvol].GetIntersection(CartPlane::XZ, val, interpts);
				for (int j = 0; j < ninter; j++) {
					aninter.y = interpts[j];
					if (InLocalBox(CartAxis::Y, ptL, ptR, aninter)) return id;
				}
				aninter.x = ptR.x; val[0] = ptR.x;
				ninter = Volumes[idvol].GetIntersection(CartPlane::XZ, val, interpts);
				for (int j = 0; j < ninter; j++) {
					aninter.y = interpts[j];
					if (InLocalBox(CartAxis::Y, ptL, ptR, aninter)) return id;
				}
			}
			if (mode == Dimension::ThreeD) {
				aninter.x = ptL.x; aninter.z = ptR.z;
				val[0] = ptL.x; val[1] = ptR.z;
				int ninter = Volumes[idvol].GetIntersection(CartPlane::XZ, val, interpts);
				for (int j = 0; j < ninter; j++) {
					aninter.y = interpts[j];
					if (InLocalBox(CartAxis::Y, ptL, ptR, aninter)) return id;
				}
				aninter.x = ptR.x; val[0] = ptR.x;
				ninter = Volumes[idvol].GetIntersection(CartPlane::XZ, val, interpts);
				for (int j = 0; j < ninter; j++) {
					aninter.y = interpts[j];
					if (InLocalBox(CartAxis::Y, ptL, ptR, aninter)) return id;
				}
			}


			// Along z axis
			if (mode == Dimension::ThreeD) {
				aninter.x = ptL.x; aninter.y = ptL.y;
				val[0] = ptL.x; val[1] = ptL.y;
				ninter = Volumes[idvol].GetIntersection(CartPlane::XY, val, interpts);
				for (int j = 0; j < ninter; j++) {
					aninter.z = interpts[j];
					if (InLocalBox(CartAxis::Z, ptL, ptR, aninter)) return id;
				}
				aninter.x = ptR.x; val[0] = ptR.x;
				ninter = Volumes[idvol].GetIntersection(CartPlane::XY, val, interpts);
				for (int j = 0; j < ninter; j++) {
					aninter.z = interpts[j];
					if (InLocalBox(CartAxis::Z, ptL, ptR, aninter)) return id;
				}
				aninter.x = ptL.x; aninter.y = ptR.y;
				val[0] = ptL.x; val[1] = ptR.y;
				int ninter = Volumes[idvol].GetIntersection(CartPlane::XY, val, interpts);
				for (int j = 0; j < ninter; j++) {
					aninter.z = interpts[j];
					if (InLocalBox(CartAxis::Z, ptL, ptR, aninter)) return id;
				}
				aninter.y = ptR.y; val[0] = ptR.y;
				ninter = Volumes[idvol].GetIntersection(CartPlane::XY, val, interpts);
				for (int j = 0; j < ninter; j++) {
					aninter.z = interpts[j];
					if (InLocalBox(CartAxis::Z, ptL, ptR, aninter)) return id;
				}
			}
		}
	}

	// When alive through above loops, this box has only one volume in it.
	// Find the volume index
	double xmid = 0.5*(ptL.x + ptR.x), ymid = 0.5*(ptL.y + ptR.y), zmid = 0.5*(ptL.z + ptR.z);
	for (int i = 0; i < nbnd; i++) {
		int idvol = VolInBound[i];
		//if (Volumes[idvol].IsInside(xmid, ymid, zmid, false)) return idvol;
		if (!Volumes[idvol].IsInside(xmid, ymid, zmid, true)) VolInBound[i] = -1;
	}
	int imax = -100; double volmin = 1.e+50;
	for (int i = 0; i < nbnd; i++) {
		int idvol = VolInBound[i];
		if (idvol > -1) {
			if (volmin > Volumes[idvol].GetVolume()) {
				imax = idvol; volmin = Volumes[idvol].GetVolume();
			}
		}
	}

	// Still alive? Error
	return imax;
}

void GeometryHandler::RecursiveSplit(double3 ptL, double3 ptR, int thisLv, splittree &thisnode) {
	double lx = lx0 / pow3[thisLv], ly = ly0, lz = lz0;
	bool lowest = false;
	int nx = 3, ny = 1, nz = 1;
	lowest = (lx / 3.0 < maxtau);
	if (mode != Dimension::OneD) {
		ly = ly0 / pow3[thisLv];
		lowest = lowest || (ly / 3.0 < maxtau);
		ny = 3;
	}
	if (mode == Dimension::ThreeD) {
		lz = lz0 / pow3[thisLv];
		lowest = lowest || (lz / 3.0 < maxtau);
		nz = 3;
	}

	int idvol = FindVolId(ptL, ptR, lowest);
	if (idvol == -1) {
		thisnode.Branching(nx*ny*nz);
		splittree* subnodes = thisnode.GetPtrSubnode();
		double3 ptL1 = ptL, ptR1 = ptR;
		lx /= 3.0; 
		if (mode != Dimension::OneD) ly /= 3.0;
		if (mode == Dimension::ThreeD) lz /= 3.0;
		int subid = 0;
		for (int k = 0; k < nz; k++) {
			if (nz > 1) { ptL1.z = ptL.z + (double)k*lz; ptR1.z = ptL.z + (double)(k + 1)*lz; }
			for (int j = 0; j < ny; j++) {
				if (ny > 1) { ptL1.y = ptL.y + (double)j*ly; ptR1.y = ptL.y + (double)(j + 1)*ly; }
				for (int i = 0; i < nx; i++) {
					ptL1.x = ptL.x + (double)i*lx; ptR1.x = ptL.x + (double)(i + 1)*lx;
					RecursiveSplit(ptL1, ptR1, (thisLv + 1), subnodes[subid]);
					subid++;
				}
			}
		}
	}
	//vector<int> ids{ idvol };
	//vector<double> vol{ lx*ly*lz };
	double3 midpt;
	midpt.x = 0.5*(ptL.x + ptR.x); midpt.y = 0.5*(ptL.y + ptR.y); midpt.z = 0.5*(ptL.z + ptR.z);
	thisnode.RecordNodeInfo(idvol, lx*ly*lz, midpt);
}

void GeometryHandler::LevelCount(int thisLv, splittree &thisnode) {
	if (thisLv > -1) nnodeLv[thisLv]++;
	splittree *subnode = thisnode.GetPtrSubnode();
	if (subnode != nullptr) {
		int nleaf = thisnode.GetNleaf();
		for (int i = 0; i < nleaf; i++) {
			LevelCount(thisLv + 1, subnode[i]);
		}
	}
}

void GeometryHandler::RecordDiscInfo(int thisLv, splittree &thisnode) {
	int nleaf = thisnode.GetNleaf();
	for (int i = 0; i < nleaf; i++) {
		splittree *subnodes = thisnode.GetPtrSubnode();
		RecordDiscInfo(thisLv + 1, subnodes[i]);
	}

	if (thisLv > -1) {
		int thisLvId = nnodeLv[thisLv];
		nnodeLv[thisLv]++;
		if (thisLv > 0) upperdivmap[thisLv][thisLvId] = nnodeLv[thisLv - 1];
		info[thisLv][thisLvId] = thisnode.GetNodeInfo();
	}
}

GeometryHandler::GeometryHandler(double origin[3], double L[3]) {
	nvol = nnode = divlevel = 0;
	SetOrdinates(origin, L);
}

void GeometryHandler::SetOrdinates(double origin[3], double L[3]) {
	issetord = true;
	x0 = origin[0]; y0 = origin[1]; z0 = origin[2];
	Lx = L[0]; Ly = L[1]; Lz = L[2];	
}

void GeometryHandler::append(UnitComp &acomp) {
	UVbuf.push(acomp);
	matidbuf.push(acomp.GetMatId());
}

void GeometryHandler::FinalizeVolumes() {
	nvol = UVbuf.size();
	Volumes.resize(nvol); imat.resize(nvol);
	BdL.resize(nvol); BdR.resize(nvol);
	for (int i = 0; i < nvol; i++) {
		Volumes[i] = UVbuf.front();
		imat[i] = matidbuf.front();
		UVbuf.pop(); matidbuf.pop();
		double2 xs, ys, zs;
		Volumes[i].GetBoundBox(xs, ys, zs);
		BdL[i].x = xs.x; BdL[i].y = ys.x; BdL[i].z = zs.x;
		BdR[i].x = xs.y; BdR[i].y = ys.y; BdR[i].z = zs.y;
	}
	if (!issetord) {
		x0 = BdL[0].x; y0 = BdL[0].y; z0 = BdL[0].z;
		Lx = BdR[0].x - x0; Ly = BdR[0].y - y0; Lz = BdR[0].z - z0;
		for (int i = 1; i < nvol; i++) {
			double xl = BdL[i].x, yl = BdL[i].y, zl = BdL[i].z;
			x0 = (x0 > xl) ? xl : x0;
			y0 = (y0 > yl) ? yl : y0;
			z0 = (z0 > zl) ? zl : z0;
			Lx = (Lx > BdR[i].x - xl) ? Lx : BdR[i].x - xl;
			Ly = (Ly > BdR[i].y - yl) ? Ly : BdR[i].y - yl;
			Lz = (Lz > BdR[i].z - zl) ? Lz : BdR[i].z - zl;
		}
		issetord = true;
	}
}

bool GeometryHandler::Discretize(Dimension Dim, double minlen, double maxlen, int multiplier) {
	mode = Dim; mintau = minlen * multiplier; maxtau = maxlen;
	for (int i = 0; i < nvol; i++) {
		double lxm = BdR[i].x - BdL[i].x, lym = BdR[i].y - BdL[i].y, lzm = BdR[i].z - BdL[i].z;
		if (lxm < maxtau) {
			cout << "  *** Warning: Given maximal length is too small. Reset to " << lxm;
			maxtau = lxm;
		}
		if (mode != Dimension::OneD && lym < maxtau) {
			cout << "  *** Warning: Given maximal length is too small. Reset to " << lym;
			maxtau = lym;
		}
		if (mode == Dimension::ThreeD && lzm < maxtau) {
			cout << "  *** Warning: Given maximal length is too small. Reset to " << lzm;
			maxtau = lzm;
		}
	}

	// Define zero-level parameters
	Nx = multiplier * (int)(Lx / mintau + 1); Ny = 1; Nz = 1;
	if (mode != Dimension::OneD) Ny = multiplier * (int)(Ly / mintau + 1);
	if (mode == Dimension::ThreeD) Nz = multiplier * (int)(Lz / mintau + 1);
	lx0 = Lx / (double)Nx; ly0 = Ly / (double)Ny; lz0 = Lz / (double)Nz;

	// Determine the potentially maximal divlevel
	double log3 = log(3.0);
	int potentdivLv = log(lx0 / maxtau) / log3 + 1;
	if (mode != Dimension::OneD) potentdivLv = max(potentdivLv, (int)(log(ly0 / maxtau) / log3) + 1);
	if (mode == Dimension::ThreeD) potentdivLv = max(potentdivLv, (int)(log(lz0 / maxtau) / log3) + 1);
	nnodeLv.resize(potentdivLv); std::fill(nnodeLv.begin(), nnodeLv.end(), 0);
	upperdivmap.resize(potentdivLv); info.resize(potentdivLv);
	for (int i = 0; i < potentdivLv; i++)	pow3.push_back(pow(3.0, i));

	// Ray tracing with the rays parallel to x,y,z axes
	splittree root(0);
	root.Branching(Nx*Ny*Nz);
	splittree *zeronodes = root.GetPtrSubnode();
	double3 ptL, ptR;
	int inode0 = 0;
	ptL.z = z0; ptR.z = z0 + lz0;
	//if (Nz == 1) ptL.z = ptR.z = z0 + 0.5*lz0;
	for (int iz = 0; iz < Nz; iz++) {
		ptL.y = y0; ptR.y = y0 + ly0;
		//if (Ny == 1) ptL.y = ptR.y = y0 + 0.5*ly0;
		for (int iy = 0; iy < Ny; iy++) {
			ptL.x = x0; ptR.x = x0 + lx0;
			//if (Nx == 1) ptL.x = ptR.x = x0 + 0.5*lx0;
			for (int ix = 0; ix < Nx; ix++) {
				RecursiveSplit(ptL, ptR, 0, zeronodes[inode0]);
				inode0++;
				ptL.x += lx0; ptR.x += lx0;
			}
			ptL.y += ly0; ptR.y += ly0;
		}
		ptL.z += lz0; ptR.z += lz0;
	}

	// Defining nnodeLv and upperdivmap
	LevelCount(-1, root);
	for (int i = 0; i < potentdivLv; i++) {
		if (nnodeLv[i] > 0) {
			upperdivmap[i].resize(nnodeLv[i]);
			info[i].resize(nnodeLv[i]);
			std::fill(upperdivmap[i].begin(), upperdivmap[i].end(), 0);
		}
		else break;
		divlevel++;
	}
	std::fill(nnodeLv.begin(), nnodeLv.end(), 0);
	RecordDiscInfo(-1, root);

	nnode = 0;
	for (auto it = info.begin(); it != info.end(); it++) {
		for (auto jt = (*it).begin(); jt != (*it).end(); jt++) {
			if (jt->idvol > -1) nnode++;
		}
	}

	cout << "Discretize Done!" << endl;

	//PrintDiscInfo();

	return true;
}

void GeometryHandler::PrintDiscInfo() const {
	int nnode = 0;
	vector<double> volval(nvol);
	std::fill(volval.begin(), volval.end(), 0.0);
	// Print-out Division Information
	ofstream DiscOut("Discretization.out");
	for (int i = 0; i < divlevel; i++) {
		for (int j = 0; j < nnodeLv[i]; j++) {
			if (info[i][j].idvol < 0) continue;
			nnode++;
			DiscOut << info[i][j].midpt.x << ' ';
			DiscOut << info[i][j].midpt.y << ' ';
			DiscOut << info[i][j].midpt.z << ' ';
			DiscOut << info[i][j].idvol << endl;
			volval[info[i][j].idvol] += info[i][j].volcm3;
		}
	}
	for (int i = 0; i < nvol; i++) {
		cout << "Vol " << i << ": " << volval[i] << endl;
	}
}

void GeometryHandler::GetSizeScalars(int3 &N, int &nnode, int &divlevel, double3 &origin, double3 &Width0) const {
	N.x = Nx; N.y = Ny; N.z = Nz; nnode = this->nnode; divlevel = this->divlevel;
	origin.x = x0; origin.y = y0; origin.z = z0;
	Width0.x = lx0, Width0.y = ly0; Width0.z = lz0;
}

}