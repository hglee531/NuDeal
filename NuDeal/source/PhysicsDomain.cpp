#include "PhysicalDomain.h"
#include "XSLibrary.h"

namespace PhysicalDomain {
void BaseDomain::Create(const GeometryHandler &rhs) {
		mode = rhs.GetDimension();

		int3 N;
		double3 origin, Width0;
		rhs.GetSizeScalars(N, nnode, divlevel, origin, Width0);

		Nx = N.x; Ny = N.y; Nz = N.z;
		x0 = origin.x; y0 = origin.y; z0 = origin.z;
		lx0 = Width0.x; ly0 = Width0.y; lz0 = Width0.x;
		nx = 3; ny = (mode != Dimension::OneD) ? 3 : 1; nz = (mode == Dimension::ThreeD) ? 3 : 1;
		nnodeLv.resize(divlevel); upperdivmap.resize(divlevel); lowerdivmap.resize(divlevel - 1);
		serialdivmap.resize(divlevel); nodeInfo.resize(nnode); innodeLv.resize(nnode);

		const vector<int>& nnodeLvRhs = rhs.GetNnodeLv();
		const vector<vector<int>>& upperdivmapRhs = rhs.GetUpperdivmap();
		const vector<vector<NodeInfo_t>>& nodeInfoRhs = rhs.GetNodeinfo();

		int inode = 0;
		for (int i = 0; i < divlevel; i++) {
			nnodeLv[i] = nnodeLvRhs[i];
			upperdivmap[i].resize(nnodeLv[i]);
			serialdivmap[i].resize(nnodeLv[i]);
			std::copy(upperdivmapRhs[i].begin(), upperdivmapRhs[i].end(), upperdivmap[i].begin());
			std::fill(serialdivmap[i].begin(), serialdivmap[i].end(), -1);
			for (int j = 0; j < nnodeLv[i]; j++) {
				if (nodeInfoRhs[i][j].idvol < 0) continue;
				serialdivmap[i][j] = inode;
				nodeInfo[inode] = nodeInfoRhs[i][j];
				innodeLv[inode] = i;
				inode++;
			}
			if (i > 0) {
				int i0 = i - 1;
				lowerdivmap[i0].resize(nnodeLv[i0] + 1);
				std::fill(lowerdivmap[i0].begin(), lowerdivmap[i0].end(), 0);
				for (int j = 0; j < nnodeLv[i]; j++) {
					int jnode = upperdivmap[i][j];
					lowerdivmap[i0][jnode + 1]++;
				}
				for (int j = 0; j < nnodeLv[i0]; j++) {
					lowerdivmap[i0][j + 1] += lowerdivmap[i0][j];
				}
			}
		}

		isalloc = true;
	}

bool ConnectedDomain::FindNeighbor(dir6 srchdir, array<int, 3> ixyz, array<int, 2> LvnId, array<int, 2> &NeighLvnId) {
		int thisLv = LvnId[0], thisId = LvnId[1];
		int ix = ixyz[0], iy = ixyz[1], iz = ixyz[2];
		int mx, my, mz;
		int jx = ix, jy = iy, jz = iz;

		if (thisLv == 0) { mx = Nx; my = Ny; mz = Nz; }
		else { mx = nx; my = ny; mz = nz; }

		// Search upward if an neighboring is not in the local box
		int upid = upperdivmap[thisLv][thisId], nowid = thisId, nowLv = thisLv;
		stack<int> jxs, jys, jzs;
		int jdoff;
		jxs.push(jx); jys.push(jy); jzs.push(jz);
		while (true) {
			bool upfurther;
			switch (srchdir) {
			case xleft:
				upfurther = (jx == 0);
				break;
			case yleft:
				upfurther = (jy == 0);
				break;
			case zleft:
				upfurther = (jz == 0);
				break;
			case xright:
				upfurther = (jx + 1 == mx);
				break;
			case yright:
				upfurther = (jy + 1 == my);
				break;
			case zright:
				upfurther = (jz + 1 == mz);
				break;
			default:
				upfurther = false;
				break;
			}
			if (!upfurther) break;
			if (nowLv == 0) return false;
			nowid = upperdivmap[nowLv][nowid];
			nowLv--;
			int jdoff = nowid;
			if (nowLv == 0) { mx = Nx; my = Ny; mz = Nz; }
			else {
				upid = upperdivmap[nowLv][nowid];
				jdoff -= lowerdivmap[nowLv - 1][upid];
			}
			jx = jdoff % mx; jy = jdoff / mx % my; jz = jdoff / mx / my;
			jxs.push(jx); jys.push(jy); jzs.push(jz);
		}
		// Search downward while moving the indices[jx,jy,jz] a step until it reaches to the same level of box
		while (!jxs.empty()) {
			jx = jxs.top(); jy = jys.top(); jz = jzs.top();
			jxs.pop(); jys.pop(); jzs.pop();
			int jxoff = 0, jyoff = 0, jzoff = 0;
			if (nowLv > 0) { mx = nx; my = ny; mz = nz; jxoff = mx; jyoff = my; jzoff = mz; }
			switch (srchdir) {
			case xleft:
				jx = (jx + jxoff - 1) % mx; break;
			case yleft:
				jy = (jy + jyoff - 1) % my; break;
			case zleft:
				jz = (jz + jzoff - 1) % mz; break;
			case xright:
				jx = (jx + 1) % mx; break;
			case yright:
				jy = (jy + 1) % my; break;
			case zright:
				jz = (jz + 1) % mz; break;
			}
			jdoff = jx + mx * (jy + my * jz);
			if (nowLv > 0) jdoff += lowerdivmap[nowLv - 1][upid];
			upid = jdoff;
			nowLv++;
			if (nowLv < divlevel) if (lowerdivmap[nowLv - 1][upid] == lowerdivmap[nowLv - 1][upid + 1]) break;
		}
		NeighLvnId[0] = nowLv - 1; NeighLvnId[1] = jdoff;
	}

void ConnectedDomain::RecursiveConnect(int &ActiveLv, int thisidonset) {
		int mx, my, mz;
		int lowid0 = 0, thisid = thisidonset;
		if (ActiveLv < divlevel - 1) lowid0 = lowerdivmap[ActiveLv][thisid];
		if (ActiveLv == 0) { mx = Nx; my = Ny; mz = Nz;	}
		else { mx = nx; my = ny; mz = nz; }
		for (int iz = 0; iz < mz; iz++) {
			for (int iy = 0; iy < my; iy++) {
				for (int ix = 0; ix < mx; ix++) {
					//if (ActiveLv == 0) cout << "Connecting [ " << ix << ',' << iy << ',' << iz << "] ..." << endl;
					int lowid1 = (ActiveLv < divlevel - 1) ? lowerdivmap[ActiveLv][thisid + 1] : 0;
					if (lowid0 != lowid1) {
						ActiveLv++;
						RecursiveConnect(ActiveLv, lowid0);
					}
					else {
						int inode = serialdivmap[ActiveLv][thisid];
						array<int, 3> ixyz = { ix,iy,iz }; array<int, 2> LvnId = { ActiveLv,thisid };
						array<int, 2> NeighLvnId;
						std::fill(connectInfo[inode].NeighLv.begin(), connectInfo[inode].NeighLv.end(), -1);
						std::fill(connectInfo[inode].NeighId.begin(), connectInfo[inode].NeighId.end(), -1);
						std::fill(connectInfo[inode].NeighInodes.begin(), connectInfo[inode].NeighInodes.end(), -1);
						if (FindNeighbor(xleft, ixyz, LvnId, NeighLvnId)) { 
							connectInfo[inode].NeighLv[0] = NeighLvnId[0];
							connectInfo[inode].NeighId[0] = NeighLvnId[1];
							connectInfo[inode].NeighInodes[0] = serialdivmap[NeighLvnId[0]][NeighLvnId[1]];
						}
						if (FindNeighbor(xright, ixyz, LvnId, NeighLvnId)) {
							connectInfo[inode].NeighLv[1] = NeighLvnId[0];
							connectInfo[inode].NeighId[1] = NeighLvnId[1];
							connectInfo[inode].NeighInodes[1] = serialdivmap[NeighLvnId[0]][NeighLvnId[1]];
						}
						if (FindNeighbor(yleft, ixyz, LvnId, NeighLvnId)) {
							connectInfo[inode].NeighLv[2] = NeighLvnId[0];
							connectInfo[inode].NeighId[2] = NeighLvnId[1];
							connectInfo[inode].NeighInodes[2] = serialdivmap[NeighLvnId[0]][NeighLvnId[1]];
						}
						if (FindNeighbor(yright, ixyz, LvnId, NeighLvnId)) {
							connectInfo[inode].NeighLv[3] = NeighLvnId[0];
							connectInfo[inode].NeighId[3] = NeighLvnId[1];
							connectInfo[inode].NeighInodes[3] = serialdivmap[NeighLvnId[0]][NeighLvnId[1]];
						}
						if (FindNeighbor(zleft, ixyz, LvnId, NeighLvnId)) {
							connectInfo[inode].NeighLv[4] = NeighLvnId[0];
							connectInfo[inode].NeighId[4] = NeighLvnId[1];
							connectInfo[inode].NeighInodes[4] = serialdivmap[NeighLvnId[0]][NeighLvnId[1]];
						}
						if (FindNeighbor(zright, ixyz, LvnId, NeighLvnId)) {
							connectInfo[inode].NeighLv[5] = NeighLvnId[0];
							connectInfo[inode].NeighId[5] = NeighLvnId[1];
							connectInfo[inode].NeighInodes[5] = serialdivmap[NeighLvnId[0]][NeighLvnId[1]];
						}
					}
					lowid0 = lowid1;
					thisid++;
				}
			}
		}
		ActiveLv--;
	}

void ConnectedDomain::CreateConnectInfo() {
		cout << "Enter Create Connect Info" << endl;

		connectInfo.resize(nnode);
		for (int i = 0; i < divlevel; i++) {
			for (int j = 0; j < nnodeLv[i]; j++) {
				int inode = serialdivmap[i][j];
				if (inode > -1) {
					connectInfo[inode].thisLv = i;
					connectInfo[inode].thisId = j;
				}
			}
		}

		int activelv = 0;
		RecursiveConnect(activelv, 0);

		cout << "Connection Done!" << endl;

		//PrintConnectInfo();
	}

void ConnectedDomain::PrintConnectInfo() const {
		ofstream ConctOut("Connection.out");
		for (int i = 0; i < nnode; i++) {
			ConctOut << connectInfo[i].thisLv << ' ' << connectInfo[i].thisId << endl;;
			ConctOut << "  xl : " << connectInfo[i].NeighLv[0] << ' ' << connectInfo[i].NeighId[0] << ' ' << connectInfo[i].NeighInodes[0] << endl;
			ConctOut << "  xr : " << connectInfo[i].NeighLv[1] << ' ' << connectInfo[i].NeighId[1] << ' ' << connectInfo[i].NeighInodes[1] << endl;
			ConctOut << "  yl : " << connectInfo[i].NeighLv[2] << ' ' << connectInfo[i].NeighId[2] << ' ' << connectInfo[i].NeighInodes[2] << endl;
			ConctOut << "  yr : " << connectInfo[i].NeighLv[3] << ' ' << connectInfo[i].NeighId[3] << ' ' << connectInfo[i].NeighInodes[3] << endl;
			ConctOut << "  zl : " << connectInfo[i].NeighLv[4] << ' ' << connectInfo[i].NeighId[4] << ' ' << connectInfo[i].NeighInodes[4] << endl;
			ConctOut << "  zr : " << connectInfo[i].NeighLv[5] << ' ' << connectInfo[i].NeighId[5] << ' ' << connectInfo[i].NeighInodes[5] << endl;
		}

		ofstream ConctOut2("SimpleConnection.out");
		for (int i = 0; i < nnode; i++) {
			for (int j = 0; j < 6; j++) {
				if (connectInfo[i].NeighInodes[j] < 0) continue;
				ConctOut2 << i << ' ' << connectInfo[i].NeighInodes[j] << endl;
			}
		}
	}

void CompiledDomain::Decompile() {
		for (int j = 0; j < compileInfo.size(); j++) {
			for (int i = 0; i < compileInfo[j].inodes.size(); i++) {
				int inode = compileInfo[j].inodes[i];
				double weight = compileInfo[j].weightnodes[i];
				compileid[inode] = j;
				compileW[inode] = weight;
			}
		}
	}

CompiledDomain::CompiledDomain(int3 block, const BaseDomain &rhs) {
		this->block = block;
		
		int3 Nxyz, nxyz;
		int nnode, divlevel;
		const vector<int>& nnodeLv = rhs.GetNnodeLv();
		const vector<vector<int>>& upperdivmap = rhs.GetUpperdivmap();
		const vector<vector<int>>& lowerdivmap = rhs.GetLowerdivmap();
		const vector<vector<int>>& serialdivmap = rhs.GetSerialdivmap();
		const vector<NodeInfo_t>& nodeInfo = rhs.GetBaseNodeInfo();

		rhs.GetBaseSizes(Nxyz, nxyz, nnode, divlevel);
	
		compileid.resize(nnode); compileW.resize(nnode);

		int Nx = Nxyz.x, Ny = Nxyz.y, Nz = Nxyz.z;
		int nx = nxyz.x, ny = nxyz.y, nz = nxyz.z;

		if (Nx % block.x) return;
		if (Ny % block.y) return;
		if (Nz % block.z) return;

		Bx = Nx / block.x, By = Ny / block.y, Bz = Nz / block.z;
		compilemap.resize(Bx * By * Bz + 1);
		std::fill(compilemap.begin(), compilemap.end(), 0);

		nblocks = 0;

		int iblock = 0;
		for (int iz = 0; iz < Bz; iz++) {
		for (int iy = 0; iy < By; iy++) {
		for (int ix = 0; ix < Bx; ix++) {
			vector<int> ivolsInBlock;
			for (int jz = 0; jz < block.z; jz++) {
			for (int jy = 0; jy < block.y; jy++) {
			for (int jx = 0; jx < block.x; jx++) {
				int id0 = (ix * block.x + jx) + Nx * (iy * block.y + jy + Ny * (iz * block.z + jz));
				int id = id0, ActiveLv = 0;

				stack<int> ids;
				ids.push(id0);

				while (ActiveLv > -1) {
					int lowid0 = 0, lowid1 = 0;
					if (ActiveLv < divlevel - 1) lowid0 = lowerdivmap[ActiveLv][id], lowid1 = lowerdivmap[ActiveLv][id + 1];

					if (lowid0 == lowid1) {
						// No lower divisions
						int inode = serialdivmap[ActiveLv][id];
						int idvol = nodeInfo[inode].idvol;
						double volcm3 = nodeInfo[inode].volcm3;

						int idInBlock = 0; bool isin = false;
						for (auto iter = ivolsInBlock.begin(); iter != ivolsInBlock.end(); iter++) {
							if (idvol == *iter) { isin = true; break; }
							idInBlock++;
						}
						if (!isin) { ivolsInBlock.push_back(idvol); compileInfo.emplace_back(); }

						compileInfo[nblocks + idInBlock].idvol = idvol;
						compileInfo[nblocks + idInBlock].inodes.push_back(inode);
						compileInfo[nblocks + idInBlock].weightnodes.push_back(volcm3);
						compileInfo[nblocks + idInBlock].volsum += volcm3;

						// Check if this level of divisions are all sweeped
						if (ActiveLv > 0) {
							int upid = ids.top();
							int uplowid1 = lowerdivmap[ActiveLv - 1][upid + 1];
							// Move to next id (id++), and compare the end of points in this division
							if (uplowid1 == ++id) {
								// If it was the end, go upward and move to upid++
								ids.pop(); id = ++upid;
								ActiveLv--;
							}
						}
						if (ActiveLv == 0) break;
					}
					// Go downward further
					else { ids.push(id);	ActiveLv++;	id = lowid0; }
				}
			}	}	}
			nblocks += ivolsInBlock.size();
			compilemap[++iblock] = nblocks;
		} } }
		// weightnodes calculation
		for (int i = 0; i < compileInfo.size(); i++) {
			for (int j = 0; j < compileInfo[i].weightnodes.size(); j++) {
				compileInfo[i].weightnodes[j] /= compileInfo[i].volsum;
			}
		}
		// Decompile
		Decompile();
		//PrintCompileInfo();
	}

CompiledDomain::CompiledDomain(int3 block, const CompiledDomain &rhs) {
		this->block = block;
		
		int Bxr, Byr, Bzr;
		int3 blockr;
		int nblocksr;
		rhs.GetCompileSizes(Bxr, Byr, Bzr, blockr, nblocksr);

		compileid.resize(nblocksr); compileW.resize(nblocksr);

		if (block.x % blockr.x) return;
		if (block.y % blockr.y) return;
		if (block.z % blockr.z) return;
		int divx = block.x / blockr.x, divy = block.y / blockr.y, divz = block.z / blockr.z;
		Bx = Bxr / divx; By = Byr / divy; Bz = Bzr / divz;
		compilemap.resize(Bx * By * Bz + 1);

		const vector<int> &compilemapr = rhs.GetCompileMap();
		const vector<CompileInfo_t> &compileInfor = rhs.GetCompileInfo();

		nblocks = 0;

		int iblock = 0;
		for (int iz = 0; iz < Bz; iz++) {
		for (int iy = 0; iy < By; iy++) {
		for (int ix = 0; ix < Bx; ix++) {
			vector<int> ivolsInBlock;
			for (int jz = 0; jz < block.z; jz++) {
			for (int jy = 0; jy < block.y; jy++) {
			for (int jx = 0; jx < block.x; jx++) {
				int id0 = (jx + ix * block.x) + Bxr * (jy + iy * block.y + Byr * (jz + iz * block.z));

				for (int inlow = compilemapr[id0]; inlow < compilemapr[id0 + 1]; inlow++) {
					double volcm3 = compileInfor[inlow].volsum;
					int idvol = compileInfor[inlow].idvol;

					int idInBlock = 0; bool isin = false;
					for (auto iter = ivolsInBlock.begin(); iter != ivolsInBlock.end(); iter++) {
						if (idvol == *iter) { isin = true; break; }
						idInBlock++;
					}
					if (!isin) { ivolsInBlock.push_back(idvol); compileInfo.emplace_back(); }

					compileInfo[nblocks + idInBlock].idvol = idvol;
					compileInfo[nblocks + idInBlock].inodes.push_back(inlow);
					compileInfo[nblocks + idInBlock].weightnodes.push_back(volcm3);
					compileInfo[nblocks + idInBlock].volsum += volcm3;
				}
			}	}	}
			nblocks += ivolsInBlock.size();
			compilemap[++iblock] = nblocks;
		} }	}
		// weightnodes calculation
		for (int i = 0; i < compileInfo.size(); i++) {
			for (int j = 0; j < compileInfo[i].weightnodes.size(); j++) {
				compileInfo[i].weightnodes[j] /= compileInfo[i].volsum;
			}
		}
		// Decompile
		Decompile();
		//PrintCompileInfo("Compile2.out");
	}

void CompiledDomain::PrintCompileInfo(string filename) const {
		ofstream compile(filename);
		for (int i = 0; i < compilemap.size() - 1; i++) {
			compile << "Box : " << i << endl;
			for (int j = compilemap[i]; j < compilemap[i + 1]; j++) {
				compile << '\t' << "Compiled block : " << j;
				compile << " Vol. ID : " << compileInfo[j].idvol;
				compile << " Volume [cm3] : " << compileInfo[j].volsum << endl;
				for (int k = 0; k < compileInfo[j].inodes.size(); k++) {
					compile << '\t' << '\t' << "Node ID : " << compileInfo[j].inodes[k];
					compile << " Weight : " << compileInfo[j].weightnodes[k] << endl;
				}
			}
		}
		cout << "Compile Info Printed Out!" << endl;
	}

void RaySegmentDomain::Initialize(int ng, int nangle_oct) { 
	bndflux.Create(ng, 2 * nangle_oct, 2 * (Nx*Ny + Ny*Nz + Nx*Nz)); 
}

void FlatSrcDomain::Initialize(int ng, int scatorder, bool isEx) {
	flux.Create(ng, nblocks, (scatorder + 1) * (scatorder + 1));
	srcS(ng, nblocks);
	srcF.Create(ng, nblocks);
	if (isEx) srcEx.Create(ng, nblocks);
}

void FlatXSDomain::InitializeMacro(int ng, int scatorder, bool isTHfeed) {
	this->ng = ng; this->scatorder = scatorder;
	isMicro = false; this->isTHfeed = isTHfeed;
	if (isTHfeed) temperature.Create(nblocks);
	xst.Create(ng, nblocks); xssm.Create(ng, ng, nblocks, scatorder+1); xsnf.Create(ng, nblocks);
	if (isTHfeed) xskf.Create(ng, nblocks);
}

void FlatXSDomain::Initialize(int ng, int scatorder, int niso) {
	this->ng = ng; this->scatorder = scatorder;
	isMicro = true; isTHfeed = true;
	idiso.Create(niso, nblocks); pnum.Create(niso, nblocks); temperature.Create(nblocks);
	xst.Create(ng, nblocks); xssm.Create(ng, ng, scatorder+1, nblocks); xsnf.Create(ng, nblocks);
	xskf.Create(ng, nblocks);
}

void FlatXSDomain::SetMacroXS(const vector<int> &imag, const XSLibrary &XS)
{
	const auto& Micro = XS.GetMicroXS();

	for (int i = 0; i < compileInfo.size(); ++i) {
		int mat = compileInfo[i].idvol;

		const auto& total = Micro[mat].total;
		const auto& scat = Micro[mat].scat;
		const auto& fis = Micro[mat].fis, &nu = Micro[mat].nu, &kappa = Micro[mat].kappa;

		std::copy(total.begin(), total.end(), &xst(0, i));
		std::copy(scat.begin(), scat.end(), &xssm(0, 0, 0, i));
		
		for (int g = 0; g < ng; ++g) {
			xsnf(g, i) = fis(g) * nu(g);
			xskf(g, i) = fis(g) * kappa(g);
		}

	}

}

}