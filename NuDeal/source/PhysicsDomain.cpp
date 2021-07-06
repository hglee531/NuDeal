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
		lx0 = Width0.x; ly0 = Width0.y; lz0 = Width0.z;
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

void RayTracingDomain::RecursiveTrace(array<queue<int2>, 8> &Qsweeploop, int Lv, int id, int &iternode, int Oct) {
	for (int iz = range[Oct][0].z; iz != range[Oct][1].z; iz += dxyz[Oct].z) {
		for (int iy = range[Oct][0].y; iy != range[Oct][1].y; iy += dxyz[Oct].y) {
			for (int ix = range[Oct][0].x; ix != range[Oct][1].x; ix += dxyz[Oct].x) {
				int thisid = id + ix + (iy + iz * ny) * nx;
				int lowid0 = 0, lowid1 = 0;
				if (Lv < divlevel - 1) lowid0 = lowerdivmap[Lv][thisid], lowid1 = lowerdivmap[Lv][thisid + 1];

				if (lowid0 < lowid1) {
					Qsweeploop[Oct].push(make_int2(iternode, 0));
					nloop[Oct]++;
					RecursiveTrace(Qsweeploop, Lv + 1, lowid0 + edgeid[Oct], iternode, Oct);
				}
				else {
					sweepnodes(iternode, Oct) = serialdivmap[Lv][thisid];
					sweepIdInBox(iternode, Oct) = thisid;
					iternode++;
				}
			} 
			Qsweeploop[Oct].push(make_int2(iternode, 1));
			nloop[Oct]++;
		}
		Qsweeploop[Oct].back().y = 2;
	}
	Qsweeploop[Oct].back().y = 3;
}

void RayTracingDomain::SetTraceOrder(array<queue<int2>, 8> &Qsweeploop) {
	std::fill(nloop, nloop + 8, 0);

	for (int Oct = 0; Oct < 8; Oct++) {
		int iternode = 0;
		Qsweeploop[Oct].push(make_int2(0, -1));
		for (int iz = range0[Oct][0].z; iz != range0[Oct][1].z; iz += dxyz[Oct].z) {
			for (int iy = range0[Oct][0].y; iy != range0[Oct][1].y; iy += dxyz[Oct].y) {
				for (int ix = range0[Oct][0].x; ix != range0[Oct][1].x; ix += dxyz[Oct].x) {
					int thisLv = 0, thisid = ix + (iy + iz * Ny) * Nx;
					int lowid0 = 0, lowid1 = 0;
					if (divlevel > 1) lowid0 = lowerdivmap[thisLv][thisid], lowid1 = lowerdivmap[thisLv][thisid + 1];

					if (lowid0 < lowid1) {
						Qsweeploop[Oct].push(make_int2(iternode, 0));
						nloop[Oct]++; 
						RecursiveTrace(Qsweeploop, thisLv + 1, lowid0, iternode, Oct);
					}
					else {
						sweepnodes(iternode, Oct) = serialdivmap[thisLv][thisid];
						sweepIdInBox(iternode, Oct) = thisid;
						iternode++;
					}
				}
				Qsweeploop[Oct].push(make_int2(iternode, 1));
				nloop[Oct]++;
			}
			Qsweeploop[Oct].back().y = 2;
		}
		Qsweeploop[Oct].back().y = 3;
	}
}

void RayTracingDomain::CreateTracingDomain() {
	array<queue<int2>, 8> Qsweeploop;

	sweepnodes.Create(nnode, 8);
	sweepIdInBox.Create(nnode, 8);

	for (int Oct = 0; Oct < 8; Oct++) {
		int dx = (Oct % 2) ? -1 : 1, dy = (Oct % 4 / 2) ? -1 : 1, dz = (Oct / 4) ? -1 : 1;
		dxyz[Oct] = make_int3(dx, dy, dz);

		int xs0[2], ys0[2], zs0[2], xs[2], ys[2], zs[2];

		if (dx > 0) { xs0[0] = 0; xs0[1] = Nx; xs[0] = 0; xs[1] = nx; }
		else { xs0[0] = Nx - 1; xs0[1] = -1; xs[0] = nx - 1; xs[1] = -1; }
		if (dy > 0) { ys0[0] = 0; ys0[1] = Ny; ys[0] = 0; ys[1] = ny; }
		else { ys0[0] = Ny - 1; ys0[1] = -1; ys[0] = ny - 1; ys[1] = -1; }
		if (dz > 0) { zs0[0] = 0; zs0[1] = Nz; zs[0] = 0; zs[1] = nz; }
		else { zs0[0] = Nz - 1; zs0[1] = -1; zs[0] = nz - 1; zs[1] = -1; }

		edgeid[Oct] = xs[0] + (ys[0] + zs[0] * Ny) * Nx;

		range0[Oct][0] = make_int3(xs0[0], ys0[0], zs0[0]); range[Oct][0] = make_int3(xs[0], ys[0], zs[0]);
		range0[Oct][1] = make_int3(xs0[1], ys0[1], zs0[1]); range[Oct][1] = make_int3(xs[1], ys[1], zs[1]);
	}

	SetTraceOrder(Qsweeploop);

	sweeploop.Create(Qsweeploop[0].size(), 8);
	for (int Oct = 0; Oct < 8; Oct++) {
		int size = Qsweeploop[Oct].size();
		for (int i = 0; i < size; i++) {
			sweeploop(i, Oct) = Qsweeploop[Oct].front();
			Qsweeploop[Oct].pop();
		}
	}

	PrintTracingInfo();
}

void RayTracingDomain::PrintTracingInfo() {
	ofstream traceout("Tracing.out");
	for (int Oct = 0; Oct < 8; Oct++) {
		traceout << "Octant : " << Oct << endl;

		int size = sweeploop.size() / 8, thisLv = 0;
		for (int iloop = 0; iloop < size - 1; iloop++) {
			int ibeg = sweeploop(iloop, Oct).x;
			int iend = sweeploop(iloop + 1, Oct).x, storedir = sweeploop(iloop + 1, Oct).y;

			if (storedir == 0) {
				for (int i = 0; i < thisLv; i++) traceout << '\t';
				thisLv++;
				traceout << "Enter(Lv) : " << thisLv << endl;
			}

			if (ibeg < iend)
				for (int i = 0; i < thisLv; i++) traceout << '\t';

			for (int inode = ibeg; inode < iend; inode++) {
				traceout << sweepnodes(inode, Oct) << ' ';
			}

			if (ibeg < iend) traceout << endl;

			switch (storedir)
			{
			case 1:
				for (int i = 0; i < thisLv; i++) traceout << '\t';
				traceout << "Record fluxes in X..." << endl;
				break;
			case 2:
				for (int i = 0; i < thisLv; i++) traceout << '\t';
				traceout << "Record fluxes in y..." << endl;
				break;
			case 3:
				for (int i = 0; i < thisLv; i++) traceout << '\t';
				traceout << "Record fluxes in z..." << endl;
				thisLv--;
				for (int i = 0; i < thisLv; i++) traceout << '\t';
				traceout << "Exit(Lv) : " << thisLv << endl;
				break;
			default:
				break;
			}
		}

		traceout << endl;
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
			for (int jz = 0; jz < divz; jz++) {
			for (int jy = 0; jy < divy; jy++) {
			for (int jx = 0; jx < divx; jx++) {
				int id0 = (jx + ix * divx) + Bxr * (jy + iy * divy + Byr * (jz + iz * divz));

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

void FluxBoundary::Initialize(int ng, int nangle_oct) {
	this->ng = ng; this->nangle_oct = nangle_oct;
	xbndflux.Create(ng, nangle_oct, Ny*Nz, 8);
	ybndflux.Create(ng, nangle_oct, Nx*Nz, 8);
	zbndflux.Create(ng, nangle_oct, Nx*Ny, 8);
}

void  FluxBoundary::VaccumFlux() {
	std::fill(xbndflux.begin(), xbndflux.end(), 0.0);
	std::fill(ybndflux.begin(), ybndflux.end(), 0.0);
	std::fill(zbndflux.begin(), zbndflux.end(), 0.0);
}

void  FluxBoundary::InitBndFluxUnity() {
	std::fill(xbndflux.begin(), xbndflux.end(), 1.0);
	std::fill(ybndflux.begin(), ybndflux.end(), 1.0);
	std::fill(zbndflux.begin(), zbndflux.end(), 1.0);
}

void FluxBoundary::UpdateBoundary(int leftcond, int rightcond) {
	if (leftcond + rightcond == 0) {
		VaccumFlux();
		return;
	}

	bool isxl, isxr, isyl, isyr, iszl, iszr;
	isxl = leftcond / 100; isxr = rightcond / 100;
	isyl = leftcond / 10 % 10; isyr = rightcond / 10 % 10;
	iszl = leftcond % 10; iszr = rightcond % 10;
		
	Array<double> bufx, bufy, bufz;
	bufx = xbndflux; bufy = ybndflux; bufz = zbndflux;
	
	if (isxl) {
		int octl[4] = { 0, 2, 4, 6 };
		int octr[4] = { 1, 3, 5, 7 };
		for (int Oct = 0; Oct < 4; Oct++) {
			int octout = octr[Oct], octin = octl[Oct];
			std::copy(&bufx(0, 0, 0, octout), &bufx(0, 0, 0, octout + 1), &xbndflux(0, 0, 0, octin));
		}
	}
	if (isxr) {
		int octr[4] = { 0, 2, 4, 6 };
		int octl[4] = { 1, 3, 5, 7 };
		for (int Oct = 0; Oct < 4; Oct++) {
			int octout = octr[Oct], octin = octl[Oct];
			std::copy(&bufx(0, 0, 0, octout), &bufx(0, 0, 0, octout + 1), &xbndflux(0, 0, 0, octin));
		}
	}

	if (isyl) {
		int octl[4] = { 0, 1, 4, 5 };
		int octr[4] = { 2, 3, 6, 7 };
		for (int Oct = 0; Oct < 4; Oct++) {
			int octout = octr[Oct], octin = octl[Oct];
			std::copy(&bufy(0, 0, 0, octout), &bufy(0, 0, 0, octout + 1), &ybndflux(0, 0, 0, octin));
		}
	}
	if (isyr) {
		int octr[4] = { 0, 1, 4, 5 };
		int octl[4] = { 2, 3, 6, 7 };
		for (int Oct = 0; Oct < 4; Oct++) {
			int octout = octr[Oct], octin = octl[Oct];
			std::copy(&bufy(0, 0, 0, octout), &bufy(0, 0, 0, octout + 1), &ybndflux(0, 0, 0, octin));
		}
	}

	if (iszl) {
		int octl[4] = { 0, 1, 2, 3 };
		int octr[4] = { 4, 5, 6, 7 };
		for (int Oct = 0; Oct < 4; Oct++) {
			int octout = octr[Oct], octin = octl[Oct];
			std::copy(&bufz(0, 0, 0, octout), &bufz(0, 0, 0, octout + 1), &zbndflux(0, 0, 0, octin));
		}
	}
	if (iszr) {
		int octr[4] = { 0, 1, 2, 3 };
		int octl[4] = { 4, 5, 6, 7 };
		for (int Oct = 0; Oct < 4; Oct++) {
			int octout = octr[Oct], octin = octl[Oct];
			std::copy(&bufz(0, 0, 0, octout), &bufz(0, 0, 0, octout + 1), &zbndflux(0, 0, 0, octin));
		}
	}
	
}

void FlatSrcDomain::Initialize(int ng, int nangle_octant, int scatorder, bool isEx) {
	this->ng = ng; this->nangle_oct = nangle_octant;
	this->scatorder = scatorder; this->isEx = isEx;
	flux.Create(ng, nblocks, (scatorder + 1) * (scatorder + 1));
	psi.Create(nblocks);
	src.Create(ng, nangle_octant, nblocks, 8);
	srcS.Create(ng, nangle_octant, nblocks, 8);
	if (isEx) srcEx.Create(ng, nangle_octant, nblocks, 8);
}

double FlatSrcDomain::UpdateFisSource(const FlatXSDomain &FXR) {
	using Math::square;

	const vector<int> &idFXR = FXR.GetDecompileId();
	const Array<double> &xsnf = FXR.GetNuFisXS();

	double norm0, norm1;
	norm0 = norm1 = 0;

	for (int i = 0; i < nblocks; i++) {
		double voli = compileInfo[i].volsum;

		//double psi0 = psi(i);
		double psi0 = psi(i) * voli;

		psi(i) = 0.0;

		int iFXR = idFXR[i];
		for (int ig = 0; ig < ng; ig++)
			psi(i) += flux(ig, i) * xsnf(ig, iFXR);
		
		//norm0 += square(psi(i));
		//norm1 += psi(i) * psi0;
		norm0 += square(psi(i) * voli);
		norm1 += psi(i) * voli * psi0;
	}

	return norm0 / norm1;
}

void FlatSrcDomain::NormalizePsi(double keff) {
	for (auto it = psi.begin(); it < psi.end(); it++) (*it) /= keff;
}

void FlatSrcDomain::AccumulateSource(const vector<int> &idFXR, const Array<double> &chi) {
	using Math::PI;
	for (int i = 0; i < nblocks; i++) {
		int iFXR = idFXR[i];
		for (int ig = 0; ig < ng; ig++) {
			double sf = psi(i) * chi(ig, iFXR) / PI * 0.25;
			for (int Oct = 0; Oct < 8; Oct++) {
				for (int iang = 0; iang < nangle_oct; iang++) {
					src(ig, iang, i, Oct) = sf + srcS(ig, iang, i, Oct);
				}
			}
		}
	}
	if (isEx)
		for (int i = 0; i < src.size(); i++) src(i) += srcEx(i);
}

void FlatSrcDomain::UpdateScatSource(const FlatXSDomain &FXR) {
	using Math::PI;
	const vector<int> &idFXR = FXR.GetDecompileId();
	const Array<double> &xssm = FXR.GetScatMatrix();

	std::fill(srcS.begin(), srcS.end(), 0.0);

	for (int i = 0; i < nblocks; i++) {
		int iFXR = idFXR[i];
		for (int ig = 0; ig < ng; ig++)
			for (int igg = 0; igg < ng; igg++) {
				double ss = flux(ig, i) * xssm(igg, ig, 0, iFXR) / PI * 0.25;

				for (int Oct = 0; Oct < 8; Oct++)
					for (int iangle = 0; iangle < nangle_oct; iangle++)
						srcS(igg, iangle, i, Oct) += ss;
			}
	}
}

void FlatSrcDomain::UpdateAngularSource(const FlatXSDomain &FXR) {
	using Math::PI;
	const vector<int> &idFXR = FXR.GetDecompileId();
	const Array<double> &xssm = FXR.GetScatMatrix();

	std::fill(srcS.begin(), srcS.end(), 0.0);

	for (int i = 0; i < nblocks; i++) {
		int iFXR = idFXR[i];

		for (int ig = 0; ig < ng; ig++)
			for (int igg = 0; igg < ng; igg++) {
				double ss = flux(ig, i) * xssm(igg, ig, 0, iFXR) / PI * 0.25;

				for (int Oct = 0; Oct < 8; Oct++)
					for (int iangle = 0; iangle < nangle_oct; iangle++)
						srcS(igg, iangle, i, Oct) += ss;
			}
	}

	const Array<double> &chi = FXR.GetFisSpectrum();
	AccumulateSource(idFXR, chi);
}

void FlatXSDomain::InitializeMacro(int ng, int scatorder, bool isTHfeed) {
	this->ng = ng; this->scatorder = scatorder;
	isMicro = false; this->isTHfeed = isTHfeed;
	if (isTHfeed) temperature.Create(nblocks);
	xst.Create(ng, nblocks); xssm.Create(ng, ng, scatorder+1, nblocks); 
	xsnf.Create(ng, nblocks); chi.Create(ng, nblocks);
	if (isTHfeed) xskf.Create(ng, nblocks);
}

void FlatXSDomain::Initialize(int ng, int scatorder, int niso) {
	this->ng = ng; this->scatorder = scatorder;
	isMicro = true; isTHfeed = true;
	idiso.Create(niso, nblocks); pnum.Create(niso, nblocks); temperature.Create(nblocks);
	xst.Create(ng, nblocks); xssm.Create(ng, ng, scatorder+1, nblocks); xsnf.Create(ng, nblocks);
	xskf.Create(ng, nblocks);
}

void FlatXSDomain::SetMacroXS(const XSLibrary &XS, const vector<int> imat)
{
	const auto& Macro = XS.GetMacroXS();

	for (int i = 0; i < compileInfo.size(); ++i) {
		int mat = imat[compileInfo[i].idvol];

		const auto& total = Macro[mat].tr;
		const auto& scat = Macro[mat].scat;
		const auto& nufis = Macro[mat].nufis, &kappafis = Macro[mat].kappafis;
		const auto& chi = Macro[mat].chi;

		std::copy(total.begin(), total.end(), &xst(0, i));
		std::copy(scat.begin(), scat.end(), &xssm(0, 0, 0, i));
		std::copy(nufis.begin(), nufis.end(), &xsnf(0, i));
		std::copy(chi.begin(), chi.end(), &this->chi(0, i));
		if (isTHfeed) std::copy(kappafis.begin(), kappafis.end(), &xskf(0, i));
	}

}

}