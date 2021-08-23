#include "UnitGeo.h"

namespace Geometry
{

const int pow10[4] = { 1, 10, 100, 1000 };

inline 
UnitSurf::Equation_t UnitSurf::CalCoeffs(const Equation_t& in, double xyzc[4][4])
{
	// xyzc : (x', y', z', 1) = g(x,y,z)
	// Convert f(x,y,z) to f(x',y',z') when (x',y',z') = g(x,y,z)

	Equation_t out;

	double c160[4][4], c161[4][4] = { { 0.0, }, };
	c160[0][0] = in[0]; c160[1][1] = in[1]; c160[2][2] = in[2]; c160[3][3] = in[9];
	c160[0][1] = in[3]; c160[0][3] = in[6];
	c160[1][2] = in[4]; c160[1][3] = in[7];
	c160[0][2] = in[5]; c160[2][3] = in[8];
	for (int i = 0; i < 4; i++) {
		for (int j = i; j < 4; j++) {
			for (int ir = 0; ir < 4; ir++) {
				for (int ic = 0; ic < 4; ic++) {
					c161[ir][ic] += c160[i][j] * xyzc[i][ir] * xyzc[j][ic];
				}
			}
		}
	}
	out[0] = c161[0][0]; out[1] = c161[1][1]; out[2] = c161[2][2]; out[9] = c161[3][3];
	out[3] = c161[0][1] + c161[1][0]; out[6] = c161[0][3] + c161[3][0];
	out[4] = c161[0][2] + c161[2][0]; out[7] = c161[1][3] + c161[3][1];
	out[5] = c161[1][2] + c161[2][1]; out[8] = c161[2][3] + c161[3][2];

	return static_cast<Equation_t>(out);
}

inline 
int UnitSurf::GetLocalExCurve(int axis, double coeff[6], double sol[2][2]) {
	// Based on Lagrange multiplier method
	// G(x) = x & F(x,y) = 0;
	// grad(F)=lambda*grad(G)

	// return single digit value
	// 0 : none for along the axis
	// 1 : left only, 2 : right only, 3 : both sides
	double gradG[2][3] = {
		{2.0*coeff[0], coeff[2], coeff[3]},
		{coeff[2], 2.0*coeff[1], coeff[4]}
	};
	int i0 = axis, i1 = (axis + 1) % 2;
	if (abs(gradG[i1][i1]) < eps_geo) return 0;
	double c1 = -gradG[i1][i0] / gradG[i1][i1], c0 = -gradG[i1][2] / gradG[i1][i1];
	double a, b, c;
	a = coeff[i0] + coeff[2] * c1 + coeff[i1] * c1 * c1; // [i0][i0] + [i0][i1] + [i1][i1]
	b = coeff[3 + i0] + coeff[2] * c0 + 2.0 * coeff[i1] * c0 * c1 + coeff[2] * c1; // [i0] + [i0][i1] + [i1]*[i1] + [i1]
	c = coeff[i1] * c0 * c0 + coeff[3 + i1] * c0 + coeff[5]; // [i1][i1] + [i1] + constant
	if (abs(a) < eps_geo) return 0;
	else {
		b /= a; c /= a;
		double det = b * b - 4.0 * c;
		if (abs(det) < eps_geo) {
			double soli0 = -b * 0.5;
			double soli1 = c1 * soli0 + c0;
			sol[0][i0] = sol[1][i0] = soli0; sol[0][i1] = sol[1][i1] = soli1;
			double lambda = gradG[i0][i0] * soli0 + gradG[i0][i1] * soli1 + gradG[i0][2];
			if (lambda < 0) return 1;
			else return 2;
		}
		else if (det > 0) {
			double soli0 = 0.5*(-b - sqrt(det));
			double soli1 = c1 * soli0 + c0;
			sol[0][i0] = soli0; sol[0][i1] = soli1;
			soli0 = 0.5*(-b + sqrt(det));
			soli1 = c1 * soli0 + c0;
			sol[1][i0] = soli0; sol[1][i1] = soli1;
			return 3;
		}
	}
}

inline 
bool UnitSurf::GetTriplePoint(double coeff0[4], double coeff1[4], double coeff2[4], double sol[3]) 
{
	double det12yz, det12xz, det12xy;
	det12yz = coeff1[1] * coeff2[2] - coeff1[2] * coeff2[1];
	det12xz = coeff1[0] * coeff2[2] - coeff1[2] * coeff2[0];
	det12xy = coeff1[0] * coeff2[1] - coeff1[1] * coeff2[0];
	double det = coeff0[0] * det12yz - coeff0[1] * det12xz + coeff0[2] * det12xy;
	if (abs(det) < eps_geo) return false;
	double det02yz, det02xz, det02xy;
	det02yz = coeff0[1] * coeff2[2] - coeff0[2] * coeff2[1];
	det02xz = coeff0[0] * coeff2[2] - coeff0[2] * coeff2[0];
	det02xy = coeff0[0] * coeff2[1] - coeff0[1] * coeff2[0];
	double det01yz, det01xz, det01xy;
	det01yz = coeff0[1] * coeff1[2] - coeff0[2] * coeff1[1];
	det01xz = coeff0[0] * coeff1[2] - coeff0[2] * coeff1[0];
	det01xy = coeff0[0] * coeff1[1] - coeff0[1] * coeff1[0];
	sol[0] = -det12yz * coeff0[3] + det02yz * coeff1[3] - det01yz * coeff2[3];
	sol[1] = det12xz * coeff0[3] - det02xz * coeff1[3] + det01xz * coeff2[3];
	sol[2] = -det12xy * coeff0[3] + det02xy * coeff1[3] - det01xy * coeff2[3];
	sol[0] /= det; sol[1] /= det;	sol[2] /= det;
	return true;

}


void UnitSurf::Transform(double invTM[4][4]) {
	eq = CalCoeffs(eq, invTM);
}

void UnitSurf::SetCurve()
{
	is_curve = false;
	for (int i = 0; i < 6; ++i) {
		if (abs(eq[i]) > eps_geo) { is_curve = true; break; }
	}
}

bool UnitSurf::GetPointVar2(CartAxis axis, double var[2], double &sol) {
	int i0 = static_cast<int>(axis), i1, i2;
	auto& coeff = eq;
	switch (axis) {
	case CartAxis::X:
		i1 = 1; i2 = 2; break;
	case CartAxis::Y:
		i1 = 2; i2 = 0; break;
	case CartAxis::Z:
		i1 = 0; i2 = 1; break;
	}
	double a, b, c;
	a = coeff[i0]; b = coeff[3 + i0] * var[i1] + coeff[3 + i2] * var[i2] + coeff[6 + i0];
	c = coeff[i1] * var[i1] * var[i1] + coeff[i2] * var[i2] * var[i2] + coeff[3 + i1] * var[i1] * var[i2] +
		coeff[6 + i1] * var[i1] + coeff[6 + i2] * var[i1] + coeff[9];
	if (abs(a) < eps_geo) {
		if (abs(b) > eps_geo) {
			sol = - b / c;
			return true;
		}
		else {
			return false;
		}
	}
	else {
		sol = - 0.5 * b / a;
		return true;
	}
}

void UnitSurf::Create(SurfType surftype, double *_coeff, CartPlane cartesianPln) {
	auto& coeff = eq;
	double sign;
	switch (surftype) {
	case SurfType::XPLN:
		sign = -1.;
		if (_coeff[0] < 0) sign = 1.;
		coeff[6] = _coeff[0];
		coeff[9] = sign*_coeff[1];
		break;
	case SurfType::YPLN:
		sign = -1.;
		if (_coeff[0] < 0) sign = 1.;
		coeff[7] = _coeff[0];
		coeff[9] = sign*_coeff[1];
		break;
	case SurfType::ZPLN:
		sign = -1.;
		if (_coeff[0] < 0) sign = 1.;
		coeff[8] = _coeff[0];
		coeff[9] = sign*_coeff[1];
		break;
	case SurfType::PLN:
		coeff[6] = _coeff[0];
		coeff[7] = _coeff[1];
		coeff[8] = _coeff[2];
		coeff[9] = _coeff[3];
		break;
	case SurfType::CIRCLE:
		switch (cartesianPln) {
		case CartPlane::XY:
			sign = 1.;
			if (_coeff[2] < 0) sign = -1.;
			coeff[0] = sign; coeff[1] = sign;
			coeff[6] = - sign * 2. * _coeff[0]; coeff[7] = - sign * 2.*_coeff[1];
			coeff[9] = sign * (sign*_coeff[0] * _coeff[0] + _coeff[1] * _coeff[1] - _coeff[2] * _coeff[2]);
			break;
		case CartPlane::YZ:
			sign = 1.;
			if (_coeff[2] < 0) sign = -1.;
			coeff[1] = sign; coeff[2] = sign;
			coeff[7] = - sign * 2. * _coeff[0]; coeff[8] =-  sign * 2.*_coeff[1];
			coeff[9] = sign * (sign*_coeff[0] * _coeff[0] + _coeff[1] * _coeff[1] - _coeff[2] * _coeff[2]);
			break;
		case CartPlane::XZ:
			sign = 1.;
			if (_coeff[2] < 0) sign = -1.;
			coeff[0] = sign; coeff[2] = sign;
			coeff[6] = - sign * 2. * _coeff[0]; coeff[8] = - sign * 2.*_coeff[1];
			coeff[9] = sign * (sign*_coeff[0] * _coeff[0] + _coeff[1] * _coeff[1] - _coeff[2] * _coeff[2]);
			break;
		}
		break;
	default: // GENERAL
		for (int i = 0; i < 10; i++) coeff[i] = _coeff[i];
		break;
	}

	SetCurve();
}

void UnitSurf::Relocate(double dx, double dy, double dz) {
	double invTM[4][4] = { { 0.0, }, };
	invTM[0][0] = 1.0;
	invTM[1][1] = 1.0;
	invTM[2][2] = 1.0;
	invTM[3][3] = 1.0;
	// Update only the fourth rows
	invTM[0][3] = -dx; invTM[1][3] = -dy; invTM[2][3] = -dz;
	Transform(invTM);
}

void UnitSurf::Rotate(double cos, double sin, CartAxis Ax) {
	double invTM[4][4] = { { 0.0, }, };
	invTM[0][0] = 1.0;
	invTM[1][1] = 1.0;
	invTM[2][2] = 1.0;
	invTM[3][3] = 1.0;

	int ir1, ir2;
	switch (Ax) {
	case CartAxis::X:
		ir1 = 1; ir2 = 2;
		break;
	case CartAxis::Y:
		ir1 = 2; ir2 = 0;
		break;
	case CartAxis::Z:
		ir1 = 0; ir2 = 1;
		break;
	default:
		break;
	}
	invTM[ir1][ir1] = cos; invTM[ir2][ir2] = cos;
	invTM[ir1][ir2] = sin; invTM[ir2][ir1] = -sin;
	Transform(invTM);
}

int UnitSurf::GetIntersection(CartPlane CartPlane, const array<double, 2>& val, array<double, 2>& sol) {
	double c2, c1, c0;
	if (is_curve) {
		switch (CartPlane) {
		case CartPlane::XY:
			c2 = c_zs;
			c1 = c_xz * val[0] + c_yz * val[1] + c_z;
			c0 = c_xs * val[0] * val[0] + c_ys * val[1] * val[1] + c_xy * val[0] * val[1] +
				c_x * val[0] + c_y * val[1] + c;
			break;
		case CartPlane::YZ:
			c2 = c_xs;
			c1 = c_xy * val[0] + c_xz * val[1] + c_x;
			c0 = c_ys * val[0] * val[0] + c_zs * val[1] * val[1] + c_yz * val[0] * val[1] +
				c_y * val[0] + c_z * val[1] + c;
			break;
		case CartPlane::XZ:
			c2 = c_ys;
			c1 = c_xy * val[0] + c_yz * val[1] + c_y;
			c0 = c_xs * val[0] * val[0] + c_zs * val[1] * val[1] + c_xz * val[0] * val[1] +
				c_x * val[0] + c_z * val[1] + c;
			break;
		default:
			break;
		}
		if (abs(c2) < eps_geo) return 0;
		double det = c1 * c1 - 4.*c2*c0;
		if (det < 0.) return 0;
		else {
			double dist = sqrt(det) / c2;
			sol[0] = (-c1 / c2 - dist) * 0.5;
			sol[1] = (-c1 / c2 + dist) * 0.5;
			if (dist < 1.e-12) {
				return 0;
			}
			else {
				return 2;
			}
		}
	}
	else {
		int rtrnval = 1;
		switch (CartPlane) {
		case CartPlane::XY:
			if (abs(c_z) < eps_geo) rtrnval = 0;
			sol[0] = -(c_x*val[0] + c_y * val[1] + c) / c_z;
			break;
		case CartPlane::YZ:
			if (abs(c_x) < eps_geo) rtrnval = 0;
			sol[0] = -(c_y*val[0] + c_z * val[1] + c) / c_x;
			break;
		case CartPlane::XZ:
			if (abs(c_y) < eps_geo) rtrnval = 0;
			sol[0] = -(c_x*val[0] + c_z * val[1] + c) / c_y;
			break;
		default:
			break;
		}
		return rtrnval;
	}
}

bool UnitSurf::IsInside(double x, double y, double z, bool includeOn) {
	double det = (c_xs*x + c_xy * y + c_xz * z + c_x)*x + (c_ys*y + c_yz * z + c_y)*y + (c_zs*z + c_z)*z + c;
	if (abs(det) < eps_geo) {
		if (includeOn) return true;
		return false;
	}
	return (det < 0);
}

int UnitSurf::GetLocalExSelf(double sol[6][3]) {
	// Based on Lagrange multiplier method
	// G(x) = x & F(x,y,z) = 0;
	// grad(F)=lambda*grad(G)

	// Notice
	// - May not be proper for cones
	// - If a corner point is not defined yet a curve line, an exact point would be
	//  defined later by the curve between this surface and {X = a|Y = a|Z = a}.

	// return three-digit value (XYZ)
	// 0 : none for along the axis
	// 1 : left only, 2 : right only, 3 : both sides

	int code = 0;

	if (IsCurve()) {
		double gradG[3][4] = {
			{2.0*c_xs, c_xy, c_xz, c_x},
			{c_xy, 2.0*c_ys, c_yz, c_y},
			{c_xz, c_yz, 2.0*c_zs, c_z} 
		};
		auto coeff1 = eq;
		for (int i = 0; i < 3; i++) {
			// (i = 0) z = c22[0][0]*x+c22[0][1] & y = c22[1][0]*x+c22[1][1]
			// (i = 1) x = c22[0][0]*y+c22[0][1] & z = c22[1][0]*y+c22[1][1]
			// (i = 2) y = c22[0][0]*z+c22[0][1] & x = c22[1][0]*z+c22[1][1]
			// (i = 0) (j,i1,i2) = (0, y, z), (1, z, y)
			// (i = 1)           = (0, z, x), (1, x, z)
			// (i = 2)           = (0, x, y), (1, y, x)
			double c22[2][2];
			int map[2] = { (i + 1) % 3,(i + 2) % 3 };

			bool iscontd = false;
			for (int j = 0; j < 2; j++) {
				// remove i1's term
				int i1 = map[j], i2 = map[(j+1)%2];
				double c0, c1, c2;
				c0 = gradG[i1][i] * gradG[i2][i1] - gradG[i2][i] * gradG[i1][i1]; // i's term
				c1 = gradG[i1][i2] * gradG[i2][i1] - gradG[i2][i2] * gradG[i1][i1]; //i2's term
				c2 = gradG[i1][3] * gradG[i2][i1] - gradG[i2][3] * gradG[i1][i1]; // constant term
				if (abs(c1) < eps_geo) iscontd = true; // *** The min or max points cannot be determined.
				c22[j][0] = c0 / c1; c22[j][1] = c2 / c1;
			}
			if (iscontd) continue;

			double xyzc[4][4] = { {0} };
			xyzc[i][i] = 1.;
			xyzc[map[0]][i] = c22[1][0]; xyzc[map[0]][3] = c22[1][1];
			xyzc[map[1]][i] = c22[0][0]; xyzc[map[1]][3] = c22[0][1];
			xyzc[3][3] = 1.0;
			coeff1 = CalCoeffs(eq, xyzc);

			double bc[2] = { coeff1[i + 6], coeff1[9] };
			if (abs(coeff1[i]) < eps_geo) {
				double soli = -bc[1] / bc[0];
				double soli1 = c22[1][0] * soli + c22[1][1], soli2 = c22[0][0] * soli + c22[0][1];
				double invlam = gradG[i][i] * soli + gradG[i][map[0]] * soli1 + gradG[i][map[1]] * soli2;
				if (abs(invlam) < eps_geo) {
					code += 3 * pow10[2 - i];
					sol[i * 2][i] = sol[i * 2 + 1][i] = soli;
					sol[i * 2][map[0]] = sol[i * 2 + 1][map[0]] = soli1;
					sol[i * 2][map[1]] = sol[i * 2 + 1][map[1]] = soli2;
				}
				else if (invlam > 0) {
					code += 2 * pow10[2 - i];
					sol[i * 2 + 1][i] = soli;
					sol[i * 2 + 1][map[0]] = soli1;
					sol[i * 2 + 1][map[1]] = soli2;
				}
				else {
					code += 1 * pow10[2 - i];
					sol[i * 2][i] = soli;
					sol[i * 2][map[0]] = soli1;
					sol[i * 2][map[1]] = soli2;
				}
				continue;
			}
			bc[0] /= coeff1[i]; bc[1] /= coeff1[i];
			double det = bc[0] * bc[0] - 4.0*bc[1];
			if (det > eps_geo) {
				code += 3 * pow10[2 - i];
				double soli = 0.5*(-bc[0] - sqrt(det));
				double soli1 = c22[1][0] * soli + c22[1][1], soli2 = c22[0][0] * soli + c22[0][1];
				sol[i * 2][i] = soli; sol[i * 2][map[0]] = soli1; sol[i * 2][map[1]] = soli2;
				soli = 0.5*(-bc[0] + sqrt(det));
				soli1 = c22[1][0] * soli + c22[1][1], soli2 = c22[0][0] * soli + c22[0][1];
				sol[i * 2 + 1][i] = soli; sol[i * 2 + 1][map[0]] = soli1; sol[i * 2 + 1][map[1]] = soli2;
			}
		}
	}
	return code;
}

int UnitSurf::GetLocalExPln(double CoeffPln[4], double sol[6][3]) {
	int code = 0, code2 = 0;
	if (IsCurve()) {
		for (int i = 0; i < 3; i++) {
			double xyzc[4][4] = { { 0.0, }, };
			for (int j = 0; j < 4; j++) xyzc[j][j] = 1.0;
			if (abs(CoeffPln[i]) > eps_geo) {
				double invc = 1. / CoeffPln[i];
				// Cancel i's term, x = ay+bz+c
				for (int j = 0; j < 4; j++) xyzc[i][j] = -CoeffPln[j] * invc;
				xyzc[i][i] = 0.;

				auto coeff1 = CalCoeffs(eq, xyzc);
				double coeff2[6] = { coeff1[(i + 1) % 3],coeff1[(i + 2) % 3],coeff1[3 + (i + 1) % 3],coeff1[6 + (i + 1) % 3],coeff1[6 + (i + 2) % 3],coeff1[9] };

				// Get local extremum for an axis
				double sol1[2][2];
				int code1 = GetLocalExCurve((i + 1) % 3, coeff2, sol1);
				if (code1 % 2 == 1) {
					sol[(i + 1) % 3 * 2][(i + 1) % 3] = sol1[0][0];
					sol[(i + 1) % 3 * 2][(i + 2) % 3] = sol1[0][1];
					sol[(i + 1) % 3 * 2][i] = -(CoeffPln[(i + 1) % 3] * sol1[0][0] + CoeffPln[(i + 2) % 3] * sol1[0][1] + CoeffPln[3]) / CoeffPln[i];
					//this->GetPointVar2(i, sol1[0], sol[(i + 1) % 3 * 2][i]);
				}
				if (code1 > 2) {
					sol[(i + 1) % 3 * 2 + 1][(i + 1) % 3] = sol1[1][0];
					sol[(i + 1) % 3 * 2 + 1][(i + 2) % 3] = sol1[1][1];
					sol[(i + 1) % 3 * 2 + 1][i] = -(CoeffPln[(i + 1) % 3] * sol1[1][0] + CoeffPln[(i + 2) % 3] * sol1[1][1] + CoeffPln[3]) / CoeffPln[i];
					//this->GetPointVar2(i, sol1[1], sol[(i + 1) % 3 * 2 + 1][i]);
				}
				code += pow10[2 - (i + 1) % 3] * code1;
				code2 += code1;
				// If (i+2)%3's term was zero, get another pair
				if (abs(CoeffPln[(i + 2) % 3]) > eps_geo) continue;
				code1 = GetLocalExCurve((i + 2) % 3, coeff2, sol1);
				if (code1 % 2 == 1) {
					sol[(i + 2) % 3 * 2][(i + 1) % 3] = sol1[0][0];
					sol[(i + 2) % 3 * 2][(i + 2) % 3] = sol1[0][1];
					sol[(i + 2) % 3 * 2][i] = -(CoeffPln[(i + 1) % 3] * sol1[0][0] + CoeffPln[(i + 2) % 3] * sol1[0][1] + CoeffPln[3]) / CoeffPln[i];
					//this->GetPointVar2(i, sol1[0], sol[(i + 2) % 3 * 2][i]);
				}
				if (code1 > 2) {
					sol[(i + 2) % 3 * 2 + 1][(i + 1) % 3] = sol1[1][0];
					sol[(i + 2) % 3 * 2 + 1][(i + 2) % 3] = sol1[1][1];
					sol[(i + 2) % 3 * 2 + 1][i] = -(CoeffPln[(i + 1) % 3] * sol1[1][0] + CoeffPln[(i + 2) % 3] * sol1[1][1] + CoeffPln[3]) / CoeffPln[i];
					//this->GetPointVar2(i, sol1[1], sol[(i + 2) % 3 * 2 + 1][i]);
				}
				code += pow10[2 - (i + 2) % 3] * code1;
				code2 += code1;
				if (code2 == 6) {
					double anyx = sol[(i + 1) % 3 * 2 + 1][0];
					double anyy = sol[(i + 1) % 3 * 2 + 1][1];
					double anyz = sol[(i + 1) % 3 * 2 + 1][2];
					sol[i * 2][0] = anyx; sol[i * 2][1] = anyy; sol[i * 2][2] = anyz;
					sol[i * 2 + 1][0] = anyx; sol[i * 2 + 1][1] = anyy; sol[i * 2 + 1][2] = anyz;
					code += pow10[2 - i] * 3;
				}
			}
		}
	}
	return code;
}

int UnitVol::OneIntersection(int idx, CartPlane CartPlane, 
	const array<double, 2>& val, array<double, 2>& sol) {
	double x, y, z;
	switch (CartPlane) {
	case CartPlane::XY:
		x = val[0]; y = val[1];
		break;
	case CartPlane::YZ:
		y = val[0]; z = val[1];
		break;
	case CartPlane::XZ:
		x = val[0]; z = val[1];
		break;
	}
	int ninter = Surfaces[idx].GetIntersection(CartPlane, val, sol);
	int minter = ninter;
	for (int i = 0; i < minter; i++) {
		double bufsol;
		if (i == 0) bufsol = sol[0];
		if (i == 1) bufsol = sol[1];
		switch (CartPlane) {
		case CartPlane::XY:
			z = bufsol;
			break;
		case CartPlane::YZ:
			x = bufsol;
			break;
		case CartPlane::XZ:
			y = bufsol;
			break;
		}
		for (int j = 0; j < Surfaces.size(); j++) {
			if (j == idx) continue;
			if (!Surfaces[j].IsInside(x, y, z, true)) {
				ninter--;
				if (ninter == 0) break;
				sol[0] = sol[1];
			}
		}
	}

	return ninter;
}

void UnitVol::Create(int _nsurf, const UnitSurf *Surfaces) {
	this->Surfaces.insert(this->Surfaces.begin(), Surfaces, Surfaces + _nsurf);
}

void UnitVol::Create(const vector<UnitSurf>& Surfaces)
{
	this->Surfaces = Surfaces;
}

void UnitVol::Create(const UnitSurf &_Surfaces) {
	this->Surfaces.emplace_back(_Surfaces);
}

void UnitVol::Create(const UnitVol& rhs)
{
	Surfaces = rhs.Surfaces;
	isbounded = rhs.isbounded;
	xlr = rhs.xlr; ylr = rhs.ylr; zlr = rhs.zlr;
	vol = rhs.vol;
}

void UnitVol::Relocate(double dx, double dy, double dz) {
	for (int i = 0; i < Surfaces.size(); i++) Surfaces[i].Relocate(dx, dy, dz);
	if (isbounded) {
		xlr.x += dx; xlr.y += dx;
		ylr.x += dy; ylr.y += dy;
		zlr.x += dz; zlr.y += dz;
	}
}

void UnitVol::Rotate(double cos, double sin, CartAxis Ax) {
	for (int i = 0; i < Surfaces.size(); i++) Surfaces[i].Rotate(cos, sin, Ax);
	if (isbounded) 
		CalBoundBox();
	
}

bool UnitVol::IsInside(double x, double y, double z, bool includeOn) {
	bool inside = true;
	for (int i = 0; i < Surfaces.size(); i++) 
		inside = inside && Surfaces[i].IsInside(x, y, z, includeOn);
	return inside;
}

int UnitVol::GetIntersection(CartPlane CartPlane, const array<double, 2>& val, vector<double>& sol) {
	vector<array<double, 2>> bufs(Surfaces.size());
	vector<int> ninters(Surfaces.size());

	int ntotint = 0;
	for (int i = 0; i < Surfaces.size(); i++) {
		ninters[i] = OneIntersection(i, CartPlane, val, bufs[i]);
		ntotint += ninters[i];
	}

	sol.resize(ntotint);

	ntotint = 0;
	for (int i = 0; i < Surfaces.size(); i++) {
		for (int j = 0; j < ninters[i]; j++) {
			sol[ntotint] = bufs[i][j];
			ntotint++;
		}
	}
	
	// Sorting the solutions
	// * not developed yet

	return ntotint;
}

void UnitVol::ResidentFilter(int &codes, int acode, double localEx[6][3], double corners[6][3]) {
	for (int axi = 0; axi < 3; axi++) {
		int globcode = codes % pow10[3 - axi] / pow10[2 - axi];
		int code0 = acode % pow10[3 - axi] / pow10[2 - axi], code1 = 0;
		if (code0 % 2 == 1) {
			if (IsInside(localEx[axi * 2][0], localEx[axi * 2][1], localEx[axi * 2][2], true)) {
				code1 += 1;
				if (globcode % 2 != 1 || localEx[axi * 2][axi] < corners[axi * 2][axi]) {
					std::copy(localEx[axi * 2], localEx[axi * 2]+3, corners[axi * 2]);
				}
			}
		}
		if (code0 > 1) {
			if (IsInside(localEx[axi * 2 + 1][0], localEx[axi * 2 + 1][1], localEx[axi * 2 + 1][2], true)) {
				code1 += 2;
				if (globcode < 2 || localEx[axi * 2 + 1][axi] > corners[axi * 2 + 1][axi]) {
					std::copy(localEx[axi * 2 + 1], localEx[axi * 2 + 1]+3, corners[axi * 2 + 1]);
				}
			}
		}
		codes += pow10[2 - axi] * (code1 - globcode);
	}
}

bool UnitVol::CalBoundBox() {
	double corners[6][3];
	int codes = 0;

	// Symbol for bound-by, (,)
	// o C : Curved surface
	// o P : Plane surface

	// [1] Curved surface local extremum : (C)
	//     1. Get local extremum
	//     2. Extract what are inside volume
	//     3. Revise codes and corner points
	for (int i = 0; i < Surfaces.size(); i++) {
		double localEx[6][3] = { { 0.0, }, };
		int acode = Surfaces[i].GetLocalExSelf(localEx);
		if (acode == 0) continue;
		ResidentFilter(codes, acode, localEx, corners);
	}

	// [2] Curve local extremum between curved surfaces and planes : (C,P)
	//     o Same procedures as above
	for (int i = 0; i < Surfaces.size(); i++) {
		if (!Surfaces[i].IsCurve()) continue;
		double localEx[6][3] = { { 0.0 }, };
		for (int j = 0; j < Surfaces.size(); j++) {
			if (i == j) continue;
			const auto& coeffj = Surfaces[j].GetEquation();
			if (Surfaces[j].IsCurve()) continue;
			double coeffPln[4] = { coeffj[6],coeffj[7],coeffj[8],coeffj[9] };
			int acode = Surfaces[i].GetLocalExPln(coeffPln, localEx);
			if (acode == 0) continue;
			ResidentFilter(codes, acode, localEx, corners);
		}
	}

	// [3] Curve local extremum between two curved surfaces : (C,C)
	//     o Same procedures as above
	//     ** Not developed yet. This will be done after C5G7 hard-code progresses.

	// From [1] to [3] steps, the local extremum points are obtained.
	// However, from [4] to [5], only the intersection points can be obtained.
	// It is impossible to tell whether each is left or right boundaries of any axes.
	// Explicit comparisons which is largest or smallest along all directions should be carried out.

	queue<double3> Inters;

	// [4] Intersection of each curved surface and two planes : (C,P,P)
	//     1. Get a line from a pair of planes attached to each curved surface
	//     2. Find two intersections to the curved surface
	//     3. Follow [1]-2 and [1]-3
	//     ** Not developed either.

	// [5] Intersection of three planes : (P,P,P)
	//     1. Get the intersection points from triple planes
	//     2. Follow [1]-2 and [1]-3
	for (int i = 0; i < Surfaces.size(); i++) {
		const auto& coeff0 = Surfaces[i].GetEquation();
		if (Surfaces[i].IsCurve()) continue;
		double coeff0_4[4] = { coeff0[6], coeff0[7], coeff0[8], coeff0[9] };
		for (int j = i + 1; j < Surfaces.size(); j++) {
			const auto& coeff1 = Surfaces[j].GetEquation();
			if (Surfaces[j].IsCurve()) continue;
			double coeff1_4[4] = { coeff1[6], coeff1[7], coeff1[8], coeff1[9] };
			for(int k = j + 1; k < Surfaces.size(); k++) {
				const auto& coeff2 = Surfaces[k].GetEquation();
				if (Surfaces[k].IsCurve()) continue;
				double coeff2_4[4] = { coeff2[6], coeff2[7], coeff2[8], coeff2[9] };
				double soltriple[3];
				if (UnitSurf::GetTriplePoint(coeff0_4, coeff1_4, coeff2_4, soltriple)) {
					double3 sol3;
					sol3.x = soltriple[0]; sol3.y = soltriple[1]; sol3.z = soltriple[2];
					Inters.push(sol3);
				}
			}
		}
	}

	// Filter residents and determine bounding values
	// The corner points from local extremum are collapsed to the queue, Inters.
	// If a plane composed of three resident points holds all the others, return false
	queue<double3> Residents;
	for (int i = 0; i < 3; i++) {
		int acode = codes % pow10[3 - i] / pow10[2 - i];
		double3 tmp; 
		tmp.x = corners[2 * i][0]; tmp.y = corners[2 * i][1]; tmp.z = corners[2 * i][2];
		if (acode % 2 == 1) Residents.push(tmp);
		tmp.x = corners[2 * i + 1][0]; tmp.y = corners[2 * i + 1][1]; tmp.z = corners[2 * i + 1][2];
		if (acode > 1) Residents.push(tmp);
	}
	while (!Inters.empty()) {
		//double aninter[3];
		//std::copy(Inters.front(), Inters.front() + 3, aninter);
		double3 aninter;
		aninter = Inters.front();
		if (IsInside(aninter.x, aninter.y, aninter.z, true)) {
			Residents.push(aninter);
		}
		Inters.pop();
	}
	if (Residents.empty()) return false;
	//double aninter0[3];
	//std::copy(Residents.front(), Residents.front() + 3, aninter0);
	double3 aninter0;
	aninter0 = Residents.front();

	double boundval[6] = { aninter0.x, aninter0.x, aninter0.y, aninter0.y, aninter0.z, aninter0.z };
	int seq = 0;
	double DirVec[3];
	while (!Residents.empty()) {
		//double aninter[3];
		//std::copy(Residents.front(), Residents.front() + 3, aninter);
		double3 aninter;
		aninter = Residents.front();
		// Boundary values update
		double vn[3] = { aninter.x, aninter.y, aninter.z };
		for (int i = 0; i < 3; i++) {
			double vl = boundval[i * 2], vr = boundval[i * 2 + 1];
			boundval[i * 2] = (vl > vn[i]) ? vn[i] : vl;
			boundval[i * 2 + 1] = (vr < vn[i]) ? vn[i] : vr;
		}
		// Boundedness check
		double DirVec1[3] = { aninter.x - aninter0.x, aninter.y - aninter0.y, aninter.z - aninter0.z };
		double det = DirVec1[0] * DirVec1[0] + DirVec1[1] * DirVec1[1] + DirVec1[2] * DirVec1[2];
		if (det > eps_geo) {
			switch (seq) {
			case 0:
				std::copy(DirVec1, DirVec1 + 3, DirVec);
				seq = 1;
				break;
			case 1:
				double3 bufdir; bufdir.x = DirVec1[0]; bufdir.y = DirVec1[1]; bufdir.z = DirVec1[2];
				DirVec1[0] = DirVec[1] * bufdir.z - DirVec[2] * bufdir.y;
				DirVec1[1] = DirVec[2] * bufdir.x - DirVec[0] * bufdir.z;
				DirVec1[2] = DirVec[0] * bufdir.y - DirVec[1] * bufdir.x;
				det = DirVec1[0] * DirVec1[0] + DirVec1[1] * DirVec1[1] + DirVec1[2] * DirVec1[2];
				if (det > eps_geo) {
					std::copy(DirVec1, DirVec1 + 3, DirVec);
					seq = 2;
				}
				break;
			case 2:
				det = DirVec[0] * DirVec1[0] + DirVec[1] * DirVec1[1] + DirVec[2] * DirVec1[2];
				if (abs(det) < eps_geo) seq = 3;
				break;
			case 3:
				break;
			}
		}
		Residents.pop();
	}
	if (seq < 3) return false;
	xlr.x = boundval[0]; xlr.y = boundval[1];
	ylr.x = boundval[2]; ylr.y = boundval[3];
	zlr.x = boundval[4]; zlr.y = boundval[5];
	return true;
}

double UnitVol::CalVolume(int nx, int ny, int nz) {
	if (!isbounded) return 0.0;
	double lx = xlr.y - xlr.x, ly = ylr.y - ylr.x, lz = zlr.y - zlr.x;
	double dx = lx / (double) nx, dy = ly / (double) ny, dz = lz / (double) nz;

	vector<double> xpts(nx), ypts(ny), zpts(nz);

	for (int ix = 0; ix < nx; ix++) xpts[ix] = xlr.x + dx * 0.5 * (double)(ix * 2 + 1);
	for (int iy = 0; iy < ny; iy++) ypts[iy] = ylr.x + dy * 0.5 * (double)(iy * 2 + 1);
	for (int iz = 0; iz < nz; iz++) zpts[iz] = zlr.x + dz * 0.5 * (double)(iz * 2 + 1);

	double volsumYZ = 0.;// , volsum0 = 0.;
	#pragma omp parallel for reduction(+:volsumYZ)
	for (int iy = 0; iy < ny; iy++) {
		//volsum0 = 0.;
		for (int iz = 0; iz < nz; iz++) {
			array<double, 2> vals = { ypts[iy], zpts[iz] };
			vector<double> sol;
			int ninter = GetIntersection(CartPlane::YZ, vals, sol);
			double volsum1 = 0.;
			for (int i = 1; i < ninter; i+=2) {
				volsum1 += (sol[i] - sol[i-1]);
			}
			//volsumYZ += 0.5*(volsum1 + volsum0);
			//volsum0 = volsum1;
			volsumYZ += volsum1;
		}
	}
	volsumYZ *= dy * dz;

	double volsumXZ = 0.;
	#pragma omp parallel for reduction(+:volsumXZ)
	for (int ix = 0; ix < nx; ix++) {
		//volsum0 = 0.;
		for (int iz = 0; iz < nz; iz++) {
			array<double, 2> vals = { xpts[ix], zpts[iz] };
			vector<double> sol;
			int ninter = GetIntersection(CartPlane::XZ, vals, sol);
			double volsum1 = 0.;
			for (int i = 1; i < ninter; i += 2) {
				volsum1 += (sol[i] - sol[i - 1]);
			}
			//volsumXZ += 0.5*(volsum1 + volsum0);
			//volsum0 = volsum1;
		  volsumXZ += volsum1;
		}
	}
	volsumXZ *= dx * dz;

	double volsumXY = 0.;
	#pragma omp parallel for reduction(+:volsumXY)
	for (int ix = 0; ix < nx; ix++) {
		//volsum0 = 0.;
		for (int iy = 0; iy < ny; iy++) {
			array<double, 2> vals = { xpts[ix],ypts[iy] };
			vector<double> sol;
			int ninter = GetIntersection(CartPlane::XY, vals, sol);
			double volsum1 = 0.;
			for (int i = 1; i < ninter; i += 2) {
				volsum1 += (sol[i] - sol[i - 1]);
			}
			//volsumXY += 0.5*(volsum1 + volsum0);
			//volsum0 = volsum1;
			volsumXY += volsum1;
		}
	}
	volsumXY *= dx * dy;

	return (volsumXY+volsumXZ+volsumYZ)/3.0;
}

void UnitVol::Finalize() {
	isbounded = CalBoundBox();
	if (!isbounded) {
		// Exception
	}
	else {
		vol = CalVolume(50,50,50);
	}
}

bool UnitVol::GetBoundBox(double2& xlr, double2& ylr, double2& zlr) const {
	xlr = this->xlr; ylr = this->ylr; zlr = this->zlr;
	return isbounded;
}

void UnitComp::Create(int imat, const UnitVol &Volume) {
	this->UnitVol::operator=(Volume);
	this->imat = imat;
	// * Calculate Bound Box of UnitComp
}

void UnitComp::Create(const UnitComp& rhs)
{
	imat = rhs.imat;
	this->UnitVol::operator=(rhs);
}

}