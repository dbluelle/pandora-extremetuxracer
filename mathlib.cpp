/* --------------------------------------------------------------------
EXTREME TUXRACER

Copyright (C) 1999-2001 Jasmin F. Patry (Tuxracer)
Copyright (C) 2004-2005 Volker Stroebel (Planetpenguin Racer)
Copyright (C) 2010 Extreme Tuxracer Team

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
---------------------------------------------------------------------*/

#ifdef HAVE_CONFIG_H
#include <etr_config.h>
#endif

#include "mathlib.h"
#include <algorithm>

ETR_DOUBLE VectorLength (const TVector3 &v) {
	return sqrt (MAG_SQD(v));
}

ETR_DOUBLE DotProduct (const TVector3& v1, const TVector3& v2) {
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

TVector3 ScaleVector (ETR_DOUBLE s, const TVector3& v) {
	TVector3 rval(
		s * v.x,
		s * v.y,
		s * v.z);
	return rval;
}

TVector3 AddVectors (const TVector3& v1, const TVector3& v2) {
	TVector3 result(
		v1.x + v2.x,
		v1.y + v2.y,
		v1.z + v2.z);
	return result;
}

TVector3 SubtractVectors (const TVector3& v1, const TVector3& v2) {
	TVector3 result(
		v1.x - v2.x,
		v1.y - v2.y,
		v1.z - v2.z);
	return result;
}

TVector3 ProjectToPlane (const TVector3& nml, const TVector3& v) {
	ETR_DOUBLE dotProd = DotProduct (nml, v);
	TVector3 nmlComp = ScaleVector (dotProd, nml);

	return SubtractVectors (v, nmlComp);
}

ETR_DOUBLE NormVector (TVector3 &v) {
	if (v.x == 0 && v.y == 0 && v.z == 0) return 0.0;
	ETR_DOUBLE denom = sqrt (v.x * v.x + v.y * v.y + v.z * v.z);
	v = ScaleVector (1.0 / denom, v);
	return denom;
}


TVector3 CrossProduct(const TVector3& u, const TVector3& v) {
	TVector3 ret(
		u.y * v.z - u.z * v.y,
		u.z * v.x - u.x * v.z,
		u.x * v.y - u.y * v.x);
	return ret;
}

TVector3 TransformVector(const TMatrix mat, const TVector3& v) {
	TVector3 r;
	r.x = v.x * mat[0][0] + v.y * mat[1][0] + v.z * mat[2][0];
	r.y = v.x * mat[0][1] + v.y * mat[1][1] + v.z * mat[2][1];
	r.z = v.x * mat[0][2] + v.y * mat[1][2] + v.z * mat[2][2];
	return r;
}

TVector3 TransformNormal(const TVector3& n, const TMatrix mat) {
	TVector3 r;
	r.x = n.x * mat[0][0] + n.y * mat[0][1] + n.z * mat[0][2];
	r.y = n.x * mat[1][0] + n.y * mat[1][1] + n.z * mat[1][2];
	r.z = n.x * mat[2][0] + n.y * mat[2][1] + n.z * mat[2][2];
	return r;
}

TVector3 TransformPoint(const TMatrix mat, const TVector3& p) {
	TVector3 r;
	r.x = p.x * mat[0][0] + p.y * mat[1][0] + p.z * mat[2][0];
	r.y = p.x * mat[0][1] + p.y * mat[1][1] + p.z * mat[2][1];
	r.z = p.x * mat[0][2] + p.y * mat[1][2] + p.z * mat[2][2];
	r.x += mat[3][0];
	r.y += mat[3][1];
	r.z += mat[3][2];
	return r;
}

TPlane MakePlane (ETR_DOUBLE nx, ETR_DOUBLE ny, ETR_DOUBLE nz, ETR_DOUBLE d) {
	TPlane tmp;
	tmp.nml.x = nx;
	tmp.nml.y = ny;
	tmp.nml.z = nz;
	tmp.d = d;

	return tmp;
}

bool IntersectPlanes (const TPlane& s1, const TPlane& s2, const TPlane& s3, TVector3 *p) {
	ETR_DOUBLE A[3][4];
	ETR_DOUBLE x[3];
	ETR_DOUBLE retval;

	A[0][0] =  s1.nml.x;
	A[0][1] =  s1.nml.y;
	A[0][2] =  s1.nml.z;
	A[0][3] = -s1.d;

	A[1][0] =  s2.nml.x;
	A[1][1] =  s2.nml.y;
	A[1][2] =  s2.nml.z;
	A[1][3] = -s2.d;

	A[2][0] =  s3.nml.x;
	A[2][1] =  s3.nml.y;
	A[2][2] =  s3.nml.z;
	A[2][3] = -s3.d;

	retval = Gauss ((ETR_DOUBLE*) A, 3, x);

	if (retval != 0) {
		return false;
	} else {
		p->x = x[0];
		p->y = x[1];
		p->z = x[2];
		return true;
	}
}

ETR_DOUBLE DistanceToPlane (const TPlane& plane, const TVector3& pt) {
	return
		plane.nml.x * pt.x +
		plane.nml.y * pt.y +
		plane.nml.z * pt.z +
		plane.d;
}

void MakeIdentityMatrix(TMatrix h) {
	for (int i= 0; i< 4; i++)
		for (int j= 0; j< 4; j++)
			h[i][j]= (i==j);
}


void MultiplyMatrices (TMatrix ret, const TMatrix mat1, const TMatrix mat2) {
	TMatrix r;

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			r[j][i] = mat1[0][i] * mat2[j][0] +
				mat1[1][i] * mat2[j][1] +
				mat1[2][i] * mat2[j][2] +
				mat1[3][i] * mat2[j][3];

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			ret[i][j] = r[i][j];
}


void TransposeMatrix (const TMatrix mat, TMatrix trans) {
	TMatrix r;

	for (int i = 0; i < 4; i++)
		for (int j = 0; j< 4; j++)
			r[j][i] = mat[i][j];

	for (int i = 0; i< 4; i++)
		for (int j = 0; j< 4; j++)
			trans[i][j] = r[i][j];
}

void MakeRotationMatrix (TMatrix mat, ETR_DOUBLE angle, char axis) {
	ETR_DOUBLE sinv, cosv;
	sinv = sin (ANGLES_TO_RADIANS (angle));
	cosv = cos (ANGLES_TO_RADIANS (angle));

	MakeIdentityMatrix (mat);

	switch (axis) {
	case 'x':
		mat[1][1] = cosv;
		mat[2][1] = -sinv;
		mat[1][2] = sinv;
		mat[2][2] = cosv;
		break;

	case 'y':
		mat[0][0] = cosv;
		mat[2][0] = sinv;
		mat[0][2] = -sinv;
		mat[2][2] = cosv;
		break;

	case 'z':
		mat[0][0] = cosv;
		mat[1][0] = -sinv;
		mat[0][1] = sinv;
		mat[1][1] = cosv;
		break;
	}
}

void MakeTranslationMatrix (TMatrix mat, ETR_DOUBLE x, ETR_DOUBLE y, ETR_DOUBLE z) {
	MakeIdentityMatrix (mat);
	mat[3][0] = x;
	mat[3][1] = y;
	mat[3][2] = z;
}

void MakeScalingMatrix (TMatrix mat, ETR_DOUBLE x, ETR_DOUBLE y, ETR_DOUBLE z) {
	MakeIdentityMatrix (mat);
	mat[0][0] = x;
	mat[1][1] = y;
	mat[2][2] = z;
}

void MakeBasisMat (TMatrix mat, const TVector3& w1, const TVector3& w2, const TVector3& w3) {
	MakeIdentityMatrix (mat);
	mat[0][0] = w1.x;
	mat[0][1] = w1.y;
	mat[0][2] = w1.z;
	mat[1][0] = w2.x;
	mat[1][1] = w2.y;
	mat[1][2] = w2.z;
	mat[2][0] = w3.x;
	mat[2][1] = w3.y;
	mat[2][2] = w3.z;
}

void MakeBasismatrix_Inv (TMatrix mat, TMatrix invMat,
		const TVector3& w1, const TVector3& w2, const TVector3& w3) {
	MakeIdentityMatrix (mat);
	mat[0][0] = w1.x;
	mat[0][1] = w1.y;
	mat[0][2] = w1.z;
	mat[1][0] = w2.x;
	mat[1][1] = w2.y;
	mat[1][2] = w2.z;
	mat[2][0] = w3.x;
	mat[2][1] = w3.y;
	mat[2][2] = w3.z;

	MakeIdentityMatrix (invMat);
	invMat[0][0] = w1.x;
	invMat[1][0] = w1.y;
	invMat[2][0] = w1.z;
	invMat[0][1] = w2.x;
	invMat[1][1] = w2.y;
	invMat[2][1] = w2.z;
	invMat[0][2] = w3.x;
	invMat[1][2] = w3.y;
	invMat[2][2] = w3.z;
}

void RotateAboutVectorMatrix (TMatrix mat, const TVector3& u, ETR_DOUBLE angle) {
	TMatrix rx, irx, ry, iry;

	ETR_DOUBLE a = u.x;
	ETR_DOUBLE b = u.y;
	ETR_DOUBLE c = u.z;

	ETR_DOUBLE d = sqrt (b*b + c*c);

	if (d < EPS) {
		if (a < 0)
			MakeRotationMatrix (mat, -angle, 'x');
		else
			MakeRotationMatrix (mat, angle, 'x');
		return;
	}

	MakeIdentityMatrix (rx);
	MakeIdentityMatrix (irx);
	MakeIdentityMatrix (ry);
	MakeIdentityMatrix (iry);

	rx[1][1] = c/d;
	rx[2][1] = -b/d;
	rx[1][2] = b/d;
	rx[2][2] = c/d;

	irx[1][1] = c/d;
	irx[2][1] = b/d;
	irx[1][2] = -b/d;
	irx[2][2] = c/d;

	ry[0][0] = d;
	ry[2][0] = -a;
	ry[0][2] = a;
	ry[2][2] = d;

	iry[0][0] = d;
	iry[2][0] = a;
	iry[0][2] = -a;
	iry[2][2] = d;

	MakeRotationMatrix (mat, angle, 'z');

	MultiplyMatrices (mat, mat, ry);
	MultiplyMatrices (mat, mat, rx);
	MultiplyMatrices (mat, iry, mat);
	MultiplyMatrices (mat, irx, mat);
}

TQuaternion MultiplyQuaternions (const TQuaternion& q, const TQuaternion& r) {
	TQuaternion res;
	res.x = q.y * r.z - q.z * r.y + r.w * q.x + q.w * r.x;
	res.y = q.z * r.x - q.x * r.z + r.w * q.y + q.w * r.y;
	res.z = q.x * r.y - q.y * r.x + r.w * q.z + q.w * r.z;
	res.w = q.w * r.w - q.x * r.x - q.y * r.y - q.z * r.z;
	return res;
}

TQuaternion AddQuaternions (const TQuaternion& q, const TQuaternion& r) {
	TQuaternion res(
		q.x + r.x,
		q.y + r.y,
		q.z + r.z,
		q.w + r.w);
	return res;
}

TQuaternion ConjugateQuaternion (const TQuaternion& q) {
	TQuaternion res(
		-1 * q.x,
		-1 * q.y,
		-1 * q.z,
		q.w);

	return res;
}

TQuaternion ScaleQuaternion (ETR_DOUBLE s, const TQuaternion& q) {
	TQuaternion res(
		s * q.x,
		s * q.y,
		s * q.z,
		s * q.w);

	return res;
}

void MakeMatrixFromQuaternion (TMatrix mat, const TQuaternion& q) {
	mat[0][0] = 1.0 - 2.0 *  (q.y * q.y + q.z * q.z);
	mat[1][0] =       2.0 *  (q.x * q.y - q.w * q.z);
	mat[2][0] =       2.0 *  (q.x * q.z + q.w * q.y);

	mat[0][1] =       2.0 *  (q.x * q.y + q.w * q.z);
	mat[1][1] = 1.0 - 2.0 *  (q.x * q.x + q.z * q.z);
	mat[2][1] =       2.0 *  (q.y * q.z - q.w * q.x);

	mat[0][2] =       2.0 *  (q.x * q.z - q.w * q.y);
	mat[1][2] =       2.0 *  (q.y * q.z + q.w * q.x);
	mat[2][2] = 1.0 - 2.0 *  (q.x * q.x + q.y * q.y);

	mat[3][0] = mat[3][1] = mat[3][2] = 0.0;
	mat[0][3] = mat[1][3] = mat[2][3] = 0.0;
	mat[3][3] = 1.0;
}

TQuaternion MakeQuaternionFromMatrix (const TMatrix m) {
	TQuaternion res;
	ETR_DOUBLE  tr, s, q[4];
	int    i, j, k;

	static int nxt[3] = {1, 2, 0};

	tr = m[0][0] + m[1][1] + m[2][2];

	if (tr > 0.0) {
		s = sqrt (tr + 1.0);
		res.w = 0.5 * s;
		s = 0.5 / s;
		res.x = (m[1][2] - m[2][1]) * s;
		res.y = (m[2][0] - m[0][2]) * s;
		res.z = (m[0][1] - m[1][0]) * s;
	} else {
		int i = 0;
		if (m[1][1] > m[0][0]) i = 1;
		if (m[2][2] > m[i][i]) i = 2;
		int j = nxt[i];
		int k = nxt[j];

		s = sqrt (m[i][i] - m[j][j] - m[k][k] + 1.0);

		q[i] = s * 0.5;

		if (s != 0.0) s = 0.5 / s;

		q[3] = (m[j][k] - m[k][j]) * s;
		q[j] = (m[i][j] + m[j][i]) * s;
		q[k] = (m[i][k] + m[k][i]) * s;

		res.x = q[0];
		res.y = q[1];
		res.z = q[2];
		res.w = q[3];
	}

	return res;
}

TQuaternion MakeRotationQuaternion (const TVector3& s, const TVector3& t) {
	TVector3 u = CrossProduct (s, t);
	ETR_DOUBLE sin2phi = NormVector (u);

	if (sin2phi < EPS) {
		return TQuaternion (0., 0., 0., 1.);
	} else {
		ETR_DOUBLE cos2phi = DotProduct (s, t);

		ETR_DOUBLE sinphi = sqrt ( (1 - cos2phi) / 2.0);
		ETR_DOUBLE cosphi = sqrt ( (1 + cos2phi) / 2.0);

		return TQuaternion(
			sinphi * u.x,
			sinphi * u.y,
			sinphi * u.z,
			cosphi);
	}
}

TQuaternion InterpolateQuaternions (const TQuaternion& q, TQuaternion r, ETR_DOUBLE t) {
	TQuaternion res;
	ETR_DOUBLE cosphi;
	ETR_DOUBLE sinphi;
	ETR_DOUBLE phi;
	ETR_DOUBLE scale0, scale1;

	cosphi = q.x * r.x + q.y * r.y + q.z * r.z + q.w * r.w;

	if (cosphi < 0.0) {
		cosphi = -cosphi;
		r.x = -r.x;
		r.y = -r.y;
		r.z = -r.z;
		r.w = -r.w;
	}

	if (1.0 - cosphi > EPS) {
		phi = acos (cosphi);
		sinphi = sin (phi);
		scale0 = sin (phi *  (1.0 - t)) / sinphi;
		scale1 = sin (phi * t) / sinphi;
	} else {
		scale0 = 1.0 - t;
		scale1 = t;
	}

	res.x = scale0 * q.x + scale1 * r.x;
	res.y = scale0 * q.y + scale1 * r.y;
	res.z = scale0 * q.z + scale1 * r.z;
	res.w = scale0 * q.w + scale1 * r.w;

	return res;
}

TVector3 RotateVector (const TQuaternion& q, const TVector3& v) {
	TQuaternion p(v.x, v.y, v.z, 1.0);

	TQuaternion qs(-q.x, -q.y, -q.z, q.w);

	TQuaternion res_q = MultiplyQuaternions (q, MultiplyQuaternions (p, qs));

	return res_q;
}

// --------------------------------------------------------------------
//				 Gauss
// --------------------------------------------------------------------

unsigned short order (ETR_DOUBLE *matrix, int n, int pivot);
void elim (ETR_DOUBLE *matrix, int n, int pivot);
void backsb (ETR_DOUBLE *matrix, int n, ETR_DOUBLE *soln);

int Gauss(ETR_DOUBLE *matrix, int n, ETR_DOUBLE *soln) {
	int pivot=0;
	unsigned short error=0;

	while ((pivot<(n-1)) && (!error)) {
		if(!(error = order(matrix,n,pivot))) {
			elim(matrix,n,pivot);
			pivot++;
		}
	}
	if (error) {
		return 1;
	} else {
		backsb(matrix, n, soln);
	}
	return 0;
}

unsigned short order (ETR_DOUBLE *matrix, int n, int pivot) {
	unsigned short error=0;

	int rmax = pivot;

	for (int row=pivot+1; row<n; row++) {
		if (fabs(*(matrix+row*(n+1)+pivot)) > fabs(*(matrix+rmax*(n+1)+pivot)))
			rmax = row;
	}

	if (fabs(*(matrix+rmax*(n+1)+pivot)) < EPS)
		error = 1;
	else if (rmax != pivot) {
		for (int k=0; k<(n+1); k++)
		{
			ETR_DOUBLE temp = *(matrix+rmax*(n+1)+k);
			*(matrix+rmax*(n+1)+k) = *(matrix+pivot*(n+1)+k);
			*(matrix+pivot*(n+1)+k) = temp;
		}
	}
	return error;
}

void elim (ETR_DOUBLE *matrix, int n, int pivot) {
	for (int row = pivot+1; row < n; row++) {
		ETR_DOUBLE factor = (*(matrix+row*(n+1)+pivot))/(*(matrix+pivot*(n+1)+pivot));
		*(matrix+row*(n+1)+pivot)=0.0;
		for (int col=pivot+1l; col<n+1; col++) {
			*(matrix+row*(n+1)+col) = *(matrix+row*(n+1)+col) -
				(*(matrix+pivot*(n+1)+col))*factor;
		}
	}
}


void backsb (ETR_DOUBLE *matrix, int n, ETR_DOUBLE *soln) {
	for (int row = n-1; row >=0; row--) {
		for (int col = n-1; col >= row+1; col--) {
			*(matrix+row*(n+1)+(n)) = *(matrix+row*(n+1)+n) -
				(*(soln+col))*(*(matrix+row*(n+1)+col));
		}
		*(soln+row) = (*(matrix+row*(n+1)+n))/(*(matrix+row*(n+1)+row));
	}
}

// ***************************************************************************
// ***************************************************************************

bool IntersectPolygon (const TPolygon& p, TVector3 *v) {
	TRay ray;
	ETR_DOUBLE d, s, nuDotProd, wec;
	ETR_DOUBLE distsq;

	TVector3 nml = MakeNormal (p, v);
	ray.pt = TVector3(0., 0., 0.);
	ray.vec = nml;

	nuDotProd = DotProduct (nml, ray.vec);
	if  (fabs(nuDotProd) < EPS)
		return false;

	d = - (nml.x * v[p.vertices[0]].x +
		nml.y * v[p.vertices[0]].y +
		nml.z * v[p.vertices[0]].z);

	if (fabs (d) > 1) return false;

	for (int i=0; i < p.num_vertices; i++) {
		TVector3 *v0, *v1;

		v0 = &v[p.vertices[i]];
		v1 = &v[p.vertices[ (i+1) % p.num_vertices ]];

		TVector3 edge_vec = SubtractVectors (*v1, *v0);
		ETR_DOUBLE edge_len = NormVector (edge_vec);

		ETR_DOUBLE t = - DotProduct (*v0, edge_vec);

		if (t < 0) {
			distsq = MAG_SQD (*v0);
		} else if (t > edge_len) {
			distsq = MAG_SQD (*v1);
		} else {
			*v0 = AddVectors (*v0, ScaleVector (t, edge_vec));
			distsq = MAG_SQD (*v0);
		}

		if (distsq <= 1) return true;
	}

	s = - (d + DotProduct (nml, TVector3(ray.pt.x, ray.pt.y, ray.pt.z))) / nuDotProd;
	TVector3 pt = AddVectors (ray.pt, ScaleVector (s, ray.vec));

	for (int i=0; i < p.num_vertices; i++) {
		TVector3 edge_nml = CrossProduct (nml,
			SubtractVectors (v[p.vertices[ (i+1) % p.num_vertices ]], v[p.vertices[i]]));

		wec = DotProduct (SubtractVectors (pt, v[p.vertices[i]]), edge_nml);
		if (wec < 0) return false;
	}
	return true;
}

bool IntersectPolyhedron (const TPolyhedron& p) {
	bool hit = false;
	for (size_t i=0; i<p.num_polygons; i++) {
		hit = IntersectPolygon (p.polygons[i], p.vertices);
		if (hit == true) break;
	}
	return hit;
}

TVector3 MakeNormal (const TPolygon& p, TVector3 *v) {
	TVector3 v1 = SubtractVectors (v[p.vertices[1]], v[p.vertices[0]]);
	TVector3 v2 = SubtractVectors (v[p.vertices[p.num_vertices-1]], v[p.vertices[0]]);
	TVector3 normal = CrossProduct (v1, v2);

	NormVector (normal);
	return normal;
}


TPolyhedron CopyPolyhedron (const TPolyhedron& ph) {
	TPolyhedron newph = ph;
	newph.vertices = new TVector3[ph.num_vertices];
#ifdef PANDORA
	TVector3* pnew =newph.vertices;
	TVector3* p =ph.vertices;
	if (ph.num_vertices > 0) {
		*pnew++ = *p;
		for (size_t i = 1; i < ph.num_vertices; ++i) {
			*pnew++ = *++p;
		}
	}
#else
	copy_n(ph.vertices, ph.num_vertices, newph.vertices);
#endif
	return newph;
}

void FreePolyhedron (const TPolyhedron& ph) {
	delete[] ph.vertices;
}

void TransPolyhedron (const TMatrix mat, const TPolyhedron& ph) {
	for (size_t i=0; i<ph.num_vertices; i++)
		ph.vertices[i] = TransformPoint (mat, ph.vertices[i]);
}

// --------------------------------------------------------------------
//					ode solver
// --------------------------------------------------------------------

const ETR_DOUBLE ode23_time_step_mat[] = { 0., 1./2., 3./4., 1. };
const ETR_DOUBLE ode23_coeff_mat[][4] = {
	{0.0, 1./2.,   0.0,  2./9.},
	{0.0,   0.0, 3./4.,  1./3.},
	{0.0,   0.0,   0.0,  4./9.},
	{0.0,   0.0,   0.0,    0.0}
};

const ETR_DOUBLE ode23_error_mat[] = {-5./72., 1./12., 1./9., -1./8. };
const ETR_DOUBLE ode23_time_step_exp = 1./3.;

int ode23_NumEstimates() {return 4; }

void ode23_InitOdeData (TOdeData *data, ETR_DOUBLE init_val, ETR_DOUBLE h) {
	data->init_val = init_val;
	data->h = h;
}

ETR_DOUBLE ode23_NextTime(TOdeData *data, int step) {
	return ode23_time_step_mat[step] * data->h;
}

ETR_DOUBLE ode23_NextValue (TOdeData *data, int step) {
	ETR_DOUBLE val = data->init_val;

	for (int i=0; i<step; i++)
		val += ode23_coeff_mat[i][step] * data->k[i];
	return val;
}

void ode23_UpdateEstimate(TOdeData *data, int step, ETR_DOUBLE val) {
	data->k[step] = data->h * val;
}

ETR_DOUBLE ode23_FinalEstimate(TOdeData *data) {
	ETR_DOUBLE val = data->init_val;

	for (int i=0; i<3; i++)
		val += ode23_coeff_mat[i][3] * data->k[i];
	return val;
}

ETR_DOUBLE ode23_EstimateError(TOdeData *data) {
	ETR_DOUBLE err=0.;

	for (int i=0; i<4; i++)
		err += ode23_error_mat[i] * data->k[i];
	return fabs(err);
}

ETR_DOUBLE ode23_TimestepExponent(){
	return ode23_time_step_exp;
}

TOdeSolver NewOdeSolver23() {
	TOdeSolver s;
	s.NumEstimates = ode23_NumEstimates;
	s.InitOdeData = ode23_InitOdeData;
	s.NextTime = ode23_NextTime;
	s.NextValue = ode23_NextValue;
	s.UpdateEstimate = ode23_UpdateEstimate;
	s.FinalEstimate = ode23_FinalEstimate;
	s.EstimateError = ode23_EstimateError;
	s.TimestepExponent = ode23_TimestepExponent;
	return s;
}

ETR_DOUBLE LinearInterp (const ETR_DOUBLE x[], const ETR_DOUBLE y[], ETR_DOUBLE val, int n) {
	int i;
	ETR_DOUBLE m, b;

	if (val < x[0]) i = 0;
	else if (val >= x[n-1]) i = n-2;
	else for (i=0; i<n-1; i++) if (val < x[i+1]) break;

	m = (y[i+1] - y[i]) / (x[i+1] - x[i]);
	b = y[i] - m * x[i];

	return m * val + b;
}

ETR_DOUBLE XRandom (float min, float max) {
	return (ETR_DOUBLE)rand () / RAND_MAX * (max - min) + min;
}

ETR_DOUBLE FRandom () {
	return (ETR_DOUBLE)rand () / RAND_MAX;
}

int IRandom (int min, int max) {
	return min + rand()%(max-min+1);
}

int ITrunc (int val, int base) {
	return (int)(val / base);
}

int IFrac (int val, int base) {
	return val - ITrunc (val, base) * base;
}
