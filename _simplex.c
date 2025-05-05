// Copyright (c) 2008, Casey Duncan (casey dot duncan at gmail dot com)
// see LICENSE.txt for details
// $Id$

#include "Python.h"
#include <math.h>
#include <float.h>
#include "_noise.h"
#include "_noise5d.h"

// 2D simplex skew factors
#define F2 0.3660254037844386f  // 0.5 * (sqrt(3.0) - 1.0)
#define G2 0.21132486540518713f // (3.0 - sqrt(3.0)) / 6.0

float 
noise2(float x, float y) 
{
	int i1, j1, I, J, c;
	float s = (x + y) * F2;
	float i = floorf(x + s);
	float j = floorf(y + s);
	float t = (i + j) * G2;

	float xx[3], yy[3], f[3];
	float noise[3] = {0.0f, 0.0f, 0.0f};
	int g[3];

	xx[0] = x - (i - t);
	yy[0] = y - (j - t);

	i1 = xx[0] > yy[0];
	j1 = xx[0] <= yy[0];

	xx[2] = xx[0] + G2 * 2.0f - 1.0f;
	yy[2] = yy[0] + G2 * 2.0f - 1.0f;
	xx[1] = xx[0] - i1 + G2;
	yy[1] = yy[0] - j1 + G2;

	I = (int) i & 255;
	J = (int) j & 255;
	g[0] = PERM[I + PERM[J]] % 12;
	g[1] = PERM[I + i1 + PERM[J + j1]] % 12;
	g[2] = PERM[I + 1 + PERM[J + 1]] % 12;

	for (c = 0; c <= 2; c++)
		f[c] = 0.5f - xx[c]*xx[c] - yy[c]*yy[c];
	
	for (c = 0; c <= 2; c++)
		if (f[c] > 0)
			noise[c] = f[c]*f[c]*f[c]*f[c] * (GRAD3[g[c]][0]*xx[c] + GRAD3[g[c]][1]*yy[c]);
	
	return (noise[0] + noise[1] + noise[2]) * 70.0f;
}

#define dot3(v1, v2) ((v1)[0]*(v2)[0] + (v1)[1]*(v2)[1] + (v1)[2]*(v2)[2])

#define ASSIGN(a, v0, v1, v2) (a)[0] = v0; (a)[1] = v1; (a)[2] = v2;

#define F3 (1.0f / 3.0f)
#define G3 (1.0f / 6.0f)

float 
noise3(float x, float y, float z) 
{
	int c, o1[3], o2[3], g[4], I, J, K;
	float f[4], noise[4] = {0.0f, 0.0f, 0.0f, 0.0f};
	float s = (x + y + z) * F3;
	float i = floorf(x + s);
	float j = floorf(y + s);
	float k = floorf(z + s);
	float t = (i + j + k) * G3;

	float pos[4][3];

	pos[0][0] = x - (i - t);
	pos[0][1] = y - (j - t);
	pos[0][2] = z - (k - t);

	if (pos[0][0] >= pos[0][1]) {
		if (pos[0][1] >= pos[0][2]) {
			ASSIGN(o1, 1, 0, 0);
			ASSIGN(o2, 1, 1, 0);
		} else if (pos[0][0] >= pos[0][2]) {
			ASSIGN(o1, 1, 0, 0);
			ASSIGN(o2, 1, 0, 1);
		} else {
			ASSIGN(o1, 0, 0, 1);
			ASSIGN(o2, 1, 0, 1);
		}
	} else {
		if (pos[0][1] < pos[0][2]) {
			ASSIGN(o1, 0, 0, 1);
			ASSIGN(o2, 0, 1, 1);
		} else if (pos[0][0] < pos[0][2]) {
			ASSIGN(o1, 0, 1, 0);
			ASSIGN(o2, 0, 1, 1);
		} else {
			ASSIGN(o1, 0, 1, 0);
			ASSIGN(o2, 1, 1, 0);
		}
	}
	
	for (c = 0; c <= 2; c++) {
		pos[3][c] = pos[0][c] - 1.0f + 3.0f * G3;
		pos[2][c] = pos[0][c] - o2[c] + 2.0f * G3;
		pos[1][c] = pos[0][c] - o1[c] + G3;
	}

	I = (int) i & 255; 
	J = (int) j & 255; 
	K = (int) k & 255;
	g[0] = PERM[I + PERM[J + PERM[K]]] % 12;
	g[1] = PERM[I + o1[0] + PERM[J + o1[1] + PERM[o1[2] + K]]] % 12;
	g[2] = PERM[I + o2[0] + PERM[J + o2[1] + PERM[o2[2] + K]]] % 12;
	g[3] = PERM[I + 1 + PERM[J + 1 + PERM[K + 1]]] % 12; 

	for (c = 0; c <= 3; c++) {
		f[c] = 0.6f - pos[c][0]*pos[c][0] - pos[c][1]*pos[c][1] - pos[c][2]*pos[c][2];
	}
	
	for (c = 0; c <= 3; c++) {
		if (f[c] > 0) {
			noise[c] = f[c]*f[c]*f[c]*f[c] * dot3(pos[c], GRAD3[g[c]]);
		}
	}
	
	return (noise[0] + noise[1] + noise[2] + noise[3]) * 32.0f;
}

static inline float
fbm_noise3(float x, float y, float z, int octaves, float persistence, float lacunarity) {
    float freq = 1.0f;
    float amp = 1.0f;
    float max = 1.0f;
    float total = noise3(x, y, z);
    int i;

    for (i = 1; i < octaves; ++i) {
        freq *= lacunarity;
        amp *= persistence;
        max += amp;
        total += noise3(x * freq, y * freq, z * freq) * amp;
    }
    return total / max;
}

#define dot4(v1, x, y, z, w) ((v1)[0]*(x) + (v1)[1]*(y) + (v1)[2]*(z) + (v1)[3]*(w))

#define F4 0.30901699437494745f /* (sqrt(5.0) - 1.0) / 4.0 */
#define G4 0.1381966011250105f /* (5.0 - sqrt(5.0)) / 20.0 */

float 
noise4(float x, float y, float z, float w) {
    float noise[5] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f};

    float s = (x + y + z + w) * F4;
    float i = floorf(x + s);
    float j = floorf(y + s);
    float k = floorf(z + s);
    float l = floorf(w + s);
    float t = (i + j + k + l) * G4;

    float x0 = x - (i - t);
    float y0 = y - (j - t);
    float z0 = z - (k - t);
    float w0 = w - (l - t);

    int c = (x0 > y0)*32 + (x0 > z0)*16 + (y0 > z0)*8 + (x0 > w0)*4 + (y0 > w0)*2 + (z0 > w0);
    int i1 = SIMPLEX[c][0]>=3;
    int j1 = SIMPLEX[c][1]>=3;
    int k1 = SIMPLEX[c][2]>=3;
    int l1 = SIMPLEX[c][3]>=3;
    int i2 = SIMPLEX[c][0]>=2;
    int j2 = SIMPLEX[c][1]>=2;
    int k2 = SIMPLEX[c][2]>=2;
    int l2 = SIMPLEX[c][3]>=2;
    int i3 = SIMPLEX[c][0]>=1;
    int j3 = SIMPLEX[c][1]>=1;
    int k3 = SIMPLEX[c][2]>=1;
    int l3 = SIMPLEX[c][3]>=1;

    float x1 = x0 - i1 + G4;
    float y1 = y0 - j1 + G4;
    float z1 = z0 - k1 + G4;
    float w1 = w0 - l1 + G4;
    float x2 = x0 - i2 + 2.0f*G4;
    float y2 = y0 - j2 + 2.0f*G4;
    float z2 = z0 - k2 + 2.0f*G4;
    float w2 = w0 - l2 + 2.0f*G4;
    float x3 = x0 - i3 + 3.0f*G4;
    float y3 = y0 - j3 + 3.0f*G4;
    float z3 = z0 - k3 + 3.0f*G4;
    float w3 = w0 - l3 + 3.0f*G4;
    float x4 = x0 - 1.0f + 4.0f*G4;
    float y4 = y0 - 1.0f + 4.0f*G4;
    float z4 = z0 - 1.0f + 4.0f*G4;
    float w4 = w0 - 1.0f + 4.0f*G4;

    int I = (int)i & 255;
    int J = (int)j & 255;
    int K = (int)k & 255;
    int L = (int)l & 255;
    int gi0 = PERM[I + PERM[J + PERM[K + PERM[L]]]] & 0x1f;
    int gi1 = PERM[I + i1 + PERM[J + j1 + PERM[K + k1 + PERM[L + l1]]]] & 0x1f; 
    int gi2 = PERM[I + i2 + PERM[J + j2 + PERM[K + k2 + PERM[L + l2]]]] & 0x1f; 
    int gi3 = PERM[I + i3 + PERM[J + j3 + PERM[K + k3 + PERM[L + l3]]]] & 0x1f; 
    int gi4 = PERM[I + 1 + PERM[J + 1 + PERM[K + 1 + PERM[L + 1]]]] & 0x1f;
    float t0, t1, t2, t3, t4;

    t0 = 0.6f - x0*x0 - y0*y0 - z0*z0 - w0*w0;
    if (t0 >= 0.0f) {
        t0 *= t0;
        noise[0] = t0 * t0 * dot4(GRAD4[gi0], x0, y0, z0, w0);
    }
    t1 = 0.6f - x1*x1 - y1*y1 - z1*z1 - w1*w1;
    if (t1 >= 0.0f) {
        t1 *= t1;
        noise[1] = t1 * t1 * dot4(GRAD4[gi1], x1, y1, z1, w1);
    }
    t2 = 0.6f - x2*x2 - y2*y2 - z2*z2 - w2*w2;
    if (t2 >= 0.0f) {
        t2 *= t2;
        noise[2] = t2 * t2 * dot4(GRAD4[gi2], x2, y2, z2, w2);
    }
    t3 = 0.6f - x3*x3 - y3*y3 - z3*z3 - w3*w3;
    if (t3 >= 0.0f) {
        t3 *= t3;
        noise[3] = t3 * t3 * dot4(GRAD4[gi3], x3, y3, z3, w3);
    }
    t4 = 0.6f - x4*x4 - y4*y4 - z4*z4 - w4*w4;
    if (t4 >= 0.0f) {
        t4 *= t4;
        noise[4] = t4 * t4 * dot4(GRAD4[gi4], x4, y4, z4, w4);
    }

    return 27.0 * (noise[0] + noise[1] + noise[2] + noise[3] + noise[4]);
}

static inline float
fbm_noise4(float x, float y, float z, float w, int octaves, float persistence, float lacunarity) {
    float freq = 1.0f;
    float amp = 1.0f;
    float max = 1.0f;
    float total = noise4(x, y, z, w);
    int i;

    for (i = 1; i < octaves; ++i) {
        freq *= lacunarity;
        amp *= persistence;
        max += amp;
        total += noise4(x * freq, y * freq, z * freq, w * freq) * amp;
    }
    return total / max;
}

#define dot5(v1, x, y, z, w, u) ((v1)[0]*(x) + (v1)[1]*(y) + (v1)[2]*(z) + (v1)[3]*(w) + (v1)[4]*(u))

#define F5 0.2898979485566356f /* (sqrt(6.0) - 1.0) / 5.0 */
#define G5 0.1183503419072274f /* (1.0 - 1.0/sqrt(6.0)) / 5.0 */

float 
noise5(float x, float y, float z, float w, float u) {
    float noise[6] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};

    float s = (x + y + z + w + u) * F5;
    float i = floorf(x + s);
    float j = floorf(y + s);
    float k = floorf(z + s);
    float l = floorf(w + s);
    float m = floorf(u + s);
    float t = (i + j + k + l + m) * G5;

    float x0 = x - (i - t);
    float y0 = y - (j - t);
    float z0 = z - (k - t);
    float w0 = w - (l - t);
    float u0 = u - (m - t);

    // int c = (x0 > y0)*32 + (x0 > z0)*16 + (y0 > z0)*8 + (x0 > w0)*4 + (y0 > w0)*2 + (z0 > w0);
    int c = (x0 > y0)*512 + (x0 > z0)*256 + (x0 > w0)*128 + (x0 > u0)*64 + (y0 > z0)*32 + (y0 > w0)*16 + (y0 > u0)*8 + (z0 > w0)*4 + (z0 > u0)*2 + (w0 > u0);
    int i1 = SIMPLEX5[c][0]>=4;
    int j1 = SIMPLEX5[c][1]>=4;
    int k1 = SIMPLEX5[c][2]>=4;
    int l1 = SIMPLEX5[c][3]>=4;
    int m1 = SIMPLEX5[c][4]>=4;
    int i2 = SIMPLEX5[c][0]>=3;
    int j2 = SIMPLEX5[c][1]>=3;
    int k2 = SIMPLEX5[c][2]>=3;
    int l2 = SIMPLEX5[c][3]>=3;
    int m2 = SIMPLEX5[c][4]>=3;
    int i3 = SIMPLEX5[c][0]>=2;
    int j3 = SIMPLEX5[c][1]>=2;
    int k3 = SIMPLEX5[c][2]>=2;
    int l3 = SIMPLEX5[c][3]>=2;
    int m3 = SIMPLEX5[c][4]>=2;
    int i4 = SIMPLEX5[c][0]>=1;
    int j4 = SIMPLEX5[c][1]>=1;
    int k4 = SIMPLEX5[c][2]>=1;
    int l4 = SIMPLEX5[c][3]>=1;
    int m4 = SIMPLEX5[c][4]>=1;

    // the formulas below need to be checked for correctness when using 5D simplex traversal lookup table; the formulas above might be edited instead
    float x1 = x0 - i1 + G5;
    float y1 = y0 - j1 + G5;
    float z1 = z0 - k1 + G5;
    float w1 = w0 - l1 + G5;
    float u1 = u0 - m1 + G5;
    float x2 = x0 - i2 + 2.0f*G5;
    float y2 = y0 - j2 + 2.0f*G5;
    float z2 = z0 - k2 + 2.0f*G5;
    float w2 = w0 - l2 + 2.0f*G5;
    float u2 = u0 - m2 + 2.0f*G5;
    float x3 = x0 - i3 + 3.0f*G5;
    float y3 = y0 - j3 + 3.0f*G5;
    float z3 = z0 - k3 + 3.0f*G5;
    float w3 = w0 - l3 + 3.0f*G5;
    float u3 = u0 - m3 + 3.0f*G5;
    float x4 = x0 - i4 + 4.0f*G5;
    float y4 = y0 - j4 + 4.0f*G5;
    float z4 = z0 - k4 + 4.0f*G5;
    float w4 = w0 - l4 + 4.0f*G5;
    float u4 = u0 - m4 + 4.0f*G5;
    float x5 = x0 - 1.0f + 5.0f*G5;
    float y5 = y0 - 1.0f + 5.0f*G5;
    float z5 = z0 - 1.0f + 5.0f*G5;
    float w5 = w0 - 1.0f + 5.0f*G5;
    float u5 = u0 - 1.0f + 5.0f*G5;

    int I = (int)i & 1023;
    int J = (int)j & 1023;
    int K = (int)k & 1023;
    int L = (int)l & 1023;
    int M = (int)m & 1023;

    int gi0 = PERM5[I + PERM5[J + PERM5[K + PERM5[L + PERM5[M]]]]] & 0x1f;
    int gi1 = PERM5[I + i1 + PERM5[J + j1 + PERM5[K + k1 + PERM5[L + l1 + PERM5[M + m1]]]]] & 0x1f; 
    int gi2 = PERM5[I + i2 + PERM5[J + j2 + PERM5[K + k2 + PERM5[L + l2 + PERM5[M + m2]]]]] & 0x1f; 
    int gi3 = PERM5[I + i3 + PERM5[J + j3 + PERM5[K + k3 + PERM5[L + l3 + PERM5[M + m3]]]]] & 0x1f; 
    int gi4 = PERM5[I + i4 + PERM5[J + j4 + PERM5[K + k4 + PERM5[L + l4 + PERM5[M + m4]]]]] & 0x1f; 
    int gi5 = PERM5[I + 1 + PERM5[J + 1 + PERM5[K + 1 + PERM5[L + 1 + PERM5[M + 1]]]]] & 0x1f;;
    float t0, t1, t2, t3, t4, t5;

    t0 = 0.7f - x0*x0 - y0*y0 - z0*z0 - w0*w0 - u0*u0;
    if (t0 >= 0.0f) {
        t0 *= t0;
        noise[0] = t0 * t0 * dot5(GRAD5[gi0], x0, y0, z0, w0, u0);
    }
    t1 = 0.7f - x1*x1 - y1*y1 - z1*z1 - w1*w1 - u1*u1;
    if (t1 >= 0.0f) {
        t1 *= t1;
        noise[1] = t1 * t1 * dot5(GRAD5[gi1], x1, y1, z1, w1, u1);
    }
    t2 = 0.7f - x2*x2 - y2*y2 - z2*z2 - w2*w2 - u2*u2;
    if (t2 >= 0.0f) {
        t2 *= t2;
        noise[2] = t2 * t2 * dot5(GRAD5[gi2], x2, y2, z2, w2, u2);
    }
    t3 = 0.7f - x3*x3 - y3*y3 - z3*z3 - w3*w3 - u3*u3;
    if (t3 >= 0.0f) {
        t3 *= t3;
        noise[3] = t3 * t3 * dot5(GRAD5[gi3], x3, y3, z3, w3, u3);
    }
    t4 = 0.7f - x4*x4 - y4*y4 - z4*z4 - w4*w4 - u4*u4;
    if (t4 >= 0.0f) {
        t4 *= t4;
        noise[4] = t4 * t4 * dot5(GRAD5[gi4], x4, y4, z4, w4, u4);
    }
    t5 = 0.7f - x5*x5 - y5*y5 - z5*z5 - w5*w5 - u5*u5;
    if (t5 >= 0.0f) {
        t5 *= t5;
        noise[5] = t5 * t5 * dot5(GRAD5[gi5], x5, y5, z5, w5, u5);
    }
// change the numeric factor in return
    return 10 * (noise[0] + noise[1] + noise[2] + noise[3] + noise[4] + noise[5]);
}

static inline float
fbm_noise5(float x, float y, float z, float w, float u, int octaves, float persistence, float lacunarity) {
    float freq = 1.0f;
    float amp = 1.0f;
    float max = 1.0f;
    float total = noise5(x, y, z, w, u);
    int i;

    for (i = 1; i < octaves; ++i) {
        freq *= lacunarity;
        amp *= persistence;
        max += amp;
        total += noise5(x * freq, y * freq, z * freq, w * freq, u * freq) * amp;
    }
    return total / max;
}

static PyObject *
py_noise2(PyObject *self, PyObject *args, PyObject *kwargs)
{
	float x, y;
	int octaves = 1;
	float persistence = 0.5f;
    float lacunarity = 2.0f;
    float repeatx = FLT_MAX;
    float repeaty = FLT_MAX;
    float z = 0.0f;
	static char *kwlist[] = {"x", "y", "octaves", "persistence", "lacunarity", 
        "repeatx", "repeaty", "base", NULL};

	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ff|ifffff:snoise2", kwlist,
		&x, &y, &octaves, &persistence, &lacunarity, &repeatx, &repeaty, &z)) {
		return NULL;
    }
    if (octaves <= 0) {
        PyErr_SetString(PyExc_ValueError, "Expected octaves value > 0");
        return NULL;
    }
	
    if (repeatx == FLT_MAX && repeaty == FLT_MAX) {
        // Flat noise, no tiling
        float freq = 1.0f;
        float amp = 1.0f;
        float max = 1.0f;
        float total = noise2(x + z, y + z);
        int i;

        for (i = 1; i < octaves; i++) {
            freq *= lacunarity;
            amp *= persistence;
            max += amp;
            total += noise2(x * freq + z, y * freq + z) * amp;
        }
        return (PyObject *) PyFloat_FromDouble((double) (total / max));
    } else { // Tiled noise
        float w = z;
        if (repeaty != FLT_MAX) {
            float yf = y * 2.0 / repeaty;
            float yr = repeaty * M_1_PI * 0.5;
            float vy = fast_sin(yf);
            float vyz = fast_cos(yf);
            y = vy * yr;
            w += vyz * yr;
            if (repeatx == FLT_MAX) {
                return (PyObject *) PyFloat_FromDouble(
                    (double) fbm_noise3(x, y, w, octaves, persistence, lacunarity));
            }
        }
        if (repeatx != FLT_MAX) {
            float xf = x * 2.0 / repeatx;
            float xr = repeatx * M_1_PI * 0.5;
            float vx = fast_sin(xf);
            float vxz = fast_cos(xf);
            x = vx * xr;
            z += vxz * xr;
            if (repeaty == FLT_MAX) {
                return (PyObject *) PyFloat_FromDouble(
                    (double) fbm_noise3(x, y, z, octaves, persistence, lacunarity));
            }
        }
        return (PyObject *) PyFloat_FromDouble(
            (double) fbm_noise4(x, y, z, w, octaves, persistence, lacunarity));
    }
}

static PyObject *
py_noise3(PyObject *self, PyObject *args, PyObject *kwargs)
{
	float x, y, z;
	int octaves = 1;
	float persistence = 0.5f;
    float lacunarity = 2.0f;

	static char *kwlist[] = {"x", "y", "z", "octaves", "persistence", "lacunarity", NULL};

	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "fff|iff:snoise3", kwlist,
		&x, &y, &z, &octaves, &persistence, &lacunarity))
		return NULL;
	
	if (octaves == 1) {
		// Single octave, return simple noise
		return (PyObject *) PyFloat_FromDouble((double) noise3(x, y, z));
	} else if (octaves > 1) {
		return (PyObject *) PyFloat_FromDouble(
            (double) fbm_noise3(x, y, z, octaves, persistence, lacunarity));
	} else {
		PyErr_SetString(PyExc_ValueError, "Expected octaves value > 0");
		return NULL;
	}
}

static PyObject *
py_noise4(PyObject *self, PyObject *args, PyObject *kwargs)
{
	float x, y, z, w;
	int octaves = 1;
	float persistence = 0.5f;
	float lacunarity = 2.0f;

	static char *kwlist[] = {"x", "y", "z", "w", "octaves", "persistence", "lacunarity", NULL};

	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ffff|iff:snoise4", kwlist,
		&x, &y, &z, &w, &octaves, &persistence, &lacunarity))
		return NULL;
	
	if (octaves == 1) {
		// Single octave, return simple noise
		return (PyObject *) PyFloat_FromDouble((double) noise4(x, y, z, w));
	} else if (octaves > 1) {
		return (PyObject *) PyFloat_FromDouble(
            (double) fbm_noise4(x, y, z, w, octaves, persistence, lacunarity));
	} else {
		PyErr_SetString(PyExc_ValueError, "Expected octaves value > 0");
		return NULL;
	}
}

static PyObject *
py_noise5(PyObject *self, PyObject *args, PyObject *kwargs)
{
	float x, y, z, w, u;
	int octaves = 1;
	float persistence = 0.5f;
	float lacunarity = 2.0f;

	static char *kwlist[] = {"x", "y", "z", "w", "u", "octaves", "persistence", "lacunarity", NULL};

	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ffff|iff:snoise4", kwlist,
		&x, &y, &z, &w, &u, &octaves, &persistence, &lacunarity))
		return NULL;
	
	if (octaves == 1) {
		// Single octave, return simple noise
		return (PyObject *) PyFloat_FromDouble((double) noise5(x, y, z, w, u));
	} else if (octaves > 1) {
		return (PyObject *) PyFloat_FromDouble(
            (double) fbm_noise5(x, y, z, w, u, octaves, persistence, lacunarity));
	} else {
		PyErr_SetString(PyExc_ValueError, "Expected octaves value > 0");
		return NULL;
	}
}

static PyMethodDef simplex_functions[] = {
	{"noise2", (PyCFunction)py_noise2, METH_VARARGS | METH_KEYWORDS, 
		"noise2(x, y, octaves=1, persistence=0.5, lacunarity=2.0, repeatx=None, repeaty=None, base=0.0) "
        "return simplex noise value for specified 2D coordinate.\n\n"
		"octaves -- specifies the number of passes, defaults to 1 (simple noise).\n\n"
		"persistence -- specifies the amplitude of each successive octave relative\n"
		"to the one below it. Defaults to 0.5 (each higher octave's amplitude\n"
		"is halved). Note the amplitude of the first pass is always 1.0.\n\n"
        "lacunarity -- specifies the frequency of each successive octave relative\n"
        "to the one below it, similar to persistence. Defaults to 2.0.\n\n"
        "repeatx, repeaty -- specifies the interval along each axis when \n"
		"the noise values repeat. This can be used as the tile size for creating \n"
		"tileable textures\n\n"
		"base -- specifies a fixed offset for the noise coordinates. Useful for\n"
		"generating different noise textures with the same repeat interval"},
	{"noise3", (PyCFunction)py_noise3, METH_VARARGS | METH_KEYWORDS, 
		"noise3(x, y, z, octaves=1, persistence=0.5, lacunarity=2.0) return simplex noise value for "
		"specified 3D coordinate\n\n"
		"octaves -- specifies the number of passes, defaults to 1 (simple noise).\n\n"
		"persistence -- specifies the amplitude of each successive octave relative\n"
		"to the one below it. Defaults to 0.5 (each higher octave's amplitude\n"
		"is halved). Note the amplitude of the first pass is always 1.0.\n\n"
        "lacunarity -- specifies the frequency of each successive octave relative\n"
        "to the one below it, similar to persistence. Defaults to 2.0."},
	{"noise4", (PyCFunction)py_noise4, METH_VARARGS | METH_KEYWORDS, 
		"noise4(x, y, z, w, octaves=1, persistence=0.5, lacunarity=2.0) return simplex noise value for "
		"specified 4D coordinate\n\n"
		"octaves -- specifies the number of passes, defaults to 1 (simple noise).\n\n"
		"persistence -- specifies the amplitude of each successive octave relative\n"
		"to the one below it. Defaults to 0.5 (each higher octave's amplitude\n"
		"is halved). Note the amplitude of the first pass is always 1.0.\n\n"
        "lacunarity -- specifies the frequency of each successive octave relative\n"
        "to the one below it, similar to persistence. Defaults to 2.0."},
    {"noise5", (PyCFunction)py_noise4, METH_VARARGS | METH_KEYWORDS, 
        "noise5(x, y, z, w, u, octaves=1, persistence=0.5, lacunarity=2.0) return simplex noise value for "
        "specified 5D coordinate\n\n"
        "octaves -- specifies the number of passes, defaults to 1 (simple noise).\n\n"
        "persistence -- specifies the amplitude of each successive octave relative\n"
        "to the one below it. Defaults to 0.5 (each higher octave's amplitude\n"
        "is halved). Note the amplitude of the first pass is always 1.0.\n\n"
        "lacunarity -- specifies the frequency of each successive octave relative\n"
        "to the one below it, similar to persistence. Defaults to 2.0."},
	{NULL}
};

PyDoc_STRVAR(module_doc, "Native-code simplex noise functions");

#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef moduledef = {
	PyModuleDef_HEAD_INIT,
	"_simplex",
	module_doc,
	-1,                 /* m_size */
	simplex_functions,  /* m_methods */
	NULL,               /* m_reload (unused) */
	NULL,               /* m_traverse */
	NULL,               /* m_clear */
	NULL                /* m_free */
};

PyObject *
PyInit__simplex(void)
{
    return PyModule_Create(&moduledef);
}

#else

void
init_simplex(void)
{
	Py_InitModule3("_simplex", simplex_functions, module_doc);
}

#endif
