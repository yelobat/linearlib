/**
 * Author: wwotz
 * Add #LINEARLIB_IMPLEMENTATION to the start of the implementation
 * file in order to add the implementation code to your project. 
 */

#ifndef LINEARLIB_H_
#define LINEARLIB_H_

//requires -lm flag to link the math library
#include <math.h>

#ifndef LINEARLIBDEF
#ifdef LINEARLIBSTATIC
#define LINEARLIBDEF static
#else /* !defined(LINEARLIBSTATIC) */
#define LINEARLIBDEF extern
#endif /* LINEARLIBSTATIC */
#endif /* LINEARLIBDEF */

/* &optional */
#define LL_USE_MATRIX

typedef enum matrix_type_t {
        LL_MATRIX_MODEL,
        LL_MATRIX_VIEW,
        LL_MATRIX_PROJECTION,
        LL_MATRIX_COUNT
} matrix_type_t;

typedef union vec2_t {
        float data[2];
        struct {
                float x;
                float y;
        };
        struct {
                float r;
                float g;
        };
        struct {
                float s;
                float t;
        };
} vec2_t;

typedef union vec3_t {
        float data[3];
        struct {
                float x;
                float y;
                float z;
        };
        struct {
                float r;
                float g;
                float b;
        };        
} vec3_t;

typedef union vec4_t {
        float data[4];
        struct {
                float x;
                float y;
                float z;
                float w;
        };
        struct {
                float r;
                float g;
                float b;
                float a;
        };
} vec4_t;

typedef union ivec2_t {
        int data[2];
        struct {
                int x;
                int y;
        };
        struct {
                int r;
                int g;
        };
        struct {
                int s;
                int t;
        };
} ivec2_t;

typedef union ivec3_t {
        int data[3];
        struct {
                int x;
                int y;
                int z;
        };
        struct {
                int r;
                int g;
                int b;
        };
} ivec3_t;

typedef union ivec4_t {
        int data[4];
        struct {
                int x;
                int y;
                int z;
                int w;
        };
        struct {
                int r;
                int g;
                int b;
                int a;
        };
} ivec4_t;

typedef union mat2_t {
        float data[4];
        struct {
                float m00, m01;
                float m10, m11;
        };
} mat2_t;

typedef union mat3_t {
        float data[9];
        struct {
                float m00, m01, m02;
                float m10, m11, m12;
                float m20, m21, m22;
        };
} mat3_t;

typedef union mat4_t {
        float data[16];
        struct {
                float m00, m01, m02, m03;
                float m10, m11, m12, m13;
                float m20, m21, m22, m23;
                float m30, m31, m32, m33;
        };
} mat4_t;

LINEARLIBDEF vec2_t
ll_vec2_create2f(float x, float y);
LINEARLIBDEF vec2_t
ll_vec2_create2fv(vec2_t vec);
LINEARLIBDEF float
ll_vec2_length2fv(vec2_t vec);
LINEARLIBDEF float
ll_vec2_length2f(float x, float y);
LINEARLIBDEF float
ll_vec2_length_squared2fv(vec2_t vec);
LINEARLIBDEF float
ll_vec2_length_squared2f(float x, float y);
LINEARLIBDEF vec2_t
ll_vec2_add2fv(vec2_t left, vec2_t right);
LINEARLIBDEF vec2_t
ll_vec2_add2f(vec2_t left, float x, float y);
LINEARLIBDEF vec2_t
ll_vec2_add1f(vec2_t left, float value);
LINEARLIBDEF vec2_t
ll_vec2_sub2fv(vec2_t left, vec2_t right);
LINEARLIBDEF vec2_t
ll_vec2_sub2f(vec2_t left, float x, float y);
LINEARLIBDEF vec2_t
ll_vec2_sub1f(vec2_t left, float value);
LINEARLIBDEF vec2_t
ll_vec2_mul2fv(vec2_t left, vec2_t right);
LINEARLIBDEF vec2_t
ll_vec2_mul2f(vec2_t left, float x, float y);
LINEARLIBDEF vec2_t
ll_vec2_mul1f(vec2_t left, float value);
LINEARLIBDEF vec2_t
ll_vec2_div2fv(vec2_t left, vec2_t right);
LINEARLIBDEF vec2_t
ll_vec2_div2f(vec2_t left, float x, float y);
LINEARLIBDEF vec2_t
ll_vec2_div1f(vec2_t left, float value);
LINEARLIBDEF float
ll_vec2_dot2fv(vec2_t left, vec2_t right);
LINEARLIBDEF float
ll_vec2_dot2f(vec2_t left, float x, float y);
LINEARLIBDEF float
ll_vec2_cross2fv(vec2_t left, vec2_t right);
LINEARLIBDEF float
ll_vec2_cross2f(vec2_t left, float x, float y);
LINEARLIBDEF vec2_t
ll_vec2_normalise2fv(vec2_t vec);
LINEARLIBDEF vec2_t
ll_vec2_normalise2f(float x, float y);

LINEARLIBDEF vec3_t
ll_vec3_create3f(float x, float y, float z);
LINEARLIBDEF vec3_t
ll_vec3_create3fv(vec3_t ivec);
LINEARLIBDEF float
ll_vec3_length3fv(vec3_t ivec);
LINEARLIBDEF float
ll_vec3_length3f(float x, float y, float z);
LINEARLIBDEF float
ll_vec3_length_squared3fv(vec3_t ivec);
LINEARLIBDEF float
ll_vec3_length_squared3f(float x, float y, float z);
LINEARLIBDEF vec3_t
ll_vec3_add3fv(vec3_t left, vec3_t right);
LINEARLIBDEF vec3_t
ll_vec3_add3f(vec3_t left, float x, float y, float z);
LINEARLIBDEF vec3_t
ll_vec3_add1f(vec3_t left, float value);
LINEARLIBDEF vec3_t
ll_vec3_sub3fv(vec3_t left, vec3_t right);
LINEARLIBDEF vec3_t
ll_vec3_sub3f(vec3_t left, float x, float y, float z);
LINEARLIBDEF vec3_t
ll_vec3_sub1f(vec3_t left, float value);
LINEARLIBDEF vec3_t
ll_vec3_mul3fv(vec3_t left, vec3_t right);
LINEARLIBDEF vec3_t
ll_vec3_mul3f(vec3_t left, float x, float y, float z);
LINEARLIBDEF vec3_t
ll_vec3_mul1f(vec3_t left, float value);
LINEARLIBDEF vec3_t
ll_vec3_div3fv(vec3_t left, vec3_t right);
LINEARLIBDEF vec3_t
ll_vec3_div3f(vec3_t left, float x, float y, float z);
LINEARLIBDEF vec3_t
ll_vec3_div1f(vec3_t left, float value);
LINEARLIBDEF float
ll_vec3_dot3fv(vec3_t left, vec3_t right);
LINEARLIBDEF float
ll_vec3_dot3f(vec3_t left, float x, float y, float z);
LINEARLIBDEF vec3_t
ll_vec3_cross3fv(vec3_t left, vec3_t right);
LINEARLIBDEF vec3_t
ll_vec3_cross3f(vec3_t left, float x, float y, float z);
LINEARLIBDEF vec3_t
ll_vec3_normalise3fv(vec3_t ivec);
LINEARLIBDEF vec3_t
ll_vec3_normalise3f(float x, float y, float z);

LINEARLIBDEF vec4_t
ll_vec4_create4f(float x, float y, float z, float w);
LINEARLIBDEF vec4_t
ll_vec4_create4fv(vec4_t ivec);
LINEARLIBDEF float
ll_vec4_length4fv(vec4_t ivec);
LINEARLIBDEF float
ll_vec4_length4f(float x, float y, float z, float w);
LINEARLIBDEF float
ll_vec4_length_squared4fv(vec4_t ivec);
LINEARLIBDEF float
ll_vec4_length_squared4f(float x, float y, float z, float w);
LINEARLIBDEF vec4_t
ll_vec4_add4fv(vec4_t left, vec4_t right);
LINEARLIBDEF vec4_t
ll_vec4_add4f(vec4_t left, float x, float y, float z, float w);
LINEARLIBDEF vec4_t
ll_vec4_add1f(vec4_t left, float value);
LINEARLIBDEF vec4_t
ll_vec4_sub4fv(vec4_t left, vec4_t right);
LINEARLIBDEF vec4_t
ll_vec4_sub4f(vec4_t left, float x, float y, float z, float w);
LINEARLIBDEF vec4_t
ll_vec4_sub1f(vec4_t left, float value);
LINEARLIBDEF vec4_t
ll_vec4_mul4fv(vec4_t left, vec4_t right);
LINEARLIBDEF vec4_t
ll_vec4_mul4f(vec4_t left, float x, float y, float z, float w);
LINEARLIBDEF vec4_t
ll_vec4_mul1f(vec4_t left, float value);
LINEARLIBDEF vec4_t
ll_vec4_div4fv(vec4_t left, vec4_t right);
LINEARLIBDEF vec4_t
ll_vec4_div4f(vec4_t left, float x, float y, float z, float w);
LINEARLIBDEF vec4_t
ll_vec4_div1f(vec4_t left, float value);
LINEARLIBDEF float
ll_vec4_dot4fv(vec4_t left, vec4_t right);
LINEARLIBDEF float
ll_vec4_dot4f(vec4_t left, float x, float y, float z, float w);
LINEARLIBDEF vec4_t
ll_vec4_normalise4fv(vec4_t ivec);
LINEARLIBDEF vec4_t
ll_vec4_normalise4f(float x, float y, float z, float w);

LINEARLIBDEF ivec2_t
ll_ivec2_create2i(int x, int y);
LINEARLIBDEF ivec2_t
ll_ivec2_create2iv(ivec2_t ivec);
LINEARLIBDEF float
ll_ivec2_length2iv(ivec2_t ivec);
LINEARLIBDEF float
ll_ivec2_length2i(int x, int y);
LINEARLIBDEF ivec2_t
ll_ivec2_add2iv(ivec2_t left, ivec2_t right);
LINEARLIBDEF ivec2_t
ll_ivec2_add2i(ivec2_t left, int x, int y);
LINEARLIBDEF ivec2_t
ll_ivec2_add1i(ivec2_t left, int value);
LINEARLIBDEF ivec2_t
ll_ivec2_sub2iv(ivec2_t left, ivec2_t right);
LINEARLIBDEF ivec2_t
ll_ivec2_sub2i(ivec2_t left, int x, int y);
LINEARLIBDEF ivec2_t
ll_ivec2_sub1i(ivec2_t left, int value);
LINEARLIBDEF ivec2_t
ll_ivec2_mul2iv(ivec2_t left, ivec2_t right);
LINEARLIBDEF ivec2_t
ll_ivec2_mul2i(ivec2_t left, int x, int y);
LINEARLIBDEF ivec2_t
ll_ivec2_mul1i(ivec2_t left, int value);
LINEARLIBDEF ivec2_t
ll_ivec2_div2iv(ivec2_t left, ivec2_t right);
LINEARLIBDEF ivec2_t
ll_ivec2_div2i(ivec2_t left, int x, int y);
LINEARLIBDEF ivec2_t
ll_ivec2_div1i(ivec2_t left, int value);
LINEARLIBDEF float
ll_ivec2_dot2iv(ivec2_t left, ivec2_t right);
LINEARLIBDEF float
ll_ivec2_dot2i(ivec2_t left, int x, int y);
LINEARLIBDEF float
ll_ivec2_cross2iv(ivec2_t left, ivec2_t right);
LINEARLIBDEF float
ll_ivec2_cross2i(ivec2_t left, int x, int y);
LINEARLIBDEF ivec2_t
ll_ivec2_normalise2iv(ivec2_t ivec);
LINEARLIBDEF ivec2_t
ll_ivec2_normalise2i(int x, int y);

LINEARLIBDEF ivec3_t
ll_ivec3_create3i(int x, int y, int z);
LINEARLIBDEF ivec3_t
ll_ivec3_create3iv(ivec3_t ivec);
LINEARLIBDEF float
ll_ivec3_length3iv(ivec3_t ivec);
LINEARLIBDEF float
ll_ivec3_length3i(int x, int y, int z);
LINEARLIBDEF float
ll_ivec3_length_squared3iv(ivec3_t ivec);
LINEARLIBDEF float
ll_ivec3_length_squared3i(int x, int y, int z);
LINEARLIBDEF ivec3_t
ll_ivec3_add3iv(ivec3_t left, ivec3_t right);
LINEARLIBDEF ivec3_t
ll_ivec3_add3i(ivec3_t left, int x, int y, int z);
LINEARLIBDEF ivec3_t
ll_ivec3_add1i(ivec3_t left, int value);
LINEARLIBDEF ivec3_t
ll_ivec3_sub3iv(ivec3_t left, ivec3_t right);
LINEARLIBDEF ivec3_t
ll_ivec3_sub3i(ivec3_t left, int x, int y, int z);
LINEARLIBDEF ivec3_t
ll_ivec3_sub1i(ivec3_t left, int value);
LINEARLIBDEF ivec3_t
ll_ivec3_mul3iv(ivec3_t left, ivec3_t right);
LINEARLIBDEF ivec3_t
ll_ivec3_mul3i(ivec3_t left, int x, int y, int z);
LINEARLIBDEF ivec3_t
ll_ivec3_mul1i(ivec3_t left, int value);
LINEARLIBDEF ivec3_t
ll_ivec3_div3iv(ivec3_t left, ivec3_t right);
LINEARLIBDEF ivec3_t
ll_ivec3_div3i(ivec3_t left, int x, int y, int z);
LINEARLIBDEF ivec3_t
ll_ivec3_div1i(ivec3_t left, int value);
LINEARLIBDEF float
ll_ivec3_dot3iv(ivec3_t left, ivec3_t right);
LINEARLIBDEF float
ll_ivec3_dot3i(ivec3_t left, int x, int y, int z);
LINEARLIBDEF ivec3_t
ll_ivec3_cross3iv(ivec3_t left, ivec3_t right);
LINEARLIBDEF ivec3_t
ll_ivec3_cross3i(ivec3_t left, int x, int y, int z);
LINEARLIBDEF ivec3_t
ll_ivec3_normalise3iv(ivec3_t ivec);
LINEARLIBDEF ivec3_t
ll_ivec3_normalise3i(int x, int y, int z);

LINEARLIBDEF ivec4_t
ll_ivec4_create4i(int x, int y, int z, int w);
LINEARLIBDEF ivec4_t
ll_ivec4_create4iv(ivec4_t ivec);
LINEARLIBDEF float
ll_ivec4_length4iv(ivec4_t ivec);
LINEARLIBDEF float
ll_ivec4_length4i(int x, int y, int z, int w);
LINEARLIBDEF float
ll_ivec4_length_squared4iv(ivec4_t ivec);
LINEARLIBDEF float
ll_ivec4_length_squared4i(int x, int y, int z, int w);
LINEARLIBDEF ivec4_t
ll_ivec4_add4iv(ivec4_t left, ivec4_t right);
LINEARLIBDEF ivec4_t
ll_ivec4_add4i(ivec4_t left, int x, int y, int z, int w);
LINEARLIBDEF ivec4_t
ll_ivec4_add1i(ivec4_t left, int value);
LINEARLIBDEF ivec4_t
ll_ivec4_sub4iv(ivec4_t left, ivec4_t right);
LINEARLIBDEF ivec4_t
ll_ivec4_sub4i(ivec4_t left, int x, int y, int z, int w);
LINEARLIBDEF ivec4_t
ll_ivec4_sub1i(ivec4_t left, int value);
LINEARLIBDEF ivec4_t
ll_ivec4_mul4iv(ivec4_t left, ivec4_t right);
LINEARLIBDEF ivec4_t
ll_ivec4_mul4i(ivec4_t left, int x, int y, int z, int w);
LINEARLIBDEF ivec4_t
ll_ivec4_mul1i(ivec4_t left, int value);
LINEARLIBDEF ivec4_t
ll_ivec4_div4iv(ivec4_t left, ivec4_t right);
LINEARLIBDEF ivec4_t
ll_ivec4_div4i(ivec4_t left, int x, int y, int z, int w);
LINEARLIBDEF ivec4_t
ll_ivec4_div1i(ivec4_t left, int value);
LINEARLIBDEF float
ll_ivec4_dot4iv(ivec4_t left, ivec4_t right);
LINEARLIBDEF float
ll_ivec4_dot4i(ivec4_t left, int x, int y, int z, int w);
LINEARLIBDEF ivec4_t
ll_ivec4_normalise4iv(ivec4_t ivec);
LINEARLIBDEF ivec4_t
ll_ivec4_normalise4i(int x, int y, int z, int w);

LINEARLIBDEF void
ll_mat4_multiply(mat4_t *left, mat4_t *right);
LINEARLIBDEF void
ll_mat4_copy(mat4_t *to, mat4_t *from);
LINEARLIBDEF void
ll_mat4_identity(mat4_t *mat);
LINEARLIBDEF void
ll_mat4_translate3f(mat4_t *mat, float dx, float dy, float dz);
LINEARLIBDEF void
ll_mat4_translate3fv(mat4_t *mat, vec3_t vec);
LINEARLIBDEF void
ll_mat4_scale3f(mat4_t *mat, float w, float h, float d);
LINEARLIBDEF void
ll_mat4_scale3fv(mat4_t *mat, vec3_t vec);
LINEARLIBDEF void
ll_mat4_rotate3f(mat4_t *mat, float x, float y, float z, float angle);
LINEARLIBDEF void
ll_mat4_rotate3fv(mat4_t *mat, vec3_t vec, float angle);
LINEARLIBDEF void
ll_mat4_orthographic(mat4_t *mat, float top, float right,
                     float bottom, float left, float near, float far);
LINEARLIBDEF void
ll_mat4_perspective(mat4_t *mat, float fovy, float aspect,
                    float near, float far);
LINEARLIBDEF void
ll_mat4_frustum(mat4_t *mat, float left, float right,
                float bottom, float top, float near, float far);

#ifdef LL_USE_MATRIX
LINEARLIBDEF void
ll_matrix_mode(matrix_type_t type);
LINEARLIBDEF void
ll_matrix_multiply(mat4_t *right);
LINEARLIBDEF void
ll_matrix_identity(void);
LINEARLIBDEF void
ll_matrix_translate3f(float dx, float dy, float dz);
LINEARLIBDEF void
ll_matrix_translate3fv(vec3_t vec);
LINEARLIBDEF void
ll_matrix_scale3f(float w, float h, float d);
LINEARLIBDEF void
ll_matrix_scale3fv(vec3_t vec);
LINEARLIBDEF void
ll_matrix_rotate3f(float x, float y, float z, float angle);
LINEARLIBDEF void
ll_matrix_rotate3fv(vec3_t vec, float angle);
LINEARLIBDEF void
ll_matrix_orthographic(float top, float right,
                       float bottom, float left, float near, float far);
LINEARLIBDEF void
ll_matrix_perspective(float fovy, float aspect,
                      float near, float far);
LINEARLIBDEF void
ll_matrix_frustum(float left, float right,
                  float bottom, float top, float near, float far);

LINEARLIBDEF mat4_t
ll_matrix_get_copy(void);

#endif /* LL_USE_MATRIX */

#ifdef LINEARLIB_IMPLEMENTATION

/**
 * @return A vec2_t structure.
 */
LINEARLIBDEF vec2_t
ll_vec2_create2f(float x, float y)
{
	return (vec2_t) { x, y };
}

/**
 * @description Essentially creates a new copy of @vec.
 *
 * @return A vec2_t structure.
 */
LINEARLIBDEF vec2_t
ll_vec2_create2fv(vec2_t vec)
{
	return vec;
}

/**
 * @return The length of @vec.
 */
LINEARLIBDEF float
ll_vec2_length2fv(vec2_t vec)
{
	return sqrtf((vec.x * vec.x) + (vec.y * vec.y));
}

/**
 * @return The length of a vector, represented by { x, y }.
 */
LINEARLIBDEF float
ll_vec2_length2f(float x, float y)
{
	return sqrtf((x * x) + (y * y));
}

/**
 * @return The squared length of @vec.
 */
LINEARLIBDEF float
ll_vec2_length_squared2fv(vec2_t vec)
{
	float length = ll_vec2_length2fv(vec);
	return length*length;
}

/**
 * @return The squared length of vector, represented by { x, y }.
 */
LINEARLIBDEF float
ll_vec2_length_squared2f(float x, float y)
{
	float length = ll_vec2_length2f(x, y);
	return length*length;
}

/**
 * @return A new vector, that is the sum of @left and @right.
 */
LINEARLIBDEF vec2_t
ll_vec2_add2fv(vec2_t left, vec2_t right)
{
	return ll_vec2_create2f(left.x + right.x, left.y + right.y);
}

/**
 * @return A new vector, that is the sum of @left and { x, y }.
 */
LINEARLIBDEF vec2_t
ll_vec2_add2f(vec2_t left, float x, float y)
{
	return ll_vec2_create2f(left.x + x, left.y + y);
}

/**
 * @return A new vector, that is the sum of @left
 * and { value, value }.
 */
LINEARLIBDEF vec2_t
ll_vec2_add1f(vec2_t left, float value)
{
	return ll_vec2_create2f(left.x + value, left.y + value);
}

/**
 * @return A new vector, that is the sum of @left and -@right.
 */
LINEARLIBDEF vec2_t
ll_vec2_sub2fv(vec2_t left, vec2_t right)
{
	return ll_vec2_add2f(left, -right.x, -right.y);
}

/**
 * @return A new vector, that is the sum of @left
 * and { -x, -y }.
 */
LINEARLIBDEF vec2_t
ll_vec2_sub2f(vec2_t left, float x, float y)
{
	return ll_vec2_add2f(left, -x, -y);
}

/**
 * @return A new vector, that is the sum of @left
 * and { -value, -value }.
 */
LINEARLIBDEF vec2_t
ll_vec2_sub1f(vec2_t left, float value)
{
	return ll_vec2_add2f(left, -value, -value);
}

/**
 * @return A new vector that is the component-wise
 * product of @left and @right.
 */
LINEARLIBDEF vec2_t
ll_vec2_mul2fv(vec2_t left, vec2_t right)
{
	return ll_vec2_create2f(left.x * right.x,
				left.y * right.y);
}

/**
 * @return A new vector that is the component-wise
 * product of @left and { x, y }.
 */
LINEARLIBDEF vec2_t
ll_vec2_mul2f(vec2_t left, float x, float y)
{
	return ll_vec2_create2f(left.x * x,
				left.y * y);
}

/**
 * @return A new vector, that is the component-wise
 * product of @left and { value, value }
 */
LINEARLIBDEF vec2_t
ll_vec2_mul1f(vec2_t left, float value)
{
	return ll_vec2_create2f(left.x * value,
				left.y * value);
}

/**
 * @return A new vector, that is the component-wise
 * division of @left and @right.
 */
LINEARLIBDEF vec2_t
ll_vec2_div2fv(vec2_t left, vec2_t right)
{
	return ll_vec2_create2f(left.x / right.x,
				left.y / right.y);
}

/**
 * @return A new vector, that is the component-wise
 * division of @left and { x, y }.
 */
LINEARLIBDEF vec2_t
ll_vec2_div2f(vec2_t left, float x, float y)
{
	return ll_vec2_create2f(left.x / x,
				left.y / y);
}

/**
 * @return A new vector, that is the component-wise
 * division of @left and { value, value }.
 */
LINEARLIBDEF vec2_t
ll_vec2_div1f(vec2_t left, float value)
{
	return ll_vec2_create2f(left.x / value,
				left.y / value);
}

/**
 * @return A new vector, that is the dot product
 * of @left and @right.
 */
LINEARLIBDEF float
ll_vec2_dot2fv(vec2_t left, vec2_t right)
{
	return left.x * right.x + left.y * right.y;
}

/**
 * @return A new vector, that is the dot product
 * of @left and { x, y }.
 */
LINEARLIBDEF float
ll_vec2_dot2f(vec2_t left, float x, float y)
{
	return left.x * x + left.y * y;
}

/**
 * @return A new vector, that is the cross
 * product of @left and @right.
 */
LINEARLIBDEF float
ll_vec2_cross2fv(vec2_t left, vec2_t right)
{
	return left.x * right.y - left.y * right.x;
}

/**
 * @return A new vector, that is the cross
 * product of @left and { x, y }.
 */
LINEARLIBDEF float
ll_vec2_cross2f(vec2_t left, float x, float y)
{
	return left.x * y - left.y * x;
}

/**
 * @return A new vector, that is normalised so
 * that the returned vector's length is 1. It is
 * the unit vector of @vec.
 */
LINEARLIBDEF vec2_t
ll_vec2_normalise2fv(vec2_t vec)
{
	float length = ll_vec2_length2fv(vec);
	return ll_vec2_create2f( vec.x / length, vec.y / length );
}

/**
 * @return A new vector, that is normalised so
 * that the returned vector's length is 1. It is
 * the unit vector of { x, y }
 */
LINEARLIBDEF vec2_t
ll_vec2_normalise2f(float x, float y)
{
	float length = ll_vec2_length2f(x, y);
	return ll_vec2_create2f( x / length, y / length );
}

/**
 * @return A vec3_t structure.
 */
LINEARLIBDEF vec3_t
ll_vec3_create3f(float x, float y, float z)
{
	return (vec3_t) { x, y, z };
}

/**
 * @return A copy of @vec.
 */
LINEARLIBDEF vec3_t
ll_vec3_create3fv(vec3_t vec)
{
	return vec;
}

/**
 * @return The length of @vec.
 */
LINEARLIBDEF float
ll_vec3_length3fv(vec3_t vec)
{
	return sqrtf(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
}

/**
 * @return The length of @vec.
 */
LINEARLIBDEF float
ll_vec3_length3f(float x, float y, float z)
{
	return sqrtf(x * x + y * y + z * z);
}

/**
 * @return The squared length of @vec.
 */
LINEARLIBDEF float
ll_vec3_length_squared3fv(vec3_t vec)
{
	float length = ll_vec3_length3fv(vec);
	return length * length;
}

/**
 * @return The squared length of the vector { x, y, z }.
 */
LINEARLIBDEF float
ll_vec3_length_squared3f(float x, float y, float z)
{
	float length = ll_vec3_length3f(x, y, z);
	return length * length;
}

/**
 * @return A new vector, that is the sum of @left
 * and @right.
 */
LINEARLIBDEF vec3_t
ll_vec3_add3fv(vec3_t left, vec3_t right)
{
	return ll_vec3_create3f(left.x + right.x, left.y + right.y,
				left.z + right.z);
}

/**
 * @return A new vector, that is the sum of @left
 * and { x, y, z }.
 */
LINEARLIBDEF vec3_t
ll_vec3_add3f(vec3_t left, float x, float y, float z)
{
	return ll_vec3_create3f(left.x + x, left.y + y,
				left.z + z);
}

/**
 * @return A new vector, that is the sum of @left
 * and { value, value, value }.
 */
LINEARLIBDEF vec3_t
ll_vec3_add1f(vec3_t left, float value)
{
	return ll_vec3_create3f(left.x + value, left.y + value,
				left.z + value);
}

/**
 * @return A new vector, that is the sum of @left
 * and -@right.
 */
LINEARLIBDEF vec3_t
ll_vec3_sub3fv(vec3_t left, vec3_t right)
{
	return ll_vec3_create3f(left.x - right.x, left.y - right.y,
				left.z - right.z);
}

/**
 * @return A new vector, that is the sum of @left
 * and { -x, -y, -z }.
 */
LINEARLIBDEF vec3_t
ll_vec3_sub3f(vec3_t left, float x, float y, float z)
{
	return ll_vec3_create3f(left.x - x, left.y - y,
				left.z - z);
}

/**
 * @return A new vector, that is the sum of @left
 * and { -value, -value, -value }.
 */
LINEARLIBDEF vec3_t
ll_vec3_sub1f(vec3_t left, float value)
{
	return ll_vec3_create3f(left.x - value, left.y - value,
				left.z - value);
}

/**
 * @return A new vector, that is the component-wise
 * product of @left and @right.
 */
LINEARLIBDEF vec3_t
ll_vec3_mul3fv(vec3_t left, vec3_t right)
{
	return ll_vec3_create3f(left.x * right.x, left.y * right.y,
				left.z * right.z);
}

/**
 * @return A new vector, that is the component-wise
 * product of @left and { x, y, z }.
 */
LINEARLIBDEF vec3_t
ll_vec3_mul3f(vec3_t left, float x, float y, float z)
{
	return ll_vec3_create3f(left.x * x, left.y * y,
				left.z * z);
}

/**
 * @return A new vector, that is the component-wise
 * product of @left and { value, value, value }.
 */
LINEARLIBDEF vec3_t
ll_vec3_mul1f(vec3_t left, float value)
{
	return ll_vec3_create3f(left.x * value, left.y * value,
				left.z * value);
}

/**
 * @return A new vector, that is the component-wise
 * division of @left and @right.
 */
LINEARLIBDEF vec3_t
ll_vec3_div3fv(vec3_t left, vec3_t right)
{
	return ll_vec3_create3f(left.x / right.x, left.y / right.y,
				left.z / right.z);
}

/**
 * @return A new vector, that is the component-wise
 * division of @left and { x, y, z }.
 */
LINEARLIBDEF vec3_t
ll_vec3_div3f(vec3_t left, float x, float y, float z)
{
	return ll_vec3_create3f(left.x / x, left.y / y,
				left.z / z);
}

/**
 * @return A new vector, that is the component-wise
 * division of @left and { value, value, value }.
 */
LINEARLIBDEF vec3_t
ll_vec3_div1f(vec3_t left, float value)
{
	return ll_vec3_create3f(left.x / value, left.y / value,
				left.z / value);
}

/**
 * @return The dot product of @left and @right.
 */
LINEARLIBDEF float
ll_vec3_dot3fv(vec3_t left, vec3_t right)
{
	return left.x * right.x + left.y * right.y + left.z * right.z;
}

/**
 * @return The dot product of @left and { x, y, z }.
 */
LINEARLIBDEF float
ll_vec3_dot3f(vec3_t left, float x, float y, float z)
{
	return left.x * x + left.y * y + left.z * z;
}

/**
 * @return The cross product of @left and @right,
 * this is essentially a new vector that is perpendicular
 * to @left and @right.
 */
LINEARLIBDEF vec3_t
ll_vec3_cross3fv(vec3_t left, vec3_t right)
{
	return ll_vec3_create3f(left.y*right.z - left.z*right.y,
				left.z*right.x - left.x*right.z,
				left.x*right.y - left.y*right.x);
}

/**
 * @return The cross product of @left and { x, y, z },
 * this is essentially a new vector that is perpendicular
 * to @left and { x, y, z }.
 */
LINEARLIBDEF vec3_t
ll_vec3_cross3f(vec3_t left, float x, float y, float z)
{
	return ll_vec3_create3f(left.y*z - left.z*y,
				left.z*x - left.x*z,
				left.x*y - left.y*x);
}

/**
 * @return A normalised vector, such that it's
 * length is 1. This is essentially the unit vector
 * of @vec.
 */
LINEARLIBDEF vec3_t
ll_vec3_normalise3fv(vec3_t vec)
{
	float length = ll_vec3_length3fv(vec);
	return ll_vec3_create3f( vec.x / length, vec.y / length,
				 vec.z / length );
}

/**
 * @return A normalised vector, such that it's
 * length is 1. This is essentially the unit vector
 * of { x, y, z }.
 */
LINEARLIBDEF vec3_t
ll_vec3_normalise3f(float x, float y, float z)
{
	float length = ll_vec3_length3f(x, y, z);
	return ll_vec3_create3f( x / length, y / length,
				 z / length );
}

/**
 * @return A vec3_t structure.
 */
LINEARLIBDEF vec4_t
ll_vec4_create4f(float x, float y, float z, float w)
{
	return (vec4_t) { x, y, z, w };
}

/**
 * @return A copy of @vec.
 */
LINEARLIBDEF vec4_t
ll_vec4_create4fv(vec4_t vec)
{
	return vec;
}

/**
 * @return The length of @vec.
 */
LINEARLIBDEF float
ll_vec4_length4fv(vec4_t vec)
{
	return sqrtf(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z
		     + vec.w * vec.w);
}

/**
 * @return The length of @vec.
 */
LINEARLIBDEF float
ll_vec4_length4f(float x, float y, float z, float w)
{
	return sqrtf(x * x + y * y + z * z + w * w);
}

/**
 * @return The squared length of @vec.
 */
LINEARLIBDEF float
ll_vec4_length_squared4fv(vec4_t vec)
{
	float length = ll_vec4_length4fv(vec);
	return length * length;
}

/**
 * @return The squared length of the vector { x, y, z, w }.
 */
LINEARLIBDEF float
ll_vec4_length_squared4f(float x, float y, float z, float w)
{
	float length = ll_vec4_length4f(x, y, z, w);
	return length * length;
}

/**
 * @return A new vector, that is the sum of @left
 * and @right.
 */
LINEARLIBDEF vec4_t
ll_vec4_add4fv(vec4_t left, vec4_t right)
{
	return ll_vec4_create4f(left.x + right.x, left.y + right.y,
				left.z + right.z, left.w + right.w);
}

/**
 * @return A new vector, that is the sum of @left
 * and { x, y, z, w }.
 */
LINEARLIBDEF vec4_t
ll_vec4_add4f(vec4_t left, float x, float y, float z, float w)
{
	return ll_vec4_create4f(left.x + x, left.y + y,
				left.z + z, left.w + w);
}

/**
 * @return A new vector, that is the sum of @left
 * and { value, value, value, value }.
 */
LINEARLIBDEF vec4_t
ll_vec4_add1f(vec4_t left, float value)
{
	return ll_vec4_create4f(left.x + value, left.y + value,
				left.z + value, left.w + value);
}

/**
 * @return A new vector, that is the sum of @left
 * and -@right.
 */
LINEARLIBDEF vec4_t
ll_vec4_sub4fv(vec4_t left, vec4_t right)
{
	return ll_vec4_create4f(left.x - right.x, left.y - right.y,
				left.z - right.z, left.w - right.w);
}

/**
 * @return A new vector, that is the sum of @left
 * and { -x, -y, -z, -w }.
 */
LINEARLIBDEF vec4_t
ll_vec4_sub4f(vec4_t left, float x, float y, float z, float w)
{
	return ll_vec4_create4f(left.x - x, left.y - y,
				left.z - z, left.w - w);
}

/**
 * @return A new vector, that is the sum of @left
 * and { -value, -value, -value, -value}.
 */
LINEARLIBDEF vec4_t
ll_vec4_sub1f(vec4_t left, float value)
{
	return ll_vec4_create4f(left.x - value, left.y - value,
				left.z - value, left.w - value);
}

/**
 * @return A new vector, that is the component-wise
 * product of @left and @right.
 */
LINEARLIBDEF vec4_t
ll_vec4_mul4fv(vec4_t left, vec4_t right)
{
	return ll_vec4_create4f(left.x * right.x, left.y * right.y,
				left.z * right.z, left.w * right.w);
}

/**
 * @return A new vector, that is the component-wise
 * product of @left and { x, y, z }.
 */
LINEARLIBDEF vec4_t
ll_vec4_mul4f(vec4_t left, float x, float y, float z, float w)
{
	return ll_vec4_create4f(left.x * x, left.y * y,
				left.z * z, left.w * w);
}

/**
 * @return A new vector, that is the component-wise
 * product of @left and { value, value, value, value }.
 */
LINEARLIBDEF vec4_t
ll_vec4_mul1f(vec4_t left, float value)
{
	return ll_vec4_create4f(left.x * value, left.y * value,
				left.z * value, left.w * value);
}

/**
 * @return A new vector, that is the component-wise
 * division of @left and @right.
 */
LINEARLIBDEF vec4_t
ll_vec4_div4fv(vec4_t left, vec4_t right)
{
	return ll_vec4_create4f(left.x / right.x, left.y / right.y,
				left.z / right.z, left.w / right.w);
}

/**
 * @return A new vector, that is the component-wise
 * division of @left and { x, y, z, w }.
 */
LINEARLIBDEF vec4_t
ll_vec4_div4f(vec4_t left, float x, float y, float z, float w)
{
	return ll_vec4_create4f(left.x / x, left.y / y,
				left.z / z, left.w / w);
}

/**
 * @return A new vector, that is the component-wise
 * division of @left and { value, value, value, value }.
 */
LINEARLIBDEF vec4_t
ll_vec4_div1f(vec4_t left, float value)
{
	return ll_vec4_create4f(left.x / value, left.y / value,
				left.z / value, left.w / value);
}

/**
 * @return The dot product of @left and @right.
 */
LINEARLIBDEF float
ll_vec4_dot4fv(vec4_t left, vec4_t right)
{
	return left.x * right.x + left.y * right.y
		+ left.z * right.z + left.w * right.w;
}

/**
 * @return The dot product of @left and { x, y, z, w }.
 */
LINEARLIBDEF float
ll_vec4_dot4f(vec4_t left, float x, float y, float z, float w)
{
	return left.x * x + left.y * y + left.z * z
		+ left.w * w;
}

/**
 * @return A normalised vector, such that it's
 * length is 1. This is essentially the unit vector
 * of @vec.
 */
LINEARLIBDEF vec4_t
ll_vec4_normalise4fv(vec4_t vec)
{
	float length = ll_vec4_length4fv(vec);
	return ll_vec4_create4f( vec.x / length, vec.y / length,
				 vec.z / length, vec.w / length);
}

/**
 * @return A normalised vector, such that it's
 * length is 1. This is essentially the unit vector
 * of { x, y, z, w }.
 */
LINEARLIBDEF vec4_t
ll_vec4_normalise4f(float x, float y, float z, float w)
{
	float length = ll_vec4_length4f(x, y, z, w);
	return ll_vec4_create4f( x / length, y / length,
				 z / length, w / length );
}

/**
 * @return A ivec2_t structure.
 */
LINEARLIBDEF ivec2_t
ll_ivec2_create2i(int x, int y)
{
	return (ivec2_t) { x, y };
}

/**
 * @description Essentially creates a new copy of @vec.
 *
 * @return A ivec2_t structure.
 */
LINEARLIBDEF ivec2_t
ll_ivec2_create2iv(ivec2_t vec)
{
	return vec;
}

/**
 * @return The length of @vec.
 */
LINEARLIBDEF float
ll_ivec2_length2iv(ivec2_t vec)
{
	return sqrtf((vec.x * vec.x) + (vec.y * vec.y));
}

/**
 * @return The length of a vector, represented by { x, y }.
 */
LINEARLIBDEF float
ll_ivec2_length2i(int x, int y)
{
	return sqrtf((x * x) + (y * y));
}

/**
 * @return The squared length of @vec.
 */
LINEARLIBDEF float
ll_ivec2_length_squared2iv(ivec2_t vec)
{
	float length = ll_ivec2_length2iv(vec);
	return length*length;
}

/**
 * @return The squared length of vector, represented by { x, y }.
 */
LINEARLIBDEF float
ll_ivec2_length_squared2i(int x, int y)
{
	float length = ll_ivec2_length2i(x, y);
	return length*length;
}

/**
 * @return A new vector, that is the sum of @left and @right.
 */
LINEARLIBDEF ivec2_t
ll_ivec2_add2iv(ivec2_t left, ivec2_t right)
{
	return ll_ivec2_create2i(left.x + right.x, left.y + right.y);
}

/**
 * @return A new vector, that is the sum of @left and { x, y }.
 */
LINEARLIBDEF ivec2_t
ll_ivec2_add2i(ivec2_t left, int x, int y)
{
	return ll_ivec2_create2i(left.x + x, left.y + y);
}

/**
 * @return A new vector, that is the sum of @left
 * and { value, value }.
 */
LINEARLIBDEF ivec2_t
ll_ivec2_add1i(ivec2_t left, int value)
{
	return ll_ivec2_create2i(left.x + value, left.y + value);
}

/**
 * @return A new vector, that is the sum of @left and -@right.
 */
LINEARLIBDEF ivec2_t
ll_ivec2_sub2iv(ivec2_t left, ivec2_t right)
{
	return ll_ivec2_add2i(left, -right.x, -right.y);
}

/**
 * @return A new vector, that is the sum of @left
 * and { -x, -y }.
 */
LINEARLIBDEF ivec2_t
ll_ivec2_sub2i(ivec2_t left, int x, int y)
{
	return ll_ivec2_add2i(left, -x, -y);
}

/**
 * @return A new vector, that is the sum of @left
 * and { -value, -value }.
 */
LINEARLIBDEF ivec2_t
ll_ivec2_sub1i(ivec2_t left, int value)
{
	return ll_ivec2_add2i(left, -value, -value);
}

/**
 * @return A new vector that is the component-wise
 * product of @left and @right.
 */
LINEARLIBDEF ivec2_t
ll_ivec2_mul2iv(ivec2_t left, ivec2_t right)
{
	return ll_ivec2_create2i(left.x * right.x,
				 left.y * right.y);
}

/**
 * @return A new vector that is the component-wise
 * product of @left and { x, y }.
 */
LINEARLIBDEF ivec2_t
ll_ivec2_mul2i(ivec2_t left, int x, int y)
{
	return ll_ivec2_create2i(left.x * x,
				 left.y * y);
}

/**
 * @return A new vector, that is the component-wise
 * product of @left and { value, value }
 */
LINEARLIBDEF ivec2_t
ll_ivec2_mul1i(ivec2_t left, int value)
{
	return ll_ivec2_create2i(left.x * value,
				 left.y * value);
}

/**
 * @return A new vector, that is the component-wise
 * division of @left and @right.
 */
LINEARLIBDEF ivec2_t
ll_ivec2_div2iv(ivec2_t left, ivec2_t right)
{
	return ll_ivec2_create2i(left.x / right.x,
				 left.y / right.y);
}

/**
 * @return A new vector, that is the component-wise
 * division of @left and { x, y }.
 */
LINEARLIBDEF ivec2_t
ll_ivec2_div2i(ivec2_t left, int x, int y)
{
	return ll_ivec2_create2i(left.x / x,
				 left.y / y);
}

/**
 * @return A new vector, that is the component-wise
 * division of @left and { value, value }.
 */
LINEARLIBDEF ivec2_t
ll_ivec2_div1i(ivec2_t left, int value)
{
	return ll_ivec2_create2i(left.x / value,
				 left.y / value);
}

/**
 * @return A new vector, that is the dot product
 * of @left and @right.
 */
LINEARLIBDEF float
ll_ivec2_dot2iv(ivec2_t left, ivec2_t right)
{
	return left.x * right.x + left.y * right.y;
}

/**
 * @return A new vector, that is the dot product
 * of @left and { x, y }.
 */
LINEARLIBDEF float
ll_ivec2_dot2i(ivec2_t left, int x, int y)
{
	return left.x * x + left.y * y;
}

/**
 * @return A new vector, that is the cross
 * product of @left and @right.
 */
LINEARLIBDEF float
ll_ivec2_cross2iv(ivec2_t left, ivec2_t right)
{
	return left.x * right.y - left.y * right.x;
}

/**
 * @return A new vector, that is the cross
 * product of @left and { x, y }.
 */
LINEARLIBDEF float
ll_ivec2_cross2i(ivec2_t left, int x, int y)
{
	return left.x * y - left.y * x;
}

/**
 * @return A new vector, that is normalised so
 * that the returned vector's length is 1. It is
 * the unit vector of @vec.
 */
LINEARLIBDEF ivec2_t
ll_ivec2_normalise2iv(ivec2_t vec)
{
	float length = ll_ivec2_length2iv(vec);
	return ll_ivec2_create2i( vec.x / length, vec.y / length );
}

/**
 * @return A new vector, that is normalised so
 * that the returned vector's length is 1. It is
 * the unit vector of { x, y }
 */
LINEARLIBDEF ivec2_t
ll_ivec2_normalise2i(int x, int y)
{
	float length = ll_ivec2_length2i(x, y);
	return ll_ivec2_create2i( x / length, y / length );
}

/**
 * @return A ivec3_t structure.
 */
LINEARLIBDEF ivec3_t
ll_ivec3_create3i(int x, int y, int z)
{
	return (ivec3_t) { x, y, z };
}

/**
 * @return A copy of @vec.
 */
LINEARLIBDEF ivec3_t
ll_ivec3_create3iv(ivec3_t vec)
{
	return vec;
}

/**
 * @return The length of @vec.
 */
LINEARLIBDEF float
ll_ivec3_length3iv(ivec3_t vec)
{
	return sqrtf(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
}

/**
 * @return The length of @vec.
 */
LINEARLIBDEF float
ll_ivec3_length3i(int x, int y, int z)
{
	return sqrtf(x * x + y * y + z * z);
}

/**
 * @return The squared length of @vec.
 */
LINEARLIBDEF float
ll_ivec3_length_squared3iv(ivec3_t vec)
{
	float length = ll_ivec3_length3iv(vec);
	return length * length;
}

/**
 * @return The squared length of the vector { x, y, z }.
 */
LINEARLIBDEF float
ll_ivec3_length_squared3i(int x, int y, int z)
{
	float length = ll_ivec3_length3i(x, y, z);
	return length * length;
}

/**
 * @return A new vector, that is the sum of @left
 * and @right.
 */
LINEARLIBDEF ivec3_t
ll_ivec3_add3iv(ivec3_t left, ivec3_t right)
{
	return ll_ivec3_create3i(left.x + right.x, left.y + right.y,
				left.z + right.z);
}

/**
 * @return A new vector, that is the sum of @left
 * and { x, y, z }.
 */
LINEARLIBDEF ivec3_t
ll_ivec3_add3i(ivec3_t left, int x, int y, int z)
{
	return ll_ivec3_create3i(left.x + x, left.y + y,
				left.z + z);
}

/**
 * @return A new vector, that is the sum of @left
 * and { value, value, value }.
 */
LINEARLIBDEF ivec3_t
ll_ivec3_add1i(ivec3_t left, int value)
{
	return ll_ivec3_create3i(left.x + value, left.y + value,
				left.z + value);
}

/**
 * @return A new vector, that is the sum of @left
 * and -@right.
 */
LINEARLIBDEF ivec3_t
ll_ivec3_sub3iv(ivec3_t left, ivec3_t right)
{
	return ll_ivec3_create3i(left.x - right.x, left.y - right.y,
				left.z - right.z);
}

/**
 * @return A new vector, that is the sum of @left
 * and { -x, -y, -z }.
 */
LINEARLIBDEF ivec3_t
ll_ivec3_sub3i(ivec3_t left, int x, int y, int z)
{
	return ll_ivec3_create3i(left.x - x, left.y - y,
				left.z - z);
}

/**
 * @return A new vector, that is the sum of @left
 * and { -value, -value, -value }.
 */
LINEARLIBDEF ivec3_t
ll_ivec3_sub1i(ivec3_t left, int value)
{
	return ll_ivec3_create3i(left.x - value, left.y - value,
				left.z - value);
}

/**
 * @return A new vector, that is the component-wise
 * product of @left and @right.
 */
LINEARLIBDEF ivec3_t
ll_ivec3_mul3iv(ivec3_t left, ivec3_t right)
{
	return ll_ivec3_create3i(left.x * right.x, left.y * right.y,
				left.z * right.z);
}

/**
 * @return A new vector, that is the component-wise
 * product of @left and { x, y, z }.
 */
LINEARLIBDEF ivec3_t
ll_ivec3_mul3i(ivec3_t left, int x, int y, int z)
{
	return ll_ivec3_create3i(left.x * x, left.y * y,
				left.z * z);
}

/**
 * @return A new vector, that is the component-wise
 * product of @left and { value, value, value }.
 */
LINEARLIBDEF ivec3_t
ll_ivec3_mul1i(ivec3_t left, int value)
{
	return ll_ivec3_create3i(left.x * value, left.y * value,
				left.z * value);
}

/**
 * @return A new vector, that is the component-wise
 * division of @left and @right.
 */
LINEARLIBDEF ivec3_t
ll_ivec3_div3iv(ivec3_t left, ivec3_t right)
{
	return ll_ivec3_create3i(left.x / right.x, left.y / right.y,
				left.z / right.z);
}

/**
 * @return A new vector, that is the component-wise
 * division of @left and { x, y, z }.
 */
LINEARLIBDEF ivec3_t
ll_ivec3_div3i(ivec3_t left, int x, int y, int z)
{
	return ll_ivec3_create3i(left.x / x, left.y / y,
				left.z / z);
}

/**
 * @return A new vector, that is the component-wise
 * division of @left and { value, value, value }.
 */
LINEARLIBDEF ivec3_t
ll_ivec3_div1i(ivec3_t left, int value)
{
	return ll_ivec3_create3i(left.x / value, left.y / value,
				left.z / value);
}

/**
 * @return The dot product of @left and @right.
 */
LINEARLIBDEF float
ll_ivec3_dot3iv(ivec3_t left, ivec3_t right)
{
	return left.x * right.x + left.y * right.y + left.z * right.z;
}

/**
 * @return The dot product of @left and { x, y, z }.
 */
LINEARLIBDEF float
ll_ivec3_dot3i(ivec3_t left, int x, int y, int z)
{
	return left.x * x + left.y * y + left.z * z;
}

/**
 * @return The cross product of @left and @right,
 * this is essentially a new vector that is perpendicular
 * to @left and @right.
 */
LINEARLIBDEF ivec3_t
ll_ivec3_cross3iv(ivec3_t left, ivec3_t right)
{
	return ll_ivec3_create3i(left.y*right.z - left.z*right.y,
				left.z*right.x - left.x*right.z,
				left.x*right.y - left.y*right.x);
}

/**
 * @return The cross product of @left and { x, y, z },
 * this is essentially a new vector that is perpendicular
 * to @left and { x, y, z }.
 */
LINEARLIBDEF ivec3_t
ll_ivec3_cross3i(ivec3_t left, int x, int y, int z)
{
	return ll_ivec3_create3i(left.y*z - left.z*y,
				left.z*x - left.x*z,
				left.x*y - left.y*x);
}

/**
 * @return A normalised vector, such that it's
 * length is 1. This is essentially the unit vector
 * of @vec.
 */
LINEARLIBDEF ivec3_t
ll_ivec3_normalise3iv(ivec3_t vec)
{
	float length = ll_ivec3_length3iv(vec);
	return ll_ivec3_create3i( vec.x / length, vec.y / length,
				 vec.z / length );
}

/**
 * @return A normalised vector, such that it's
 * length is 1. This is essentially the unit vector
 * of { x, y, z }.
 */
LINEARLIBDEF ivec3_t
ll_ivec3_normalise3i(int x, int y, int z)
{
	float length = ll_ivec3_length3i(x, y, z);
	return ll_ivec3_create3i( x / length, y / length,
				  z / length );
}

/**
 * @return A vec3_t structure.
 */
LINEARLIBDEF ivec4_t
ll_ivec4_create4i(int x, int y, int z, int w)
{
	return (ivec4_t) { x, y, z, w };
}

/**
 * @return A copy of @vec.
 */
LINEARLIBDEF ivec4_t
ll_ivec4_create4iv(ivec4_t vec)
{
	return vec;
}

/**
 * @return The length of @vec.
 */
LINEARLIBDEF float
ll_ivec4_length4iv(ivec4_t vec)
{
	return sqrtf(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z
		     + vec.w * vec.w);
}

/**
 * @return The length of @vec.
 */
LINEARLIBDEF float
ll_ivec4_length4i(int x, int y, int z, int w)
{
	return sqrtf(x * x + y * y + z * z + w * w);
}

/**
 * @return The squared length of @vec.
 */
LINEARLIBDEF float
ll_ivec4_length_squared4iv(ivec4_t vec)
{
	float length = ll_ivec4_length4iv(vec);
	return length * length;
}

/**
 * @return The squared length of the vector { x, y, z, w }.
 */
LINEARLIBDEF float
ll_ivec4_length_squared4i(int x, int y, int z, int w)
{
	float length = ll_ivec4_length4i(x, y, z, w);
	return length * length;
}

/**
 * @return A new vector, that is the sum of @left
 * and @right.
 */
LINEARLIBDEF ivec4_t
ll_ivec4_add4iv(ivec4_t left, ivec4_t right)
{
	return ll_ivec4_create4i(left.x + right.x, left.y + right.y,
				left.z + right.z, left.w + right.w);
}

/**
 * @return A new vector, that is the sum of @left
 * and { x, y, z, w }.
 */
LINEARLIBDEF ivec4_t
ll_ivec4_add4i(ivec4_t left, int x, int y, int z, int w)
{
	return ll_ivec4_create4i(left.x + x, left.y + y,
				left.z + z, left.w + w);
}

/**
 * @return A new vector, that is the sum of @left
 * and { value, value, value, value }.
 */
LINEARLIBDEF ivec4_t
ll_ivec4_add1i(ivec4_t left, int value)
{
	return ll_ivec4_create4i(left.x + value, left.y + value,
				left.z + value, left.w + value);
}

/**
 * @return A new vector, that is the sum of @left
 * and -@right.
 */
LINEARLIBDEF ivec4_t
ll_ivec4_sub4iv(ivec4_t left, ivec4_t right)
{
	return ll_ivec4_create4i(left.x - right.x, left.y - right.y,
				left.z - right.z, left.w - right.w);
}

/**
 * @return A new vector, that is the sum of @left
 * and { -x, -y, -z, -w }.
 */
LINEARLIBDEF ivec4_t
ll_ivec4_sub4i(ivec4_t left, int x, int y, int z, int w)
{
	return ll_ivec4_create4i(left.x - x, left.y - y,
				left.z - z, left.w - w);
}

/**
 * @return A new vector, that is the sum of @left
 * and { -value, -value, -value, -value}.
 */
LINEARLIBDEF ivec4_t
ll_ivec4_sub1i(ivec4_t left, int value)
{
	return ll_ivec4_create4i(left.x - value, left.y - value,
				left.z - value, left.w - value);
}

/**
 * @return A new vector, that is the component-wise
 * product of @left and @right.
 */
LINEARLIBDEF ivec4_t
ll_ivec4_mul4iv(ivec4_t left, ivec4_t right)
{
	return ll_ivec4_create4i(left.x * right.x, left.y * right.y,
				left.z * right.z, left.w * right.w);
}

/**
 * @return A new vector, that is the component-wise
 * product of @left and { x, y, z }.
 */
LINEARLIBDEF ivec4_t
ll_ivec4_mul4i(ivec4_t left, int x, int y, int z, int w)
{
	return ll_ivec4_create4i(left.x * x, left.y * y,
				left.z * z, left.w * w);
}

/**
 * @return A new vector, that is the component-wise
 * product of @left and { value, value, value, value }.
 */
LINEARLIBDEF ivec4_t
ll_ivec4_mul1i(ivec4_t left, int value)
{
	return ll_ivec4_create4i(left.x * value, left.y * value,
				left.z * value, left.w * value);
}

/**
 * @return A new vector, that is the component-wise
 * division of @left and @right.
 */
LINEARLIBDEF ivec4_t
ll_ivec4_div4iv(ivec4_t left, ivec4_t right)
{
	return ll_ivec4_create4i(left.x / right.x, left.y / right.y,
				left.z / right.z, left.w / right.w);
}

/**
 * @return A new vector, that is the component-wise
 * division of @left and { x, y, z, w }.
 */
LINEARLIBDEF ivec4_t
ll_ivec4_div4i(ivec4_t left, int x, int y, int z, int w)
{
	return ll_ivec4_create4i(left.x / x, left.y / y,
				left.z / z, left.w / w);
}

/**
 * @return A new vector, that is the component-wise
 * division of @left and { value, value, value, value }.
 */
LINEARLIBDEF ivec4_t
ll_ivec4_div1i(ivec4_t left, int value)
{
	return ll_ivec4_create4i(left.x / value, left.y / value,
				left.z / value, left.w / value);
}

/**
 * @return The dot product of @left and @right.
 */
LINEARLIBDEF float
ll_ivec4_dot4iv(ivec4_t left, ivec4_t right)
{
	return left.x * right.x + left.y * right.y
		+ left.z * right.z + left.w * right.w;
}

/**
 * @return The dot product of @left and { x, y, z, w }.
 */
LINEARLIBDEF float
ll_ivec4_dot4i(ivec4_t left, int x, int y, int z, int w)
{
	return left.x * x + left.y * y + left.z * z
		+ left.w * w;
}

/**
 * @return A normalised vector, such that it's
 * length is 1. This is essentially the unit vector
 * of @vec.
 */
LINEARLIBDEF ivec4_t
ll_ivec4_normalise4iv(ivec4_t vec)
{
	float length = ll_ivec4_length4iv(vec);
	return ll_ivec4_create4i( vec.x / length, vec.y / length,
				 vec.z / length, vec.w / length);
}

/**
 * @return A normalised vector, such that it's
 * length is 1. This is essentially the unit vector
 * of { x, y, z, w }.
 */
LINEARLIBDEF ivec4_t
ll_ivec4_normalise4i(int x, int y, int z, int w)
{
	float length = ll_ivec4_length4i(x, y, z, w);
	return ll_ivec4_create4i( x / length, y / length,
				  z / length, w / length );
}

/* performs matrix multiplcation on the matrices @l and @r
   matrix multiplication is not commutative, so order matters
   @left corresponds to the matrix on the left.
   @right corresponds to the matrix on the right. */
void
ll_mat4_multiply(mat4_t *left, mat4_t *right)
{
        mat4_t final;
        if (!left || !right) return;
        for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                        float sum = 0.0;
                        for (int k = 0; k < 4; k++)
                                sum += left->data[i * 4 + k] * right->data[k * 4 + j];
                        final.data[i * 4 + j] = sum;
                }
        }
        ll_mat4_copy(left, &final);
}

/* copy the matrix contents of @from into @to */
void
ll_mat4_copy(mat4_t *to,  mat4_t *from)
{
        if (!to || !from) return;
        for (int i = 0; i < 16; i++)
                to->data[i] = from->data[i];
}

/* stores the identity matrix into @m */
void
ll_mat4_identity(mat4_t *mat)
{
        if (!mat) return;
        for (int i = 0; i < 16; i++)
                mat->data[i] = 0;
        mat->m00 = 1.0;
        mat->m11 = 1.0;
        mat->m22 = 1.0;
        mat->m33 = 1.0;
}

/* stores a translation matrix inside of @mat with @dx, @dy and @dz */
void
ll_mat4_translate3f(mat4_t *mat, float dx, float dy, float dz)
{
        if (!mat) return;
        ll_mat4_identity(mat);
        mat->m30 = dx;
        mat->m31 = dy;
        mat->m32 = dz;
}

/* same as above, but instead allows to supply a 3-dimensional
   vector as the 2nd argument, containing (dx, dy, dz) */
void
ll_mat4_translate3fv(mat4_t *mat, vec3_t vec)
{
        ll_mat4_translate3f(mat, vec.x, vec.y, vec.z);
}

/* multiplies @m by a scaling matrix with components (w, h, d)
   storing the result back into @m */
void
ll_mat4_scale3f(mat4_t *mat, float w, float h, float d)
{
        if (!mat) return;
        ll_mat4_identity(mat);
        mat->m00 = w;
        mat->m11 = h;
        mat->m22 = d;
}

/* same as above, but instead allows to supply a 3-dimensional
   vector as the 2nd arugment, containing (w, h, d) */
void
ll_mat4_scale3fv(mat4_t *mat, vec3_t vec)
{
        ll_mat4_scale3f(mat, vec.x, vec.y, vec.z);
}

/* found whilst reading over rougier/freetype-gl implementation.
   that can be found here: https://github.com/rougier/freetype-gl */
void
ll_mat4_rotate3f(mat4_t *mat, float x, float y, float z, float theta)
{
        float c, s, norm;
        if (!mat) return;

        c =    (float) cos(M_PI * theta/180.0);
        s =    (float) sin(M_PI * theta/180.0);
        norm = (float) sqrt(x*x+y*y+z*z);
        
        x /= norm;
        y /= norm;
        z /= norm;

        ll_mat4_identity(mat);
        
        mat->m00 = x*x*(1-c)+c;
        mat->m10 = y*x*(1-c)-z*s;
        mat->m20 = z*x*(1-c)+y*s;
        mat->m01 = x*y*(1-c)+z*s;
        mat->m11 = y*y*(1-c)+c;
        mat->m21 = z*y*(1-c)-x*s;
        mat->m02 = x*z*(1-c)-y*s;
        mat->m12 = y*z*(1-c)+x*s;
        mat->m22 = z*z*(1-c)+c;
}

/* same as above, but instead allows to supply a 3-dimensional
   vector as the 2nd argument, containing (x, y, z) */
void
ll_mat4_rotate3fv(mat4_t *mat, vec3_t vec, float theta)
{
        ll_mat4_rotate3f(mat, vec.x, vec.y, vec.z, theta);
}

/* stores the orthographic matrix into @m, details of how this works can be found online or
   found at https://github.com/wwotz/linearlib/README.md */
void
ll_mat4_orthographic(mat4_t *mat, float top, float right,
                     float bottom, float left, float near, float far)
{
        if (!mat || left == right || near == far || bottom == top) return;
        ll_mat4_identity(mat);
        mat->m00 = 2.0f/(right-left);
        mat->m30 = -(right+left)/(right-left);
        mat->m11 = 2.0f/(top-bottom);
        mat->m31 = -(top+bottom)/(top-bottom);
        mat->m22 = -2.0f/(far-near);
        mat->m32 = -(far+near)/(far-near);
}

/* stores the perspective matrix into @m, details of how this works can be found online or
   found at https://github.com/wwotz/linearlib/README.md */
void
ll_mat4_perspective(mat4_t *mat, float fovy, float aspect,
                    float near, float far)
{
        float w, h;
        if (!mat || near == far) return;
        
        h = (float) tan(fovy / 360.0 * M_PI) * near;
        w = h * aspect;
        ll_mat4_frustum(mat, -w, w, -h, h, near, far);
}

void
ll_mat4_frustum(mat4_t *mat, float top, float right,
                float bottom, float left, float near, float far)
{
        if (!mat || left == right || bottom == top || near == far) return;
        ll_mat4_identity(mat);
        mat->m00 = (2.0f*near)/(right-left);
        mat->m20 = (right+left)/(right-left);
        mat->m11 = (2.0f*near)/(top-bottom);
        mat->m21 = (top+bottom)/(top-bottom);
        mat->m22 = -(far+near)/(far-near);
        mat->m32 = -(2.0f*far*near)/(far-near);
        mat->m23 = -1.0f;
        mat->m33 = 0.0;
}

#ifdef LL_USE_MATRIX

static mat4_t ll_matrices[LL_MATRIX_COUNT];
static int ll_matrices_idx;

void
ll_matrix_mode(matrix_type_t type)
{
        if (type >= 0 && type < LL_MATRIX_COUNT)
                ll_matrices_idx = type;
}

void
ll_matrix_multiply(mat4_t *right)
{
        ll_mat4_multiply(ll_matrices+ll_matrices_idx, right);
}

void
ll_matrix_identity(void)
{
        ll_mat4_identity(ll_matrices+ll_matrices_idx);
}

void
ll_matrix_translate3f(float dx, float dy, float dz)
{
        mat4_t m;
        ll_mat4_translate3f(&m, dx, dy, dz);
        ll_matrix_multiply(&m);
}

void
ll_matrix_translate3fv(vec3_t vec)
{
        ll_matrix_translate3f(vec.x, vec.y, vec.z);
}

void
ll_matrix_scale3f(float w, float h, float d)
{
        mat4_t m;
        ll_mat4_scale3f(&m, w, h, d);
        ll_matrix_multiply(&m);
}

void
ll_matrix_scale3fv(vec3_t vec)
{
        ll_matrix_scale3f(vec.x, vec.y, vec.z);
}

void
ll_matrix_rotate3f(float x, float y, float z, float angle)
{
        mat4_t m;
        ll_mat4_rotate3f(&m, x, y, z, angle);
        ll_matrix_multiply(&m);
}

void
ll_matrix_rotate3fv(vec3_t vec, float angle)
{
        ll_matrix_rotate3f(vec.x, vec.y, vec.z, angle);
}

void
ll_matrix_orthographic(float top, float right,
                       float bottom, float left, float near, float far)
{
        ll_mat4_orthographic(ll_matrices+ll_matrices_idx, top,
                             right, bottom, left, near, far);
}

void
ll_matrix_perspective(float fovy, float aspect,
                      float near, float far)
{
        ll_mat4_perspective(ll_matrices+ll_matrices_idx,
                            fovy, aspect, near, far);
}

void
ll_matrix_frustum(float left, float right,
                  float bottom, float top, float near, float far)
{
        ll_mat4_frustum(ll_matrices+ll_matrices_idx, left, right,
                        bottom, top, near, far);
}

mat4_t
ll_matrix_get_copy(void)
{
        return ll_matrices[ll_matrices_idx];
}

#endif /* LL_USE_MATRIX */
#endif /* LINEARLIB_IMPLEMENTATION */
#endif /* LINEARLIB_H_ */
