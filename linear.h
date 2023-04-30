/*
  Author: wwotz
  Add #LINEARLIB_IMPLEMENTATION to the start of the implementation
  file in order to add the implementation code to your project. 
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

typedef union vec2_t {
        float data[2];
        struct {
                float x, y;
        };
        struct {
                float r, g;
        };
} vec2_t;

typedef union vec3_t {
        float data[3];
        struct {
                float x, y, z;
        };
        struct {
                float r, g, b;
        };        
} vec3_t;

typedef union vec4_t {
        float data[4];
        struct {
                float x, y, z, w;
        };
        struct {
                float r, g, b, a;
        };
} vec4_t;

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

LINEARLIBDEF vec2_t ll_vec2_create(float x, float y);
LINEARLIBDEF float  ll_vec2_length(vec2_t vec);
LINEARLIBDEF vec2_t ll_vec2_add(vec2_t left, vec2_t right);
LINEARLIBDEF vec2_t ll_vec2_sub(vec2_t left, vec2_t right);
LINEARLIBDEF vec2_t ll_vec2_mul(vec2_t left, vec2_t right);
LINEARLIBDEF vec2_t ll_vec2_div(vec2_t left, vec2_t right);
LINEARLIBDEF float  ll_vec2_dot(vec2_t left, vec2_t right);
LINEARLIBDEF float  ll_vec2_cross(vec2_t left, vec2_t right);
LINEARLIBDEF vec2_t ll_vec2_normalise(vec2_t vec);

LINEARLIBDEF vec3_t ll_vec3_create(float x, float y, float z);
LINEARLIBDEF float  ll_vec3_length(vec3_t vec);
LINEARLIBDEF vec3_t ll_vec3_add(vec3_t left, vec3_t right);
LINEARLIBDEF vec3_t ll_vec3_sub(vec3_t left, vec3_t right);
LINEARLIBDEF vec3_t ll_vec3_mul(vec3_t left, vec3_t right);
LINEARLIBDEF vec3_t ll_vec3_div(vec3_t left, vec3_t right);
LINEARLIBDEF float  ll_vec3_dot(vec3_t left, vec3_t right);
LINEARLIBDEF vec3_t ll_vec3_cross(vec3_t left, vec3_t right);
LINEARLIBDEF vec3_t ll_vec3_normalise(vec3_t vec);

LINEARLIBDEF vec4_t ll_vec4_create(float x, float y, float z, float w);
LINEARLIBDEF float  ll_vec4_length(vec4_t vec);
LINEARLIBDEF vec4_t ll_vec4_add(vec4_t left, vec4_t right);
LINEARLIBDEF vec4_t ll_vec4_sub(vec4_t left, vec4_t right);
LINEARLIBDEF vec4_t ll_vec4_mul(vec4_t left, vec4_t right);
LINEARLIBDEF vec4_t ll_vec4_div(vec4_t left, vec4_t right);
LINEARLIBDEF float  ll_vec4_dot(vec4_t left, vec4_t right);
LINEARLIBDEF vec4_t ll_vec4_normalise(vec4_t vec);

LINEARLIBDEF void ll_mat4_multiply(mat4_t *left, mat4_t *right);
LINEARLIBDEF void ll_mat4_copy(mat4_t *to, mat4_t *from);
LINEARLIBDEF void ll_mat4_identity(mat4_t *mat);
LINEARLIBDEF void ll_mat4_translate3f(mat4_t *mat, float dx, float dy, float dz);
LINEARLIBDEF void ll_mat4_translate3fv(mat4_t *mat, vec3_t vec);
LINEARLIBDEF void ll_mat4_scale3f(mat4_t *mat, float w, float h, float d);
LINEARLIBDEF void ll_mat4_scale3fv(mat4_t *mat, vec3_t vec);
LINEARLIBDEF void ll_mat4_rotate3f(mat4_t *mat, float x, float y, float z, float angle);
LINEARLIBDEF void ll_mat4_rotate3fv(mat4_t *mat, vec3_t vec, float angle);
LINEARLIBDEF void ll_mat4_orthographic(mat4_t *mat, float top, float right,
                                 float bottom, float left, float near, float far);
LINEARLIBDEF void  ll_mat4_perspective(mat4_t *mat, float fovy, float aspect,
                                 float near, float far);
LINEARLIBDEF void ll_mat4_frustum(mat4_t *mat, float left, float right,
                            float bottom, float top, float near, float far);

#ifdef LINEARLIB_IMPLEMENTATION

/* returns a 2-dimensional vector with (x, y) as it's components */
vec2_t ll_vec2_create(float x, float y)
{
        return (vec2_t) {x, y};
}

/* returns the length of @vec  */
float ll_vec2_length(vec2_t vec)
{
        return sqrtf((vec.x * vec.x) + (vec.y * vec.y));
}

/* returns a new vector containing the sum of first and second,
   component-wise sum */
vec2_t ll_vec2_add(vec2_t first, vec2_t second)
{
        return ll_vec2_create(first.x + second.x, first.y + second.y);
}

/* returns a new vector containing the subtraction of second from first,
   component-wise difference  */
vec2_t ll_vec2_sub(vec2_t first, vec2_t second)
{
        return ll_vec2_create(first.x - second.x, first.y - second.y);
}

/* returns a new vector containing the product of first and second,
   component-wise multiplication */
vec2_t ll_vec2_mul(vec2_t first, vec2_t second)
{
        return ll_vec2_create(first.x * second.x, first.y * second.y);
}

/* returns a new vector containing the division of first by second,
   component-wise division */
vec2_t ll_vec2_div(vec2_t first, vec2_t second)
{
        return ll_vec2_create(first.x / second.x, first.y / second.y);
}

/* returns the dot-product of @first and @second */
float ll_vec2_dot(vec2_t first, vec2_t second)
{
        return first.x * second.x + first.y * second.y;
}

/* returns the cross-product of @first and @second */
float ll_vec2_cross(vec2_t first, vec2_t second)
{
        return first.x * second.y - first.y * second.x;
}

/* returns a normalised vector of @vec */
vec2_t ll_vec2_normalise(vec2_t vec)
{
        float length = ll_vec2_length(vec);
        return (vec2_t) {vec.x / length, vec.y / length};
}

/* returns a 3-dimensional vector with (x, y, z) as it's components */
vec3_t ll_vec3_create(float x, float y, float z)
{
        return (vec3_t) {x, y, z};
}

/* returns the length of @vec  */
float ll_vec3_length(vec3_t vec)
{
        return sqrtf(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
}

/* returns a new vector containing the sum of first and second,
   component-wise sum */
vec3_t ll_vec3_add(vec3_t first, vec3_t second)
{
        return ll_vec3_create(first.x + second.x, first.y + second.y,
                              first.z + second.z);
}

/* returns a new vector containing the subtraction of second from first,
   component-wise difference  */
vec3_t ll_vec3_sub(vec3_t first, vec3_t second)
{
        return ll_vec3_create(first.x - second.x, first.y - second.y,
                              first.z - second.z);
}

/* returns a new vector containing the product of first and second,
   component-wise multiplication */
vec3_t ll_vec3_mul(vec3_t first, vec3_t second)
{
        return ll_vec3_create(first.x * second.x, first.y * second.y,
                              first.z * second.z);
}

/* returns a new vector containing the division of first by second,
   component-wise division */
vec3_t ll_vec3_div(vec3_t first, vec3_t second)
{
        return ll_vec3_create(first.x / second.x, first.y / second.y,
                              first.z / second.z);
}

/* returns the dot-product of @first and @second */
float ll_vec3_dot(vec3_t first, vec3_t second)
{
        return first.x * second.x + first.y * second.y + first.z * second.z;
}

/* returns the cross-product of @first and @second */
vec3_t ll_vec3_cross(vec3_t first, vec3_t second)
{
        return ll_vec3_create(first.y*second.z - first.z*second.y,
                              first.z*second.x - first.x*second.z,
                              first.x*second.y - first.y*second.x);
}

/* returns a normalised vector of @vec */
vec3_t ll_vec3_normalise(vec3_t vec)
{
        float length = ll_vec3_length(vec);
        return (vec3_t) {vec.x / length, vec.y / length, vec.z / length};
}

/* returns a 4-dimensional vector with (x, y, z, w) as it's components */
vec4_t ll_vec4_create(float x, float y, float z, float w)
{
        return (vec4_t) {x, y, z, w};
}

/* returns the length of @vec  */
float ll_vec4_length(vec4_t vec)
{
        return sqrtf(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z + vec.w * vec.w);
}

/* returns a new vector containing the sum of first and second,
   component-wise sum */
vec4_t ll_vec4_add(vec4_t first, vec4_t second)
{
        return ll_vec4_create(first.x + second.x, first.y + second.y,
                             first.z + second.z, first.w + second.w);
}

/* returns a new vector containing the subtraction of second from first,
   component-wise difference  */
vec4_t ll_vec4_sub(vec4_t first, vec4_t second)
{
        return ll_vec4_create(first.x - second.x, first.y - second.y,
                             first.z - second.z, first.w - second.w);
}

/* returns a new vector containing the product of first and second,
   component-wise multiplication */
vec4_t ll_vec4_mul(vec4_t first, vec4_t second)
{
        return ll_vec4_create(first.x * second.x, first.y * second.y,
                             first.z * second.z, first.w * second.w);
}

/* returns a new vector containing the division of first by second,
   component-wise division */
vec4_t ll_vec4_div(vec4_t first, vec4_t second)
{
        return ll_vec4_create(first.x / second.x, first.y / second.y,
                             first.z / second.z, first.w / second.w);
}

/* return the dot-product of @first and @second */
float ll_vec4_dot(vec4_t first, vec4_t second)
{
        return first.x * second.x + first.y * second.y
                + first.z * second.z + first.w * second.w;
}

/* returns a normalised vector of @vec */
vec4_t ll_vec4_normalise(vec4_t vec)
{
        float length = ll_vec4_length(vec);
        return (vec4_t) {vec.x / length, vec.y / length, vec.z / length, vec.w / length};
}

/* performs matrix multiplcation on the matrices @l and @r
   matrix multiplication is not commutative, so order matters
   @left corresponds to the matrix on the left.
   @right corresponds to the matrix on the right. */
void ll_mat4_multiply(mat4_t *left, mat4_t *right)
{
        mat4_t final;
        if (!left || !right) return;
        for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                        float sum = 0.0;
                        for (int k = 0; k < 4; k++)
                                sum += left->data[j * 4 + k] * right->data[i * 4 + k];
                        final.data[j * 4 + i] = sum;
                }
        }
        ll_mat4_copy(left, &final);
}

/* copy the matrix contents of @from into @to */
void ll_mat4_copy(mat4_t *to,  mat4_t *from)
{
        if (!to || !from) return;
        for (int i = 0; i < 16; i++)
                to->data[i] = from->data[i];
}

/* stores the identity matrix into @m */
void ll_mat4_identity(mat4_t *mat)
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
void ll_mat4_translate3f(mat4_t *mat, float dx, float dy, float dz)
{
        if (!mat) return;
        ll_mat4_identity(mat);
        mat->m30 = dx;
        mat->m31 = dy;
        mat->m32 = dz;
}

/* same as above, but instead allows to supply a 3-dimensional
   vector as the 2nd argument, containing (dx, dy, dz) */
void ll_mat4_translate3fv(mat4_t *mat, vec3_t vec)
{
        ll_mat4_translate3f(mat, vec.x, vec.y, vec.z);
}

/* multiplies @m by a scaling matrix with components (w, h, d)
   storing the result back into @m */
void ll_mat4_scale3f(mat4_t *mat, float w, float h, float d)
{
        if (!mat) return;
        ll_mat4_identity(mat);
        mat->m00 = w;
        mat->m11 = h;
        mat->m22 = d;
}

/* same as above, but instead allows to supply a 3-dimensional
   vector as the 2nd arugment, containing (w, h, d) */
void ll_mat4_scale3fv(mat4_t *mat, vec3_t vec)
{
        ll_mat4_scale3f(mat, vec.x, vec.y, vec.z);
}

/* found whilst reading over rougier/freetype-gl implementation.
   that can be found here: https://github.com/rougier/freetype-gl */
void ll_mat4_rotate3f(mat4_t *mat, float x, float y, float z, float theta)
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
void ll_mat4_rotate3fv(mat4_t *mat, vec3_t vec, float theta)
{
        ll_mat4_rotate3f(mat, vec.x, vec.y, vec.z, theta);
}

/* stores the orthographic matrix into @m, details of how this works can be found online or
   found at https://github.com/wwotz/linearlib/README.md */
void ll_mat4_orthographic(mat4_t *mat, float top, float right,
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
void ll_mat4_perspective(mat4_t *mat, float fovy, float aspect, float near, float far)
{
        float w, h;
        if (!mat || near == far) return;
        
        h = (float) tan(fovy / 360.0 * M_PI) * near;
        w = h * aspect;
        ll_mat4_frustum(mat, -w, w, -h, h, near, far);
}

void ll_mat4_frustum(mat4_t *mat, float top, float right,
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

#endif /* LINEARLIB_IMPLEMENTATION */
#endif // LINEARLIB_H_
