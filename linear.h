#ifndef LINEAR_H_
#define LINEAR_H_

//requires -lm flag to link the math library
#include <math.h>

typedef struct vec2_t {
        float x, y;
} vec2_t;

typedef struct vec3_t {
        float x, y, z;
} vec3_t;

typedef struct vec4_t {
        float x, y, z, w;
} vec4_t;

extern vec2_t vec2_create(float x, float y);
extern float  vec2_length(vec2_t vec);
extern vec2_t vec2_add(vec2_t first, vec2_t second);
extern vec2_t vec2_sub(vec2_t first, vec2_t second);
extern vec2_t vec2_mul(vec2_t first, vec2_t second);
extern vec2_t vec2_div(vec2_t first, vec2_t second);
extern float  vec2_dot(vec2_t first, vec2_t second);
extern float  vec2_cross(vec2_t first, vec2_t second);
extern vec2_t vec2_normalise(vec2_t vec);

extern vec3_t vec3_create(float x, float y, float z);
extern float  vec3_length(vec3_t vec);
extern vec3_t vec3_add(vec3_t first, vec3_t second);
extern vec3_t vec3_sub(vec3_t first, vec3_t second);
extern vec3_t vec3_mul(vec3_t first, vec3_t second);
extern vec3_t vec3_div(vec3_t first, vec3_t second);
extern float  vec3_dot(vec3_t first, vec3_t second);
extern vec3_t vec3_cross(vec3_t first, vec3_t second);
extern vec3_t vec3_normalise(vec3_t vec);

extern vec4_t vec4_create(float x, float y, float z, float w);
extern float  vec4_length(vec4_t vec);
extern vec4_t vec4_add(vec4_t first, vec4_t second);
extern vec4_t vec4_sub(vec4_t first, vec4_t second);
extern vec4_t vec4_mul(vec4_t first, vec4_t second);
extern vec4_t vec4_div(vec4_t first, vec4_t second);
extern float  vec4_dot(vec4_t first, vec4_t second);
extern vec4_t vec4_normalise(vec4_t vec);

extern int  matrix_multiply(float l[16], float r[16]);
extern void matrix_copy(float to[16], float from[16]);
extern int  matrix_identity(float m[16]);
extern int  matrix_translate3f(float m[16], float dx, float dy, float dz);
extern int  matrix_translate3fv(float m[16], vec3_t vec);
extern int  matrix_scale3f(float m[16], float w, float h, float d);
extern int  matrix_scale3fv(float m[16], vec3_t vec);
extern int  matrix_rotate3f(float m[16], float x, float y, float z, float angle);
extern int  matrix_rotate3fv(float m[16], vec3_t, float angle);
extern int  matrix_orthographic(float m[16], float t, float r, float b, float l, float n, float f);
extern int  matrix_perspective(float m[16], float t, float r, float b, float l, float n, float f);
extern int matrix_lookat(float m[16], vec3_t at, vec3_t eye, vec3_t up);

#ifdef LINEARLIB_IMPLEMENTATION

/* returns a 2-dimensional vector with (x, y) as it's components */
inline vec2_t vec2_create(float x, float y)
{
        return (vec2_t) {x, y};
}

/* returns the length of @vec  */
inline float vec2_length(vec2_t vec)
{
        return sqrtf((vec.x * vec.x) + (vec.y * vec.y));
}

/* returns a new vector containing the sum of first and second,
   component-wise sum */
inline vec2_t vec2_add(vec2_t first, vec2_t second)
{
        return vec2_create(first.x + second.x, first.y + second.y);
}

/* returns a new vector containing the subtraction of second from first,
   component-wise difference  */
inline vec2_t vec2_sub(vec2_t first, vec2_t second)
{
        return vec2_create(first.x - second.x, first.y - second.y);
}

/* returns a new vector containing the product of first and second,
   component-wise multiplication */
inline vec2_t vec2_mul(vec2_t first, vec2_t second)
{
        return vec2_create(first.x * second.x, first.y * second.y);
}

/* returns a new vector containing the division of first by second,
   component-wise division */
inline vec2_t vec2_div(vec2_t first, vec2_t second)
{
        return vec2_create(first.x / second.x, first.y / second.y);
}

/* returns the dot-product of @first and @second */
inline float vec2_dot(vec2_t first, vec2_t second)
{
        return first.x * second.x + first.y * second.y;
}

/* returns the cross-product of @first and @second */
inline float vec2_cross(vec2_t first, vec2_t second)
{
        return first.x * second.y - first.y * second.x;
}

/* returns a normalised vector of @vec */
inline vec2_t vec2_normalise(vec2_t vec)
{
        float length = vec2_length(vec);
        return (vec2_t) {vec.x / length, vec.y / length};
}

/* returns a 3-dimensional vector with (x, y, z) as it's components */
inline vec3_t vec3_create(float x, float y, float z)
{
        return (vec3_t) {x, y, z};
}

/* returns the length of @vec  */
inline float vec3_length(vec3_t vec)
{
        return sqrtf(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
}

/* returns a new vector containing the sum of first and second,
   component-wise sum */
inline vec3_t vec3_add(vec3_t first, vec3_t second)
{
        return vec3_create(first.x + second.x, first.y + second.y,
                             first.z + second.z);
}

/* returns a new vector containing the subtraction of second from first,
   component-wise difference  */
inline vec3_t vec3_sub(vec3_t first, vec3_t second)
{
        return vec3_create(first.x - second.x, first.y - second.y,
                             first.z - second.z);
}

/* returns a new vector containing the product of first and second,
   component-wise multiplication */
inline vec3_t vec3_mul(vec3_t first, vec3_t second)
{
        return vec3_create(first.x * second.x, first.y * second.y,
                             first.z * second.z);
}

/* returns a new vector containing the division of first by second,
   component-wise division */
inline vec3_t vec3_div(vec3_t first, vec3_t second)
{
        return vec3_create(first.x / second.x, first.y / second.y,
                             first.z / second.z);
}

/* returns the dot-product of @first and @second */
inline float vec3_dot(vec3_t first, vec3_t second)
{
        return first.x * second.x + first.y * second.y + first.z * second.z;
}

/* returns the cross-product of @first and @second */
inline vec3_t vec3_cross(vec3_t first, vec3_t second)
{
        return vec3_create(first.y*second.z - first.z*second.y,
                             first.z*second.x - first.x*second.z,
                             first.x*second.y - first.y*second.x);
}

/* returns a normalised vector of @vec */
inline vec3_t vec3_normalise(vec3_t vec)
{
        float length = vec3_length(vec);
        return (vec3_t) {vec.x / length, vec.y / length, vec.z / length};
}

/* returns a 4-dimensional vector with (x, y, z, w) as it's components */
inline vec4_t vec4_create(float x, float y, float z, float w)
{
        return (vec4_t) {x, y, z, w};
}

/* returns the length of @vec  */
inline float vec4_length(vec4_t vec)
{
        return sqrtf(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z + vec.w * vec.w);
}

/* returns a new vector containing the sum of first and second,
   component-wise sum */
inline vec4_t vec4_add(vec4_t first, vec4_t second)
{
        return vec4_create(first.x + second.x, first.y + second.y,
                             first.z + second.z, first.w + second.w);
}

/* returns a new vector containing the subtraction of second from first,
   component-wise difference  */
inline vec4_t vec4_sub(vec4_t first, vec4_t second)
{
        return vec4_create(first.x - second.x, first.y - second.y,
                             first.z - second.z, first.w - second.w);
}

/* returns a new vector containing the product of first and second,
   component-wise multiplication */
inline vec4_t vec4_mul(vec4_t first, vec4_t second)
{
        return vec4_create(first.x * second.x, first.y * second.y,
                             first.z * second.z, first.w * second.w);
}

/* returns a new vector containing the division of first by second,
   component-wise division */
inline vec4_t vec4_div(vec4_t first, vec4_t second)
{
        return vec4_create(first.x / second.x, first.y / second.y,
                             first.z / second.z, first.w / second.w);
}

/* return the dot-product of @first and @second */
inline float vec4_dot(vec4_t first, vec4_t second)
{
        return first.x * second.x + first.y * second.y
                + first.z * second.z + first.w * second.w;
}

/* returns a normalised vector of @vec */
inline vec4_t vec4_normalise(vec4_t vec)
{
        float length = vec4_length(vec);
        return (vec4_t) {vec.x / length, vec.y / length, vec.z / length, vec.w / length};
}

/* performs matrix multiplcation on the matrices @l and @r
   matrix multiplication is not commutative, so order matters
   @l corresponds to the matrix on the left.
   @r corresponds to the matrix on the right. */
int matrix_multiply(float l[16], float r[16])
{
        float final[16];
        for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                        float sum = 0.0;
                        for (int k = 0; k < 4; k++)
                                sum += l[j * 4 + k] * r[i * 4 + k];
                        final[j * 4 + i] = sum;
                }
        }

        matrix_copy(l, final);
        return 0;
}

/* copy the matrix contents of @from into @to */
void matrix_copy(float to[16], float from[16])
{
        for (int i = 0; i < 16; i++)
                to[i] = from[i];
}

/* stores the identity matrix into @m */
int matrix_identity(float m[16])
{
        for (int i = 0; i < 16; i++)
                m[i] = 0;
        m[0]  = 1;
        m[5]  = 1;
        m[10] = 1;
        m[15] = 1;
        return 0;
}

/* multiplies @m by a translation matrix with components (dx, dy, dz)
   storing the result back into @m */
int matrix_translate3f(float m[16], float dx, float dy, float dz)
{
        float other[16] = {
                1.0, 0.0, 0.0, dx,
                0.0, 1.0, 0.0, dy,
                0.0, 0.0, 1.0, dz,
                0.0, 0.0, 0.0, 1.0,
        };

        return matrix_multiply(m, other);
}

/* same as above, but instead allows to supply a 3-dimensional
   vector as the 2nd argument, containing (dx, dy, dz) */
int matrix_translate3fv(float m[16], vec3_t vec)
{
        return matrix_translate3f(m, vec.x, vec.y, vec.z);
}

/* multiplies @m by a scaling matrix with components (w, h, d)
   storing the result back into @m */
int matrix_scale3f(float m[16], float w, float h, float d)
{
        float other[16] = {
                w, 0, 0, 0,
                0, h, 0, 0,
                0, 0, d, 0,
                0, 0, 0, 1,
        };

        return matrix_multiply(m, other);
}

/* same as above, but instead allows to supply a 3-dimensional
   vector as the 2nd arugment, containing (w, h, d) */
int matrix_scale3fv(float m[16], vec3_t vec)
{
        return matrix_scale3f(m, vec.x, vec.y, vec.z);
}

/* needs to some testing */
int matrix_rotate3f(float m[16], float x, float y, float z, float theta)
{
        float other[16] = {
                cos(theta) + x*x*(1-cos(theta)),   x*y*(1-cos(theta)) - z*sin(theta), x*z*(1-cos(theta)) + y*sin(theta), 0.0,
                y*x*(1-cos(theta)) + z*sin(theta), cos(theta) + y*y*(1-cos(theta)),   y*z*(1-cos(theta)) - x*sin(theta), 0.0,
                z*x*(1-cos(theta)) - y*sin(theta), y*z*(1-cos(theta)) + x*sin(theta), cos(theta) + x*x*(1-cos(theta)),   0.0,
                0.0,                               0.0,                               0.0,                               1.0,
        };

        return matrix_multiply(m, other);
}

int matrix_rotatex3f(float m[16], float theta)
{
        float other[16] = {
                1.0, 0.0, 0.0, 0.0,
                0.0, cos(theta), -sin(theta), 0.0,
                0.0, sin(theta), cos(theta), 0.0,
                0.0, 0.0, 0.0, 1.0
        };

        return matrix_multiply(m, other);
}

int matrix_rotatey3f(float m[16], float theta)
{
        float other[16] = {
                cos(theta), 0.0, sin(theta), 0.0,
                0.0, 1.0, 0.0, 0.0,
                -sin(theta), 0.0, cos(theta), 0.0,
                0.0, 0.0, 0.0, 1.0
        };

        return matrix_multiply(m, other);
}

/* same as above, but instead allows to supply a 3-dimensional
   vector as the 2nd argument, containing (x, y, z) */
int matrix_rotate3fv(float m[16], vec3_t vec, float theta)
{
        return matrix_rotate3f(m, vec.x, vec.y, vec.z, theta);
}

/* stores the orthographic matrix into @m, details of how this works can be found online or
   found at https://github.com/wwotz/linearlib/README.md */
int matrix_orthographic(float m[16], float t, float r, float b, float l, float n, float f)
{
        float other[16] = {
                2.0/(r-l), 0.0,       0.0,        -(r+l)/(r-l),
                0.0,       2.0/(t-b), 0.0,        -(t+b)/(t-b),
                0.0,       0.0,       -2.0/(f-n), -(f+n)/(f-n),
                0.0,       0.0,       0.0,        1.0
        };

        matrix_identity(m);
        return matrix_multiply(m, other);
}

/* stores the perspective matrix into @m, details of how this works can be found online or
   found at https://github.com/wwotz/linearlib/README.md */
int matrix_perspective(float m[16], float t, float r, float b, float l, float n, float f)
{
        float other[16] = {
                2.0*n/(r-l), 0.0, 2.0*(l+r)/(r-l), 0.0,
                0.0, 2.0*n/(b-t), 2.0*(t+b)/(b-t), 0.0,
                0.0, 0.0, -(f+n)/(f-n), 2.0*(f*n)/(f-n),
                0.0, 0.0, -1.0, 0.0
        };

        matrix_identity(m);
        return matrix_multiply(m, other);
}

int matrix_lookat(float m[16], vec3_t at, vec3_t eye, vec3_t up)
{
        vec3_t cam_direction = vec3_normalise(vec3_sub(at, eye));
        vec3_t cam_right = vec3_normalise(vec3_cross(up, cam_direction));
        vec3_t cam_up = vec3_cross(cam_direction, cam_right);

        float other[16] = {
                cam_right.x, cam_right.y, cam_right.z, 0.0,
                cam_up.x, cam_up.y, cam_up.z, 0.0,
                cam_direction.x, cam_direction.y, cam_direction.z, 0.0,
                0.0, 0.0, 0.0, 1.0,
        };

        return matrix_multiply(m, other);
}


#endif /* LINEARLIB_IMPLEMENTATION */

#endif // LINEAR_H_
