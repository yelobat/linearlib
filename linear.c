#include "linear.h"

inline vec2_t create_vec2_t(float x, float y)
{
        return (vec2_t) {x, y};
}

inline float length_vec2_t(vec2_t vec)
{
        return (float) sqrt((vec.x * vec.x) + (vec.y * vec.y));
}

vec2_t normalise_vec2_t(vec2_t vec)
{
        float length = length_vec2_t(vec);
        return (vec2_t) {vec.x / length, vec.y / length};
}

inline vec3_t create_vec3_t(float x, float y, float z)
{
        return (vec3_t) {x, y, z};
}

inline float length_vec3_t(vec3_t vec)
{
        return (float) sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
}

inline vec3_t normalise_vec3_t(vec3_t vec)
{
        float length = length_vec3_t(vec);
        return (vec3_t) {vec.x / length, vec.y / length, vec.z / length};
}

inline vec4_t create_vec4_t(float x, float y, float z, float w)
{
        return (vec4_t) {x, y, z, w};
}

inline float length_vec4_t(vec4_t vec)
{
        return (float) sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z + vec.w * vec.w);
}

inline vec4_t normalise_vec4_t(vec4_t vec)
{
        float length = length_vec4_t(vec);
        return (vec4_t) {vec.x / length, vec.y / length, vec.z / length, vec.w / length};
}

int multiply_matrix(float l[16], float r[16])
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

        for (int i = 0; i < 16; i++)
                l[i] = final[i];
        return 0;
}

int identity_matrix(float m[16])
{
        for (int i = 0; i < 16; i++)
                m[i] = 0;
        m[0]  = 1;
        m[5]  = 1;
        m[10] = 1;
        m[15] = 1;
        return 0;
}

int translate_matrix3f(float m[16], float dx, float dy, float dz)
{
        float other[16] = {
                1.0, 0.0, 0.0, dx,
                0.0, 1.0, 0.0, dy,
                0.0, 0.0, 1.0, dz,
                0.0, 0.0, 0.0, 1.0,
        };

        return multiply_matrix(m, other);
}

int translate_matrix3fv(float m[16], vec3_t vec)
{
        return translate_matrix3f(m, vec.x, vec.y, vec.z);
}

int scale_matrix3f(float m[16], float w, float h, float d)
{
        float other[16] = {
                w, 0, 0, 0,
                0, h, 0, 0,
                0, 0, d, 0,
                0, 0, 0, 1,
        };

        return multiply_matrix(m, other);
}

int scale_matrix3fv(float m[16], vec3_t vec)
{
        return scale_matrix3f(m, vec.x, vec.y, vec.z);
}

int rotate_matrix3f(float m[16], float x, float y, float z, float theta)
{
        float other[16] = {
                cos(theta) + x*x*(1-cos(theta)),   x*y*(1-cos(theta)) - z*sin(theta), x*z*(1-cos(theta)) + y*sin(theta), 0.0,
                y*x*(1-cos(theta)) + z*sin(theta), cos(theta) + y*y*(1-cos(theta)),   y*z*(1-cos(theta)) - x*sin(theta), 0.0,
                z*x*(1-cos(theta)) - y*sin(theta), y*z*(1-cos(theta)) + x*sin(theta), cos(theta) + x*x*(1-cos(theta)),   0.0,
                0.0,                               0.0,                               0.0,                               1.0,
        };

        return multiply_matrix(m, other);
}

int rotate_matrix3fv(float m[16], vec3_t vec, float theta)
{
        return rotate_matrix3f(m, vec.x, vec.y, vec.z, theta);
}

int orthographic_matrix(float m[16], float t, float r, float b, float l, float n, float f)
{
        float other[16] = {
                2.0/(r-l), 0.0,       0.0,        -(r+l)/(r-l),
                0.0,       2.0/(t-b), 0.0,        -(t+b)/(t-b),
                0.0,       0.0,       -2.0/(f-n), -(f+n)/(f-n),
                0.0,       0.0,       0.0,        1.0
        };

        identity_matrix(m);
        return multiply_matrix(m, other);
}

int perspective_matrix(float m[16], float t, float r, float b, float l, float n, float f)
{
        float other[16] = {
                2.0*n/(r-l), 0.0,         (r+l)/(r-l),        0.0,
                0.0,         2.0*n/(t-b), (t+b)/(t-b),        0.0,
                0.0,         0.0,         -(f+n)/(f-n),       -2*(f*n)/(f-n),
                0.0,         0.0,         -1.0,               0.0
        };

        identity_matrix(m);
        return multiply_matrix(m, other);
}
