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

vec2_t create_vec2_t(float x, float y);
float length_vec2_t(vec2_t vec);
float dot_vec2_t(vec2_t first, vec2_t second);
float cross_vec2_t(vec2_t first, vec2_t second);
vec2_t normalise_vec2_t(vec2_t vec);

vec3_t create_vec3_t(float x, float y, float z);
float length_vec3_t(vec3_t vec);
float dot_vec3_t(vec3_t first, vec3_t second);
vec3_t cross_vec3_t(vec3_t first, vec3_t second);
vec3_t normalise_vec3_t(vec3_t vec);

vec4_t create_vec4_t(float x, float y, float z, float w);
float length_vec4_t(vec4_t vec);
float dot_vec4_t(vec4_t first, vec4_t second);
vec4_t normalise_vec4_t(vec4_t vec);

int multiply_matrix(float l[16], float r[16]);
int identity_matrix(float m[16]);
int translate_matrix3f(float m[16], float dx, float dy, float dz);
int translate_matrix3fv(float m[16], vec3_t vec);
int scale_matrix3f(float m[16], float w, float h, float d);
int scale_matrix3fv(float m[16], vec3_t vec);
int rotate_matrix3f(float m[16], float x, float y, float z, float angle);
int rotate_matrix3fv(float m[16], vec3_t, float angle);
int orthographic_matrix(float m[16], float t, float r, float b, float l, float n, float f);
int perspective_matrix(float m[16], float t, float r, float b, float l, float n, float f);

#endif // LINEAR_H_
