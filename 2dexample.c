/* very simple maze generation demo using linearlib, SDL2 and OpenGL */

#define LINEARLIB_IMPLEMENTATION
#include "linear.h"
#include <SDL2/SDL.h>
#include <GL/glew.h>
#include <SDL2/SDL_opengl.h>
#include <sys/time.h>

#define WINDOW_NAME "2d Example"
#define WINDOW_X SDL_WINDOWPOS_UNDEFINED
#define WINDOW_Y SDL_WINDOWPOS_UNDEFINED
#define WINDOW_W 1200
#define WINDOW_H 800
#define WINDOW_FLAGS (SDL_WINDOW_SHOWN | SDL_WINDOW_OPENGL)

static SDL_Window *window = 0;
static SDL_GLContext *context = 0;
static int running = 0;
static SDL_Event event = {0};

float model[16], projection[16], view[16];

static const char *vert_src =
        "#version 330 core\n"
        "uniform mat4 model;\n"
        "uniform mat4 projection;\n"
        "uniform mat4 view;\n"
        "layout (location = 0) in vec2 vertex;\n"
        "void main(){\n"
        "gl_Position = projection * view * model * vec4(vertex, 0.0, 1.0);\n"
        "}";

static const char *frag_src =
        "#version 330 core\n"
        "uniform vec4 colour;\n"
        "void main(){\n"
        "gl_FragColor = colour;"
        "}";

static GLuint shader_create(const char *src, GLenum type)
{
        GLuint shader;
        shader = glCreateShader(type);
        glShaderSource(shader, 1, (const GLchar **)&src, NULL);
        glCompileShader(shader);
        GLuint status;
        glGetShaderiv(shader, GL_COMPILE_STATUS, &status);
        if (status != GL_TRUE) {
                fprintf(stderr, "Failed to create shader!\n");
                exit(EXIT_FAILURE);
        }
        return shader;
}

static GLuint program_create(const char *vert_src, const char *frag_src)
{
        GLuint vert_shader, frag_shader, program;
        vert_shader = shader_create(vert_src, GL_VERTEX_SHADER);
        frag_shader = shader_create(frag_src, GL_FRAGMENT_SHADER);
        program = glCreateProgram();
        glAttachShader(program, vert_shader);
        glAttachShader(program, frag_shader);
        glLinkProgram(program);
        glDeleteShader(vert_shader);
        glDeleteShader(frag_shader);

        GLuint status;
        glGetProgramiv(program, GL_LINK_STATUS, &status);
        if (status != GL_TRUE) {
                fprintf(stderr, "Failed to create program!\n");
                exit(EXIT_FAILURE);
        }
        return program;
}

typedef enum {
        COLOUR_R = 0,
        COLOUR_G,
        COLOUR_B,
        COLOUR_A,
        COLOUR_MASK_COUNT
} COLOUR_MASK;

static float colour_mask(GLuint colour, COLOUR_MASK mask, GLboolean normalize)
{
        float col = 0.0;
        switch (mask) {
        case COLOUR_R:
                col = (colour & 0xff000000) >> 24;
                return (normalize ? col/255.0 : col);
        case COLOUR_G:
                col = (colour & 0x00ff0000) >> 16;
                return (normalize ? col/255.0 : col);
        case COLOUR_B:
                col = (colour & 0x0000ff00) >> 8;
                return (normalize ? col/255.0 : col);
        default:
        case COLOUR_A:
                col = (colour & 0x000000ff) >> 0;
                return (normalize ? col/255.0 : col);
        }
}

typedef struct {
        GLuint vbo, vao, ebo, program, colour;
        GLfloat x, y, w, h;
} rect_t;

static void rect_create(rect_t *rect, GLfloat x, GLfloat y, GLfloat w, GLfloat h,
                        GLuint program, GLuint colour)
{
        static vec2_t vertices[4] = {
                { 0.0, 0.0 },
                { 1.0, 0.0 },
                { 1.0, 1.0 },
                { 0.0, 1.0 }
        };

        static GLuint indices[6] = {
                0, 1, 2,
                2, 3, 0
        };

        rect->x = x;
        rect->y = y;
        rect->w = w;
        rect->h = h;
        rect->colour = colour;

        rect->program = program;
        glUseProgram(program);

        glGenVertexArrays(1, &rect->vao);
        glBindVertexArray(rect->vao);

        glGenBuffers(1, &rect->vbo);
        glBindBuffer(GL_ARRAY_BUFFER, rect->vbo);
        glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

        glGenBuffers(1, &rect->ebo);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, rect->ebo);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray(0);

        glBindVertexArray(0);
        glUseProgram(0);
}

static void rect_render(rect_t *rect)
{
        glUseProgram(rect->program);
        matrix_identity(model);
        matrix_scale3f(model, rect->w, rect->h, 1.0);
        matrix_translate3f(model, rect->x, rect->y, 0.0);
        glUniformMatrix4fv(glGetUniformLocation(rect->program, "model"), 1, GL_FALSE, model);
        glUniformMatrix4fv(glGetUniformLocation(rect->program, "view"), 1, GL_FALSE, view);
        glUniformMatrix4fv(glGetUniformLocation(rect->program, "projection"), 1, GL_FALSE, projection);
        glUniform4f(glGetUniformLocation(rect->program, "colour"),
                    colour_mask(rect->colour, COLOUR_R, GL_TRUE),
                    colour_mask(rect->colour, COLOUR_G, GL_TRUE),
                    colour_mask(rect->colour, COLOUR_B, GL_TRUE),
                    colour_mask(rect->colour, COLOUR_A, GL_TRUE));
        glBindVertexArray(rect->vao);
        glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
        glBindVertexArray(0);
        glUseProgram(0);
}

static void rect_free(rect_t *rect)
{
        glDeleteBuffers(1, &rect->vbo);
        glDeleteBuffers(1, &rect->ebo);
        glDeleteVertexArrays(1, &rect->vao);
        glDeleteProgram(rect->program);
        rect->x = 0;
        rect->y = 0;
        rect->w = 0;
        rect->h = 0;
}

static struct timeval start, curr;

void timer_reset(void)
{
        gettimeofday(&start, NULL);
        gettimeofday(&curr, NULL);
}

void timer_start(void)
{
        timer_reset();
}

double timer_get_time(void)
{
        gettimeofday(&curr, NULL);
        double time = (curr.tv_sec - start.tv_sec)
                + (double) (curr.tv_usec - start.tv_usec)/1000000;
        return time;
}

int main(int argc, char **argv)
{
        SDL_Init(SDL_INIT_EVERYTHING);
        window = SDL_CreateWindow(WINDOW_NAME,
                                  WINDOW_X,
                                  WINDOW_Y,
                                  WINDOW_W,
                                  WINDOW_H,
                                  WINDOW_FLAGS);
        context = SDL_GL_CreateContext(window);
        glewInit();

        glClearColor(0.1, 0.1, 0.1, 1.0);
        glViewport(0, 0, WINDOW_W, WINDOW_H);

        matrix_identity(model);
        matrix_identity(view);
        matrix_orthographic(projection, 0.0, WINDOW_W, WINDOW_H, 0.0, 10.0, -10.0);

        timer_start();

        GLuint program = program_create(vert_src, frag_src);
        rect_t rect[120 * 80];
        for (int i = 0; i < 120; i++) {
                for (int j = 0; j < 80; j++) {
                        rect_create(rect+i+(j*120), i*10.0, j*10.0, 10.0, 10.0, program, 0x000000ff);
                }
        }

        srand(60);

        //select random start point
        int x = (int) 119.0* (float) rand()/RAND_MAX;
        int y = (int) 79.0 * (float) rand()/RAND_MAX;

        GLboolean visited[120 * 80] = {0};
        GLuint stack_capacity = 120 * 80;
        GLuint stack[120 * 80] = {0};
        GLuint direction[120 * 80] = {0};
        GLuint stack_size = 0;

        typedef enum {
                LEFT = 0,
                RIGHT,
                UP,
                DOWN,
                DIRECTION_COUNT,
        } DIRECTION;

        stack[stack_size++] = (x + y*120);

        running = 1;
        while (running) {
                while (SDL_PollEvent(&event)) {
                        if (event.type == SDL_QUIT) {
                                running = 0;
                                break;
                        }
                }

                //perform new search each frame
                if (stack_size != 0) {
                        int index = stack[--stack_size];
                        int dir = direction[stack_size];
                        visited[index] = 1;
                        int y = index / 120;
                        int x = index % 120;
                        rect[index].colour = 0xff0000ff;
                        if (dir == LEFT) {
                                rect[(x+1) + y*120].colour = 0xff0000ff;
                                visited[(x+1) + y*120] = 1;
                        } else if (dir == RIGHT) {
                                rect[(x-1) + y*120].colour = 0xff0000ff;
                                visited[(x-1) + y*120] = 1;
                        } else if (dir == UP) {
                                rect[x + (y-1)*120].colour = 0xff0000ff;
                                visited[x + (y-1)*120] = 1;
                        } else {
                                rect[x + (y+1)*120].colour = 0xff0000ff;
                                visited[x + (y+1)*120] = 1;
                        }
                        float random = (float) rand() / RAND_MAX;
                        if (random < 0.25 && x-2 >= 0 && !visited[(x-2) + y*120]) {
                                stack[stack_size++] = (x-2) + y*120;
                                direction[stack_size-1] = LEFT;

                                if (x+2 <= 119 && !visited[(x+2) + y*120]) {
                                        stack[stack_size++] = (x+2) + y*120;
                                        direction[stack_size-1] = RIGHT;
                                }

                                if (y-2 >= 0 && !visited[x + (y-2)*120]) {
                                        stack[stack_size++] = (x) + (y-2)*120;
                                        direction[stack_size-1] = DOWN;
                                }

                                if (y+2 <= 79 && !visited[(x) + (y+2)*120]) {
                                        stack[stack_size++] = (x) + (y+2)*120;
                                        direction[stack_size-1] = UP;
                                }
                        } else if (random < 0.5 && x+2 <= 119 && !visited[(x+2) + y*120]) {
                                stack[stack_size++] = (x+2) + y*120;
                                direction[stack_size-1] = RIGHT;
                                if (x-2 >= 0 && !visited[(x-2) + y*120]) {
                                        stack[stack_size++] = (x-2) + y*120;
                                        direction[stack_size-1] = LEFT;
                                }

                                if (y-2 >= 0 && !visited[x + (y-2)*120]) {
                                        stack[stack_size++] = (x) + (y-2)*120;
                                        direction[stack_size-1] = DOWN;
                                }

                                if (y+2 <= 79 && !visited[(x) + (y+2)*120]) {
                                        stack[stack_size++] = (x) + (y+2)*120;
                                        direction[stack_size-1] = UP;
                                }
                        } else if (random < 0.75 && y-2 >= 0 && !visited[(x) + (y-2)*120]) {
                                stack[stack_size++] = (x) + (y-2)*120;
                                direction[stack_size-1] = DOWN;
                                if (x-2 >= 0 && !visited[(x-2) + y*120]) {
                                        stack[stack_size++] = (x-2) + y*120;
                                        direction[stack_size-1] = LEFT;
                                }

                                if (x+2 <= 119 && !visited[(x+2) + y*120]) {
                                        stack[stack_size++] = (x+2) + y*120;
                                        direction[stack_size-1] = RIGHT;
                                }

                                if (y+2 <= 79 && !visited[(x) + (y+2)*120]) {
                                        stack[stack_size++] = (x) + (y+2)*120;
                                        direction[stack_size-1] = UP;
                                }
                        } else if (y+2 <= 79 && !visited[(x) + (y+2)*120]) {
                                stack[stack_size++] = (x) + (y+2)*120;
                                direction[stack_size-1] = UP;

                                if (x-2 >= 0 && !visited[(x-2) + y*120]) {
                                        stack[stack_size++] = (x-2) + y*120;
                                        direction[stack_size-1] = LEFT;
                                }

                                if (x+2 <= 119 && !visited[(x+2) + y*120]) {
                                        stack[stack_size++] = (x+2) + y*120;
                                        direction[stack_size-1] = RIGHT;
                                }

                                if (y-2 >= 0 && !visited[x + (y-2)*120]) {
                                        stack[stack_size++] = (x) + (y-2)*120;
                                        direction[stack_size-1] = DOWN;
                                }
                        }
                }

                glClear(GL_COLOR_BUFFER_BIT);
                for (int i = 0; i < 120 * 80; i++)
                        rect_render(rect+i);
                SDL_GL_SwapWindow(window);
        }

        SDL_DestroyWindow(window);
        SDL_GL_DeleteContext(context);
        SDL_Quit();
        return 0;
}
