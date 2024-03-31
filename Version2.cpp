#include "TXLib.h"
#include <stdio.h>

//#define TIME_MEASUREMENT

static const float window_width  = 800.f;
static const float window_height = 600.f;

static const float dx = 1/window_width;
static const float dy = 1/window_width;

typedef RGBQUAD (&scr_t)[(unsigned long long)window_height][(unsigned long long)window_width];

static const float x_centre = -1.325f;
static const float y_centre = 0;

static const int N_iterations_max = 256;
static const float R2_max = 100.f;

const size_t running_num = 100;

const size_t points_num = 8;

volatile float N_iter = 0;

extern "C" uint64_t get_time();

void view_regulation (float* x_shift, float* y_shift);
void draw_pixels(float* N_iterations, RGBQUAD scr[(size_t)window_height][(size_t)window_width], size_t ix, size_t iy);
void count_mandelbrot(float x_shift, float y_shift, RGBQUAD scr[(size_t)window_height][(size_t)window_width]);

int main() {

    Win32::_fpreset();

    float x_shift = 0.f;
    float y_shift = 0.f;

    #ifndef TIME_MEASUREMENT

    txCreateWindow(window_width, window_height);
    txBegin();
    #endif

    scr_t scr = (scr_t) *txVideoMemory();

    #ifndef TIME_MEASUREMENT

    for (;;) {

        if (txGetAsyncKeyState (VK_ESCAPE))
            break;

        view_regulation(&x_shift, &y_shift);

        #else

        uint64_t general_time = 0;
        for (size_t run_num = 0; run_num < running_num; run_num++) {

            uint64_t start_time = get_time();
            #endif

            count_mandelbrot( x_shift, y_shift, scr);

            #ifdef TIME_MEASUREMENT
            uint64_t time_taken = get_time() - start_time;
            general_time += time_taken;
        }

        printf("%llu\n", general_time/running_num);
        #else

        printf("\t\t\rFPS: %.0lf", txGetFPS());
        txSleep();
    }
    #endif

}

void view_regulation (float* x_shift, float* y_shift) {

    if (txGetAsyncKeyState (VK_RIGHT))
        *x_shift += dx * (txGetAsyncKeyState (VK_SHIFT) ? 100.f : 10.f);

    if (txGetAsyncKeyState (VK_LEFT))
        *x_shift -= dx * (txGetAsyncKeyState (VK_SHIFT) ? 100.f : 10.f);

    if (txGetAsyncKeyState (VK_UP))
        *y_shift += dy * (txGetAsyncKeyState (VK_SHIFT) ? 100.f : 10.f);

    if (txGetAsyncKeyState (VK_DOWN))
        *y_shift -= dy * (txGetAsyncKeyState (VK_SHIFT) ? 100.f : 10.f);


}


void draw_pixels(float* N_iterations, RGBQUAD scr[(size_t)window_height][(size_t)window_width], size_t ix, size_t iy) {

    for (size_t i = 0; i < points_num; i++) {

        float I = sqrtf (sqrtf ((float)N_iterations[i] / (float)N_iterations_max)) * 255.f;

        BYTE c = (BYTE)I;

        RGBQUAD color = (N_iterations[i] < N_iterations_max) ? RGBQUAD {(BYTE) (255 - c), (BYTE) (c%2 * 64), c} : RGBQUAD {};
        scr[iy][ix+i] = color;
    }

}

void count_mandelbrot(float x_shift, float y_shift, RGBQUAD scr[(size_t)window_height][(size_t)window_width]) {

    for (size_t iy = 0; iy < (size_t)window_height; iy++) {

        if (txGetAsyncKeyState (VK_ESCAPE))
            break;

        float x0 = ((          -  (window_width / 2)) * dx + x_centre + x_shift);
        float y0 = (((float)iy - (window_height / 2)) * dy + y_centre + y_shift);

        for (size_t ix = 0; ix < (size_t)window_width; ix += points_num, x0 += points_num*dx) {

            float X0[points_num] = {x0, x0 + dx, x0 + 2*dx, x0 + 3*dx, x0 + 4*dx, x0 + 5*dx, x0 + 6*dx, x0 + 7*dx},
                  Y0[points_num] = {y0, y0, y0, y0, y0, y0, y0, y0};

            float X[points_num] = {}; for (size_t i = 0; i < points_num; i++) X[i] = X0[i];
            float Y[points_num] = {}; for (size_t i = 0; i < points_num; i++) Y[i] = Y0[i];

            float N_iterations[points_num] = {0, 0, 0, 0, 0, 0, 0, 0};


            for (size_t n_iter = 0; n_iter < N_iterations_max; n_iter++) {

                float X2[points_num] = {}; for (size_t i = 0; i < points_num; i++) X2[i] = X[i] * X[i];
                float Y2[points_num] = {}; for (size_t i = 0; i < points_num; i++) Y2[i] = Y[i] * Y[i];
                float XY[points_num] = {}; for (size_t i = 0; i < points_num; i++) XY[i] = X[i] * Y[i];

                float R2[points_num] = {}; for (size_t i = 0; i < points_num; i++) R2[i] = X2[i] + Y2[i];

                int cmp[points_num] = {};
                for (size_t i = 0; i < points_num; i++) {
                    if (R2[i] <= R2_max)
                        cmp[i] = 1;
                }

                int mask = 0;
                for (size_t i = 0; i < points_num; i++) {
                    mask |= cmp[i] << i;
                    if (!mask)
                        break;
                }

                for (size_t i = 0; i < points_num; i++) N_iterations[i] += (float)cmp[i];

                for (size_t i = 0; i < points_num; i++) X[i] = X2[i] - Y2[i] + X0[i];
                for (size_t i = 0; i < points_num; i++) Y[i] = XY[i] + XY[i] + Y0[i];
            }
            N_iter = N_iterations[0];

            #ifndef TIME_MEASUREMENT
            draw_pixels(N_iterations,  scr, ix, iy);
            #endif

        }
    }
}
