#include "TXLib.h"
#include <stdio.h>
#include <stdint.h>
#include <stddef.h>

#define TIME_MEASUREMENT

static const float window_width  = 800.f;
static const float window_height = 600.f;

static const float dx = 1/window_width;
static const float dy = 1/window_width;

typedef RGBQUAD (&scr_t)[(unsigned long long)window_height][(unsigned long long)window_width];

static const float x_centre = -1.325f;
static const float y_centre = 0;

static const int N_iterations_max = 256;
static const float r2_max = 100.f;

const size_t running_num = 100;

volatile int N_iterations = 0;

extern "C" uint64_t get_time();

void view_regulation (float* x_shift, float* y_shift);
void draw_pixels(int N_iterations, RGBQUAD scr[(size_t)window_height][(size_t)window_width], size_t ix, size_t iy);
void count_mandelbrot(float x_shift, float y_shift, RGBQUAD scr[(size_t)window_height][(size_t)window_width]);



int main() {

    float x_shift = 0.f;
    float y_shift = 0.f;

    Win32::_fpreset();

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

            count_mandelbrot(x_shift, y_shift, scr);

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

void draw_pixels(int N_iterations, RGBQUAD scr[(size_t)window_height][(size_t)window_width], size_t ix, size_t iy) {

    float I = sqrtf (sqrtf ((float)N_iterations / (float)N_iterations_max)) * 255.f;
    //float I = (N % 2) * 255.f;

    BYTE c = (BYTE)I;

    RGBQUAD color = (N_iterations < N_iterations_max) ? RGBQUAD {(BYTE) (255 - c), (BYTE) (c%2 * 64), c} : RGBQUAD {};
    //RGBQUAD color = (N < nMax) ? RGBQUAD {(BYTE) c*6, 0, c*10} : RGBQUAD {};
    //RGBQUAD color = (N < nMax) ? RGBQUAD {c, c, c} : RGBQUAD {};
    scr[iy][ix] = color;

}

void count_mandelbrot(float x_shift, float y_shift, RGBQUAD scr[(size_t)window_height][(size_t)window_width]) {

    for (size_t iy = 0; iy < (size_t)window_height; iy++) {

        if (txGetAsyncKeyState (VK_ESCAPE))
            break;

        float x0 = ((          -  (window_width / 2)) * dx + x_centre + x_shift);
        float y0 = (((float)iy - (window_height / 2)) * dy + y_centre + y_shift);

        for (size_t ix = 0; ix < (size_t)window_width; ix++, x0 += dx) {

            float x = x0,
                  y = y0;

            N_iterations = 0;

            for (;N_iterations < N_iterations_max; N_iterations++) {

                float x2 = x*x,
                      y2 = y*y,
                      xy = x*y;

                float r2 = x2 + y2;

                if (r2 >= r2_max)
                    break;

                x = x2 - y2 + x0;
                y = xy + xy + y0;
            }

            #ifndef TIME_MEASUREMENT
            draw_pixels(N_iterations, scr, ix, iy);
            #endif

        }
    }
}



