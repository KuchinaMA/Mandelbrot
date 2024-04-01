#include "TXLib.h"
#include <stdio.h>

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

const size_t points_num = 8;

volatile int N_iter = 0;

extern "C" uint64_t get_time();

void view_regulation (float* x_shift, float* y_shift);
void draw_pixels(int* N_iterations, RGBQUAD scr[(size_t)window_height][(size_t)window_width], size_t ix, size_t iy);
void count_mandelbrot(float x_shift, float y_shift, RGBQUAD scr[(size_t)window_height][(size_t)window_width], float* multipliers, float* R2_max);

inline void mm256_set_ps (float res[points_num], float val7, float val6, float val5, float val4, float val3, float val2, float val1, float val0) __attribute__((always_inline));
inline void mm256_set1_ps (float res[points_num], float val) __attribute__((always_inline));
inline void mm256_setzero_si256 (int res[points_num]) __attribute__((always_inline));

inline void mm256_cpy_ps (float res[points_num], const float src[points_num]) __attribute__((always_inline));

inline void mm256_add_ps (float res[points_num], const float reg1[points_num], const float reg2[points_num]) __attribute__((always_inline));
inline void mm256_sub_ps (float res[points_num], const float reg1[points_num], const float reg2[points_num]) __attribute__((always_inline));
inline void mm256_mul_ps (float res[points_num], const float reg1[points_num], const float reg2[points_num]) __attribute__((always_inline));

inline void mm256_sub_epi32 (int res[points_num], const int reg1[points_num], const float reg2[points_num]) __attribute__((always_inline));

inline void mm256_cmple_ps (float cmp[points_num], const float reg1[points_num], const float reg2[points_num]) __attribute__((always_inline));

inline int mm256_movemask_ps (const float cmp[points_num]) __attribute__((always_inline));

int main() {

    Win32::_fpreset();

    float x_shift = 0.f;
    float y_shift = 0.f;

    float R2_max[points_num] = {};
    mm256_set1_ps(R2_max, r2_max);

    float multipliers[points_num] = {};
    mm256_set_ps(multipliers, 7.f, 6.f, 5.f, 4.f, 3.f, 2.f, 1.f, 0.f);

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

            count_mandelbrot(x_shift, y_shift, scr, multipliers, R2_max);

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

void draw_pixels(int* N_iterations, RGBQUAD scr[(size_t)window_height][(size_t)window_width], size_t ix, size_t iy) {

    for (size_t i = 0; i < points_num; i++) {

        float I = sqrtf (sqrtf ((float)N_iterations[i] / (float)N_iterations_max)) * 255.f;

        BYTE c = (BYTE)I;

        RGBQUAD color = (N_iterations[i] < N_iterations_max) ? RGBQUAD {(BYTE) (255 - c), (BYTE) (c%2 * 64), c} : RGBQUAD {};
        scr[iy][ix+i] = color;
    }

}

void count_mandelbrot(float x_shift, float y_shift, RGBQUAD scr[(size_t)window_height][(size_t)window_width], float* multipliers, float* R2_max) {


    for (size_t iy = 0; iy < (size_t)window_height; iy++) {

        if (txGetAsyncKeyState (VK_ESCAPE))
            break;

        float x0 = ((          -  (window_width / 2)) * dx + x_centre + x_shift);
        float y0 = (((float)iy - (window_height / 2)) * dy + y_centre + y_shift);

        for (size_t ix = 0; ix < (size_t)window_width; ix += points_num, x0 += points_num*dx) {

            float DX[points_num] = {};
            mm256_set1_ps(DX, dx);

            mm256_mul_ps(DX, DX, multipliers);

            float X0[points_num] = {};
            mm256_set1_ps(X0, x0);

            mm256_add_ps(X0, X0, DX);

            float X[points_num] = {};
            mm256_cpy_ps(X, X0);


            float Y0[points_num] = {};
            mm256_set1_ps(Y0, y0);

            float Y[points_num] = {};
            mm256_cpy_ps(Y, Y0);


            int N_iterations[points_num] = {};
            mm256_setzero_si256(N_iterations);

            for (size_t n_iter = 0; n_iter < N_iterations_max; n_iter++) {

                float X2[points_num] = {};
                mm256_mul_ps(X2, X, X);

                float Y2[points_num] = {};
                mm256_mul_ps(Y2, Y, Y);

                float XY[points_num] = {};
                mm256_mul_ps(XY, X, Y);


                float R2[points_num] = {};
                mm256_add_ps(R2, X2, Y2);

                float cmp[points_num] = {};
                mm256_cmple_ps(cmp, R2, R2_max);

                int mask = mm256_movemask_ps(cmp);
                if (!mask)
                    break;

                mm256_sub_epi32(N_iterations, N_iterations, cmp);

                mm256_sub_ps(X, X2, Y2);
                mm256_add_ps(X, X, X0);

                mm256_add_ps(Y, XY, XY);
                mm256_add_ps(Y, Y, Y0);

            }
            N_iter = N_iterations[0];

            #ifndef TIME_MEASUREMENT
            draw_pixels(N_iterations, scr, ix, iy);
            #endif
        }
    }
}


inline void mm256_set_ps (float res[points_num], float val7, float val6, float val5, float val4,
                                       float val3, float val2, float val1, float val0) {


    res[0] = val0;
    res[1] = val1;
    res[2] = val2;
    res[3] = val3;
    res[4] = val4;
    res[5] = val5;
    res[6] = val6;
    res[7] = val7;

}


inline void mm256_set1_ps (float res[points_num], float val) {

    for (size_t i = 0; i < points_num; i++)
        res[i] = val;
}

inline void mm256_setzero_si256 (int res[points_num]) {

    for (size_t i = 0; i < points_num; i++)
        res[i] = 0;
}

inline void mm256_cpy_ps (float res[points_num], const float src[points_num]) {

     for (size_t i = 0; i < points_num; i++)
        res[i] = src[i];

}

inline void mm256_add_ps (float res[points_num], const float reg1[points_num], const float reg2[points_num]) {

    for (size_t i = 0; i < points_num; i++)
        res[i] = reg1[i] + reg2[i];

}

inline void mm256_sub_ps (float res[points_num], const float reg1[points_num], const float reg2[points_num]) {

    for (size_t i = 0; i < points_num; i++)
        res[i] = reg1[i] - reg2[i];

}

inline void mm256_mul_ps (float res[points_num], const float reg1[points_num], const float reg2[points_num]) {

    for (size_t i = 0; i < points_num; i++)
        res[i] = reg1[i] * reg2[i];

}

inline void mm256_sub_epi32 (int res[points_num], const int reg1[points_num], const float reg2[points_num]) {

    for (size_t i = 0; i < points_num; i++)
        res[i] = reg1[i] - (int)reg2[i];

}

inline void mm256_cmple_ps (float cmp[points_num], const float reg1[points_num], const float reg2[points_num]) {

    for (size_t i = 0; i < points_num; i++) {

        if (reg1[i] < reg2[i])
            cmp[i] = -1;

        else
            cmp[i] = 0;
    }
}


inline int mm256_movemask_ps (const float cmp[points_num]) {

    int mask = 0;

    for (size_t i = 0; i < points_num; i++)
        mask |= (int)cmp[i];

    return mask;
}




















