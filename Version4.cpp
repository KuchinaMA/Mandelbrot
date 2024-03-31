#include "TXLib.h"
#include <immintrin.h>

//#define TIME_MEASUREMENT

static const float window_width  = 800.f;
static const float window_height = 600.f;

static const float dx = 1/window_width;
static const float dy = 1/window_width;

typedef RGBQUAD (&scr_t)[(unsigned long long)window_height][(unsigned long long)window_width];

static const float x_centre = -1.325f;
static const float y_centre = 0;

static const int N_iterations_max = 256;
static const float r2_max = 100.f;

static const __m256 N_max = _mm256_set1_ps(N_iterations_max);
static const __m256 _255 = _mm256_set1_ps(255.f);
static const __m256 multipliers = _mm256_set_ps(7.f, 6.f, 5.f, 4.f, 3.f, 2.f, 1.f, 0.f);
static const __m256 R2_max = _mm256_set1_ps(r2_max);

const size_t running_num = 100;

const size_t points_num = 8;

volatile __m256i N_iterations;

extern "C" uint64_t get_time();

void view_regulation (float* x_shift, float* y_shift, float dx, float dy);
void draw_pixels(__m256i N_iterations, RGBQUAD scr[(size_t)window_height][(size_t)window_width], size_t ix, size_t iy);
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

        view_regulation(&x_shift, &y_shift, dx, dy);

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



void view_regulation (float* x_shift, float* y_shift, float dx, float dy) {

    if (txGetAsyncKeyState (VK_RIGHT))
        *x_shift += dx * (txGetAsyncKeyState (VK_SHIFT) ? 100.f : 10.f);

    if (txGetAsyncKeyState (VK_LEFT))
        *x_shift -= dx * (txGetAsyncKeyState (VK_SHIFT) ? 100.f : 10.f);

    if (txGetAsyncKeyState (VK_UP))
        *y_shift += dy * (txGetAsyncKeyState (VK_SHIFT) ? 100.f : 10.f);

    if (txGetAsyncKeyState (VK_DOWN))
        *y_shift -= dy * (txGetAsyncKeyState (VK_SHIFT) ? 100.f : 10.f);


}

void draw_pixels(__m256i N_iterations, RGBQUAD scr[(size_t)window_height][(size_t)window_width], size_t ix, size_t iy) {

    __m256 I = _mm256_mul_ps(_mm256_sqrt_ps (_mm256_sqrt_ps (_mm256_div_ps (_mm256_cvtepi32_ps (N_iterations), N_max))), _255);

    for (size_t i = 0; i < points_num; i++) {

        int* pn = (int*) &N_iterations;
        float* pI = (float*) &I;

        BYTE c = (BYTE)pI[i];

        RGBQUAD color = (pn[i] < N_iterations_max) ? RGBQUAD {(BYTE) (255 - c), (BYTE) (c%2 * 64), c} : RGBQUAD {};
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

            __m256 X0 = _mm256_add_ps(_mm256_set1_ps(x0), _mm256_mul_ps(multipliers, _mm256_set1_ps(dx)));

            __m256 Y0 = _mm256_set1_ps(y0);

            __m256 X = X0;
            __m256 Y = Y0;

            N_iterations = _mm256_setzero_si256();

            for (int n_iter = 0; n_iter < N_iterations_max; n_iter++) {

                __m256 X2 = _mm256_mul_ps(X, X);
                __m256 Y2 = _mm256_mul_ps(Y, Y);
                __m256 XY = _mm256_mul_ps(X, Y);

                __m256 R2 = _mm256_add_ps(X2, Y2);

                __m256 cmp = _mm256_cmp_ps(R2, R2_max, _CMP_LT_OQ);

                int mask = _mm256_movemask_ps(cmp);
                if (!mask)
                    break;

                N_iterations = _mm256_sub_epi32(N_iterations, _mm256_castps_si256(cmp));

                X = _mm256_add_ps(_mm256_sub_ps(X2, Y2), X0);

                Y = _mm256_add_ps(_mm256_add_ps(XY, XY), Y0);
            }

            #ifndef TIME_MEASUREMENT

            draw_pixels(N_iterations, scr, ix, iy);

            #endif
        }
    }
}




//1 1810934   415921944
//2 1834192   187651032  1029456696 187737443
//3 1767614   259166399  282390538  276580043
//4 1963043   57402160




















