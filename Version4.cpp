#include "TXLib.h"
#include <immintrin.h>

#define TIME_MEASUREMENT

extern "C" uint64_t get_time();
const size_t running_num = 100;

void view_regulation (float* x_shift, float* y_shift, float dx, float dy);

int main() {

    const size_t points_num = 8;

    const float x_centre = -1.325f;
    const float y_centre = 0;

    static const float width  = 800.f;
    static const float height = 600.f;

    const float dx = 1/width;
    const float dy = 1/width;

    const int N_iterations_max = 256;
    const __m256 N_max = _mm256_set1_ps(N_iterations_max);

    const float r2_max = 100.f;
    __m256 R2_max = _mm256_set1_ps(r2_max);

    float x_shift = 0.f;
    float y_shift = 0.f;

    const __m256 _255 = _mm256_set1_ps(255.f);
    const __m256 multipliers = _mm256_set_ps(7.f, 6.f, 5.f, 4.f, 3.f, 2.f, 1.f, 0.f);

    #ifndef TIME_MEASUREMENT

    txCreateWindow(width, height);
    Win32::_fpreset();
    txBegin();

    typedef RGBQUAD (&scr_t)[(unsigned long long)height][(unsigned long long)width];
    scr_t scr = (scr_t) *txVideoMemory();

    for (;;) {

        if (txGetAsyncKeyState (VK_ESCAPE))
            break;

        view_regulation(&x_shift, &y_shift, dx, dy);

        #else

        uint64_t general_time = 0;
        for (size_t run_num = 0; run_num < running_num; run_num++) {

            uint64_t start_time = get_time();
            #endif

            for (size_t iy = 0; iy < (size_t)height; iy++) {

                if (txGetAsyncKeyState (VK_ESCAPE))
                    break;

                float x0 = ((          -  (width / 2)) * dx + x_centre + x_shift);
                float y0 = (((float)iy - (height / 2)) * dy + y_centre + y_shift);

                for (size_t ix = 0; ix < (size_t)width; ix += points_num, x0 += points_num*dx) {

                    __m256 X0 = _mm256_add_ps(_mm256_set1_ps(x0), _mm256_mul_ps(multipliers, _mm256_set1_ps(dx)));

                    __m256 Y0 = _mm256_set1_ps(y0);

                    __m256 X = X0;
                    __m256 Y = Y0;

                    __m256i N_iterations = _mm256_setzero_si256();

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

                    __m256 I = _mm256_mul_ps(_mm256_sqrt_ps (_mm256_sqrt_ps (_mm256_div_ps (_mm256_cvtepi32_ps (N_iterations), N_max))), _255);

                    for (size_t i = 0; i < points_num; i++) {

                        int* pn = (int*) &N_iterations;
                        float* pI = (float*) &I;

                        BYTE c = (BYTE)pI[i];

                        RGBQUAD color = (pn[i] < N_iterations_max) ? RGBQUAD {(BYTE) (255 - c), (BYTE) (c%2 * 64), c} : RGBQUAD {};
                        scr[iy][ix+i] = color;

                    }
                    #endif
                }
            }
            #ifdef TIME_MEASUREMENT
            uint64_t time_taken = get_time() - start_time;
            general_time += time_taken;
        }

        printf("%llu\n", general_time/running_num);
        #else

        printf("\t\r%.0lf", txGetFPS());
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




//1 1810934
//2 1886058
//3 1888461
//4 1961616




















