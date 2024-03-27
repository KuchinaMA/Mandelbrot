#include <stdio.h>
#include "TXLib.h"

#define TIME_MEASUREMENT

extern "C" uint64_t get_time();

int main() {

    const float x_centre = -1.325f;
    const float y_centre = 0;

    static const float width  = 800.f;
    static const float height = 600.f;

    const float dx = 1/width;
    const float dy = 1/width;

    const int N_iterations_max = 256;
    const float R2_max = 100.f;

    float x_shift = 0.f;
    float y_shift = 0.f;

    const size_t points_num = 8;

    #ifndef TIME_MEASUREMENT

    txCreateWindow(width, height);
    Win32::_fpreset();
    txBegin();

    typedef RGBQUAD (&scr_t)[(unsigned long long)height][(unsigned long long)width];
    scr_t scr = (scr_t) *txVideoMemory();

    for (;;) {

        if (txGetAsyncKeyState (VK_ESCAPE))
            break;

        if (txGetAsyncKeyState (VK_RIGHT))
            x_shift += dx * (txGetAsyncKeyState (VK_SHIFT) ? 100.f : 10.f);

        if (txGetAsyncKeyState (VK_LEFT))
            x_shift -= dx * (txGetAsyncKeyState (VK_SHIFT) ? 100.f : 10.f);

        if (txGetAsyncKeyState (VK_UP))
            y_shift += dy * (txGetAsyncKeyState (VK_SHIFT) ? 100.f : 10.f);

        if (txGetAsyncKeyState (VK_DOWN))
            y_shift -= dy * (txGetAsyncKeyState (VK_SHIFT) ? 100.f : 10.f);

        /*if (txGetAsyncKeyState ('A'))
            scale += dx * (txGetAsyncKeyState (VK_SHIFT) ? 100.f : 10.f);

        if (txGetAsyncKeyState ('Z'))
            scale -= dx * (txGetAsyncKeyState (VK_SHIFT) ? 100.f : 10.f);*/


        #else
        uint64_t start_time = get_time();
        #endif

        for (size_t iy = 0; iy < (size_t)height; iy++) {

            if (txGetAsyncKeyState (VK_ESCAPE))
                break;

            /*float x0 = ((          -  (width / 2)) * dx + ROI_X + xC) * scale,
                  y0 = (((float)iy - (height / 2)) * dy + ROI_Y + yC) * scale;*/
            float x0 = ((          -  (width / 2)) * dx + x_centre + x_shift);
            float y0 = (((float)iy - (height / 2)) * dy + y_centre + y_shift);

            for (size_t ix = 0; ix < (size_t)width; ix += points_num, x0 += points_num*dx) {

                float X0[points_num] = {x0, x0 + dx, x0 + 2*dx, x0 + 3*dx, x0 + 4*dx, x0 + 5*dx, x0 + 6*dx, x0 + 7*dx},
                      Y0[points_num] = {y0, y0, y0, y0, y0, y0, y0, y0};

                float X[points_num] = {}; for (size_t i = 0; i < points_num; i++) X[i] = X0[i];
                float Y[points_num] = {}; for (size_t i = 0; i < points_num; i++) Y[i] = Y0[i];

                float N_iterations[8] = {0, 0, 0, 0, 0, 0, 0, 0};

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

                #ifndef TIME_MEASUREMENT

                for (size_t i = 0; i < points_num; i++) {

                    float I = sqrtf (sqrtf ((float)N_iterations[i] / (float)N_iterations_max)) * 255.f;

                    BYTE c = (BYTE)I;

                    RGBQUAD color = (N_iterations[i] < N_iterations_max) ? RGBQUAD {(BYTE) (255 - c), (BYTE) (c%2 * 64), c} : RGBQUAD {};
                    scr[iy][ix+i] = color;
                }
                #endif

            }
        }
        #ifdef TIME_MEASUREMENT

        uint64_t end_time = get_time() - start_time;
        printf("%llu\n", end_time);

        #else

        printf("\t\r%.0lf", txGetFPS());
        txSleep();
    }
    #endif

}
