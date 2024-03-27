#include <stdio.h>
#include <stdint.h>
#include "TXLib.h"

#define TIME_MEASUREMENT

extern "C" uint64_t get_time();
const size_t running_num = 500;

int main() {

    const float x_centre = -1.325f;
    const float y_centre = 0;

    static const float width  = 800.f;
    static const float height = 600.f;

    const float dx = 1/width;
    const float dy = 1/width;

    const int N_iterations_max = 256;
    const float r2_max = 100.f;

    float x_shift = 0.f;
    float y_shift = 0.f;
    //float scale = 1.f;

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
            scale -= dx * (txGetAsyncKeyState (VK_SHIFT) ? 100.f : 10.f);   */
        #else
        //uint64_t start_time = get_time();
        uint64_t general_time = 0;
        for (size_t run_num = 0; run_num < running_num; run_num++) {

            uint64_t start_time = get_time();
            #endif

            for (size_t iy = 0; iy < (size_t)height; iy++) {

                if (txGetAsyncKeyState (VK_ESCAPE))
                    break;

                /*float x0 = ((          -  (width / 2)) * dx + ROI_X + xC) * scale,
                        y0 = (((float)iy - (height / 2)) * dy + ROI_Y + yC) * scale;*/
                float x0 = ((          -  (width / 2)) * dx + x_centre + x_shift);
                float y0 = (((float)iy - (height / 2)) * dy + y_centre + y_shift);

                for (size_t ix = 0; ix < (size_t)width; ix++, x0 += dx) {

                    float x = x0,
                          y = y0;

                    int N_iterations = 0;
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
                float I = sqrtf (sqrtf ((float)N_iterations / (float)N_iterations_max)) * 255.f;
                //float I = (N % 2) * 255.f;

                BYTE c = (BYTE)I;

                RGBQUAD color = (N_iterations < N_iterations_max) ? RGBQUAD {(BYTE) (255 - c), (BYTE) (c%2 * 64), c} : RGBQUAD {};
                //RGBQUAD color = (N < nMax) ? RGBQUAD {(BYTE) c*6, 0, c*10} : RGBQUAD {};
                //RGBQUAD color = (N < nMax) ? RGBQUAD {c, c, c} : RGBQUAD {};
                scr[iy][ix] = color;
                #endif

                }
            }
            #ifdef TIME_MEASUREMENT
            uint64_t time_taken = get_time() - start_time;
            general_time += time_taken;

        }

        //#ifdef TIME_MEASUREMENT

        //uint64_t end_time = get_time() - start_time;
        //printf("%llu\n", end_time);
        printf("%llu\n", general_time/running_num);

        #else
        printf("\t\r%.0lf", txGetFPS());
        txSleep();
    }
    #endif

}
