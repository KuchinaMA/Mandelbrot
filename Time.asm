; nasm    "Time.asm" -f win64 -l "Time.lst"
; g++ Time.obj TimeTest.cpp -o Time.exe
; g++ Time.obj Version1.cpp -o Version1.exe

global get_time

section .text

get_time:
            rdtsc
            shl rdx, 32
            add rax, rdx

            ret
