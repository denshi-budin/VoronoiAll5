#include <stdio.h>
#include <stdlib.h>
#include "atom.h"

#ifdef WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <sys/time.h>
#endif

#ifdef WIN32
double PCFreq = 1000.0;
__int64 timerStart = 0;
#else
struct timeval timerStart;
#endif

void StartTimer()
{
	#ifdef WIN32
	LARGE_INTEGER li;
	if (!QueryPerformanceFrequency(&li))
		throw "QueryPerformanceFrequency failed!";

	PCFreq = (double)li.QuadPart / 1000.0;

	QueryPerformanceCounter(&li);
	timerStart = li.QuadPart;
	#else
	    gettimeofday(&timerStart, NULL);
	#endif
}

// time elapsed in ms
double GetTimer()
{
	#ifdef WIN32
	LARGE_INTEGER li;
	QueryPerformanceCounter(&li);
	return (double)(li.QuadPart - timerStart) / PCFreq;
	#else
	struct timeval timerStop, timerElapsed;
	gettimeofday(&timerStop, NULL);
	timersub(&timerStop, &timerStart, &timerElapsed);
	return timerElapsed.tv_sec*1000.0+timerElapsed.tv_usec/1000.0;
	#endif
}

void PrintTimer(double timer)
{
	int mins = 0, hour = 0, day = 0;
    timer /= 1000;
    if (timer >= 60){
        mins = (int)(timer / 60);
        timer = timer - (mins * 60);
        if (mins >= 60){
            hour = mins / 60;
            mins = mins % 60;
            if (hour >= 24){
                day = hour / 24;
                hour = hour % 24;
                printf("%d Days %d Hours %d Minutes %g Seconds\n", day, hour, mins, timer);
            }
            else printf("%d Hours %d Minutes %g Seconds\n", hour, mins, timer);
        }
        else printf("%d Minutes %g Seconds\n", mins, timer);
    }
    else printf("%g Seconds\n", timer);
}

#endif
