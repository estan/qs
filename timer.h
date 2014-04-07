#ifndef TIMER_H
#define TIMER_H

#ifdef ENABLE_TIMER
    /*
     * If ENABLE_TIMER is defined we define a clock and two macros START and
     * STOP. Use START() to start a timer and STOP("some message") to stop
     * it and print the time elapsed since START was called in ms to stdout.
     */
    clock_t timer;
    #define START() timer = std::clock();
    #define STOP(msg) \
        std::cout << msg << " in ";\
        std::cout << (1000.0 * (std::clock() - timer)/CLOCKS_PER_SEC);\
        std::cout << " ms" << std::endl;
#else
    // Else the macros are defined as no-ops.
    #define START(x)
    #define STOP(x)
#endif // ENABLE_TIMER

#endif // TIMER_H
