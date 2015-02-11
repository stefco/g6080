/*  2/10/15  Stefan Countryman
 *
 *  A program for timing various looped operations on Unix-like machines.
 *
 *  Most timing functionality is copied straight from RDM's real_time_clock.
 *
 *  The program will print the results in CSV format. The user should run the
 *  program and append each new line to the file timing-results-o0.csv if
 *  compiled with -O0 optimization or timing-results-o3.csv if compiled with
 *  -O3 optimization.
 *
 *  The Mac OS X implementation came from gist/github.com/jbenet/1087739
 *  This code uses the real time clock feature.  These clocks return times in
 *  nanoseconds.
 *
 *  On Linux, compile with
 *    gcc -o real_time_clock real_time_clock.c -lrt
 *
 *  (The -lrt flag links in the real-time libraries.)

 *  On Mac, compile with
 *    gcc -o real_time_clock real_time_clock.c
 *
 */

/* 
author: jbenet
os x, compile with: gcc -o testo test.c 
linux, compile with: gcc -o testo test.c -lrt
*/
 
#include <time.h>
#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif
 
 
void current_utc_time(struct timespec *ts) {
 
#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
  clock_serv_t cclock;
  mach_timespec_t mts;
  host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
  clock_get_time(cclock, &mts);
  mach_port_deallocate(mach_task_self(), cclock);
  ts->tv_sec = mts.tv_sec;
  ts->tv_nsec = mts.tv_nsec;
#else
  clock_gettime(CLOCK_REALTIME, ts);
#endif
 
}

/*  The size of the a & b arrays, for section (c) */

const int absize = 20020;

/*  The number of operations in each loop */

const int nops = 10000;

/*  The number of results to print out */

const int nprint = 8;

/*  The seed to start the random number sequence */

const int seed = 1234;

int main()
{

  double a[absize], b[absize], c[nops];
  int i;

  struct timespec tstart, tstop;

  struct timespec ts;
  current_utc_time(&ts);

  /*  Initialize random floating point values */

  srand48(seed);

  for ( i = 0; i < nops; i++ ) {

    a[i] = fabs(drand48()); // make a[i]>0 for (h)
    b[i] = drand48();

  }

  /* Initialize timing variables */
  
  tstart.tv_sec = 0;
  tstart.tv_nsec = 0;
  tstop.tv_sec = 0;
  tstop.tv_nsec = 0;

  /* Time each operation */

  // a
  current_utc_time(&tstart);

  for ( i = 0; i < nops; i++ ) {

    c[i] = a[i] * b[i]; 

  }

  current_utc_time(&tstop);

  printf("%lu,", tstop.tv_nsec - tstart.tv_nsec);

  // b
  current_utc_time(&tstart);

  for ( i = 0; i < nops; i++ ) {

    c[i] += a[i] * b[i]; 

  }

  current_utc_time(&tstop);

  printf("%lu,", tstop.tv_nsec - tstart.tv_nsec);

  // c
  current_utc_time(&tstart);

  for ( i = 0; i < nops; i++ ) {

    c[i] += a[i] * b[2*i + 20]; 

  }

  current_utc_time(&tstop);

  printf("%lu,", tstop.tv_nsec - tstart.tv_nsec);

  // d
  current_utc_time(&tstart);

  for ( i = 0; i < nops; i++ ) {

    c[i] = a[i] / b[i]; 

  }

  current_utc_time(&tstop);

  printf("%lu,", tstop.tv_nsec - tstart.tv_nsec);

  // e
  current_utc_time(&tstart);

  for ( i = 0; i < nops; i++ ) {

    c[i] += a[i] / b[i]; 

  }

  current_utc_time(&tstop);

  printf("%lu,", tstop.tv_nsec - tstart.tv_nsec);

  // f
  current_utc_time(&tstart);

  for ( i = 0; i < nops; i++ ) {

    c[i] = sin(a[i]); 

  }

  current_utc_time(&tstop);

  printf("%lu,", tstop.tv_nsec - tstart.tv_nsec);

  // g
  current_utc_time(&tstart);

  for ( i = 0; i < nops; i++ ) {

    c[i] = exp(a[i]); 

  }

  current_utc_time(&tstop);

  printf("%lu,", tstop.tv_nsec - tstart.tv_nsec);
  
  // h
  current_utc_time(&tstart);

  for ( i = 0; i < nops; i++ ) {

    c[i] = sqrt(a[i]);

  }

  current_utc_time(&tstop);

  /* Replace final comma with newline character */
  printf("%lu\n", tstop.tv_nsec - tstart.tv_nsec);

  return 0;
}
