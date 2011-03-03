/* Note that some values may need changing. */
/* Code taken from Fedora core */
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <unistd.h>

long i1mach_(long *i)
{
	switch(*i){
	  case 1:  return STDIN_FILENO;	/* standard input */
	  case 2:  return STDOUT_FILENO;	/* standard output */
	  case 3:  return 7;	/* standard punch */
	  case 4:  return STDERR_FILENO;	/* standard error */
	  case 5:  return 8*sizeof(int);	/* bits per integer */
	  case 6:  return sizeof(int);
	  case 7:  return 2;	/* base for integers */
	  case 8:  return 8*sizeof(int)-1;	/* digits of integer base */
	  case 9:  return LONG_MAX;
	  case 10: return FLT_RADIX;
	  case 11: return FLT_MANT_DIG;
	  case 12: return FLT_MIN_EXP;
	  case 13: return FLT_MAX_EXP;
	  case 14: return DBL_MANT_DIG;
	  case 15: return DBL_MIN_EXP;
	  case 16: return DBL_MAX_EXP;
	  }
	fprintf(stderr, "invalid argument: i1mach(%ld)\n", *i);
	exit(1);return 0; /* some compilers demand return values */
}
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
double d1mach_(long *i)
{
	switch(*i){
	  case 1: return DBL_MIN;
	  case 2: return DBL_MAX;
	  case 3: return DBL_EPSILON/FLT_RADIX;
	  case 4: return DBL_EPSILON;
	  case 5: return log10(FLT_RADIX);
	  }
	fprintf(stderr, "invalid argument: d1mach(%ld)\n", *i);
	exit(1); return 0; /* some compilers demand return values */
}

#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
float r1mach_(long *i)
{
	switch(*i){
	  case 1: return FLT_MIN;
	  case 2: return FLT_MAX;
	  case 3: return FLT_EPSILON/FLT_RADIX;
	  case 4: return FLT_EPSILON;
	  case 5: return log10(FLT_RADIX);
	  }
	fprintf(stderr, "invalid argument: r1mach(%ld)\n", *i);
	exit(1); return 0; /* else complaint of missing return value */
}


