/* Return in array double precision system parameters */
#include <float.h>

void mfloat(values)
float values[];
{
	values[0] = FLT_EPSILON;
	values[1] = FLT_MIN;
	values[2] = FLT_MAX;
	values[3] = FLT_MIN_10_EXP;
	values[4] = FLT_MAX_10_EXP;
}
void mdble_(values)
double values[];
{
	values[0] = DBL_EPSILON;
	values[1] = DBL_MIN;
	values[2] = DBL_MAX;
	values[3] = DBL_MIN_10_EXP;
	values[4] = DBL_MAX_10_EXP;
}
void mdble(values)
double values[];
{ mdble_(values);
}
