#include "src/PACO.h"

static float upsilon = 0.95;

static char * pattern_buffer = 0;
static float * autocorrelate_data = 0;
static float * delta = 0;
static int delta_length = 0;
static int length_of_the_pattern = 0;
static float psratio = 0;
static char * repeating_pattern = 0;

void cleanup_buffers()
{
    if(pattern_buffer != 0) {
        delete[] pattern_buffer;
        pattern_buffer = 0;
    }
    if(autocorrelate_data != 0) {
        delete[] autocorrelate_data;
        autocorrelate_data = 0;
    }
    if(delta != 0) {
        delete[] delta;
        delta = 0;
    }
    if(repeating_pattern != 0) {
        delete[] repeating_pattern;
        repeating_pattern = 0;
    }
}

int determine_pattern(double * x, double * y,  double * vx, double * vy, int number_of_points)
{
    length_of_the_pattern = 0;
    pattern_buffer = new char[number_of_points];
    
    PatternConstruct( x, y, vx, vy, number_of_points, pattern_buffer, &length_of_the_pattern) ;
    pattern_buffer[length_of_the_pattern] = 0;
    
    autocorrelate_data = new float[length_of_the_pattern];
    
	int length_of_repeating_pattern = AutoCorrelate(pattern_buffer, autocorrelate_data, length_of_the_pattern, upsilon, 1);
    
    delta_length = length_of_the_pattern - length_of_repeating_pattern;
    delta = new float[delta_length];
    Delta(delta, delta_length, pattern_buffer,length_of_the_pattern,length_of_repeating_pattern);
    
    psratio = PSratio(delta,length_of_the_pattern - 2 * length_of_repeating_pattern); 
    	
	int check = 0;
    int i = 0; 
    
    while ( check < length_of_repeating_pattern && i < length_of_the_pattern - length_of_repeating_pattern ) {
        check++;
        check *= delta[i];
        i++;
    } 
    
	if( i < length_of_the_pattern - length_of_repeating_pattern && i >= length_of_repeating_pattern ) { 
        repeating_pattern = new char[length_of_repeating_pattern + 1];
        for(int j = 0; j < length_of_repeating_pattern; j++) {
            repeating_pattern[j] = pattern_buffer[j + (i - length_of_repeating_pattern)];
        }
        repeating_pattern[length_of_repeating_pattern] = 0;
        return 0;
    } 
	else {
        return -1;
    }
}

int get_upsilon(float * value)
{
    *value = upsilon;
    return 0;
}

int set_upsilon(float value)
{
    upsilon = value;
    return 0;
}

int get_psratio(float * output) 
{
    *output = psratio;
    return 0;
}

static char * empty_string = "";

int get_repeating_pattern(char ** output)
{
   
    if(repeating_pattern == 0) {
        *output = empty_string;
    } else {
        *output = repeating_pattern;
    }
    return 0;
}

int get_pattern_buffer(char ** output)
{ 
    if(pattern_buffer == 0) {
        *output = empty_string;
    } else {
        *output = pattern_buffer;
    }
    return 0;
}

int get_delta_length(int * output)
{
    *output = delta_length;
    return 0;
}

int get_delta(int * output, int n)
{
    if( n > delta_length) {
        return -1;
    }
    for(int i = 0; i < n; i++) {
        output[i] = delta[i];
    }
    return 0;
}

int get_length_of_the_pattern(int * output)
{
    *output = length_of_the_pattern;
    return 0;
}


int get_autocorrelate_data(float * output, int n)
{
    if( n > length_of_the_pattern) {
        return -1;
    }
    for(int i = 0; i < n; i++) {
        output[i] = autocorrelate_data[i];
    }
    return 0;
}
