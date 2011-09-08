#include <stdio.h>
#include <stdint.h>


/*
 * dummy-code.c
 *
 *  Created on: Aug 28, 2011
 *      Author: Niels Drost
 */

void arpc_set_input_buffers(int *header, int *integers, int *longs) {
	fprintf(stderr, "dummy code: input buffers set\n");

}

//sets buffers to be used for output
void arpc_set_output_buffers(int32_t *header, int32_t *integers) {
	fprintf(stderr, "dummy code: output buffer set\n");
}

//actual call to code. expects the code to take the input from the input
//buffers, and return the output (or any error) in the output buffers
void arpc_call() {
	fprintf(stderr, "dummy code: code called\n");
}
