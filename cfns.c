#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "mkl.h"

MKL_LONG status = 0;
DFTI_DESCRIPTOR_HANDLE FORWARD_DH, BACKWARD_DH, FORWARD_PADDED_DH, BACKWARD_PADDED_DH;

void create_descriptor_handles(int Nx, int Ny) {
    MKL_LONG N[2] = {Nx, Ny};
    MKL_LONG M[2] = {3*Nx/2, 3*Ny/2};
    MKL_LONG R_strides[3], C_strides[3];                        // Strides for real and Fourier space arrays
    MKL_LONG R_strides_pad[3], C_strides_pad[3];                // Strides for real and Fourier space padded arrays

    /* Set array sizes and strides */
    R_strides[0] = 0; R_strides[1] = 1; R_strides[2] = Nx;
    C_strides[0] = 0; C_strides[1] = 1; C_strides[2] = Nx;
    R_strides_pad[0] = 0; R_strides_pad[1] = 1; R_strides_pad[2] = 3*Nx/2;
    C_strides_pad[0] = 0; C_strides_pad[1] = 1; C_strides_pad[2] = 3*Nx/2;

    /* Create descriptor handles used for transforms */
    DftiCreateDescriptor(&FORWARD_DH, DFTI_DOUBLE, DFTI_REAL, 2, N);
    DftiSetValue(FORWARD_DH, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
    DftiSetValue(FORWARD_DH, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    DftiSetValue(FORWARD_DH, DFTI_INPUT_STRIDES, R_strides);
    DftiSetValue(FORWARD_DH, DFTI_OUTPUT_STRIDES, C_strides);
    DftiSetValue(FORWARD_DH, DFTI_FORWARD_SCALE, 1./(Nx*Ny));
    DftiCommitDescriptor(FORWARD_DH);

    DftiCreateDescriptor(&BACKWARD_DH, DFTI_DOUBLE, DFTI_REAL, 2, N);
    DftiSetValue(BACKWARD_DH, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
    DftiSetValue(BACKWARD_DH, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    DftiSetValue(BACKWARD_DH, DFTI_INPUT_STRIDES, C_strides);
    DftiSetValue(BACKWARD_DH, DFTI_OUTPUT_STRIDES, R_strides);
    DftiCommitDescriptor(BACKWARD_DH);

    DftiCreateDescriptor(&FORWARD_PADDED_DH, DFTI_DOUBLE, DFTI_REAL, 2, M);
    DftiSetValue(FORWARD_PADDED_DH, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
    DftiSetValue(FORWARD_PADDED_DH, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    DftiSetValue(FORWARD_PADDED_DH, DFTI_INPUT_STRIDES, R_strides_pad);
    DftiSetValue(FORWARD_PADDED_DH, DFTI_OUTPUT_STRIDES, C_strides_pad);
    DftiSetValue(FORWARD_PADDED_DH, DFTI_FORWARD_SCALE, 1./(M[0]*M[1]));
    DftiCommitDescriptor(FORWARD_PADDED_DH);

    DftiCreateDescriptor(&BACKWARD_PADDED_DH, DFTI_DOUBLE, DFTI_REAL, 2, M);
    DftiSetValue(BACKWARD_PADDED_DH, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
    DftiSetValue(BACKWARD_PADDED_DH, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    DftiSetValue(BACKWARD_PADDED_DH, DFTI_INPUT_STRIDES, C_strides_pad);
    DftiSetValue(BACKWARD_PADDED_DH, DFTI_OUTPUT_STRIDES, R_strides_pad);
    DftiCommitDescriptor(BACKWARD_PADDED_DH);
}

void free_descriptor_handles() {
    DftiFreeDescriptor(&FORWARD_DH);
    DftiFreeDescriptor(&BACKWARD_DH);
    DftiFreeDescriptor(&FORWARD_PADDED_DH);
    DftiFreeDescriptor(&BACKWARD_PADDED_DH);
}

void fft_forward(double *in, MKL_Complex16 *out) {
    status = DftiComputeForward(FORWARD_DH, in, out);
}

void fft_backward(MKL_Complex16 *in, double *out) {
    status = DftiComputeBackward(BACKWARD_DH, in, out);
}

void fft_padded_forward(double *in, MKL_Complex16 *out) {
    status = DftiComputeForward(FORWARD_PADDED_DH, in, out);
}

void fft_padded_backward(MKL_Complex16 *in, double *out) {
    status = DftiComputeBackward(BACKWARD_PADDED_DH, in, out);
}