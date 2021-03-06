#include <stdio.h>

// We move the body of collide into a separate device function collide_body to allow us
// to pass collide_body to various differentiation methods. This is necessary as differentiation
// can only be done on device, not global kernel functions.
__device__
void collide_body(float* src, float* dst) {
    size_t idx = threadIdx.x;
    if (idx < 100) {
        dst[idx] += src[idx] * src[idx] - 3 * src[idx];
    }
}

// GPU Kernel
__global__
void collide(float* src, float* dst) {
    collide_body(src, dst);
}

// Wrapper CPU function which calls kernel
__attribute__((noinline))
void kern(float* src, float* dst) {
    collide<<<1, 100>>>(src, dst);
}

// // Main CPU code that calls wrapper function
// void iter(int nTimeSteps, float* src, float* dst) {
//     for (unsigned int i=0; i<nTimeSteps/2; i++) {
//         kern(src, dst);
//         kern(dst, src);
//     }
// }
// 
// template <typename... Args>
// void __enzyme_autodiff(Args...);
// 
// void grad_iter(int nTimeSteps, float* src, float* dsrc, float* dst, float* ddst) {
//     __enzyme_autodiff(iter, nTimeSteps, src, dsrc, dst, ddst);
// }

// A function similar to __enzyme_autodiff, except it only calls the augmented forward pass, returning
// a tape structure to hold any values that may be overwritten and needed for the reverse.
template <typename... Args>
__device__ void* __enzyme_augmentfwd(Args...);

// A function similar to __enzyme_autodiff, except it only calls the revese pass, taking in the tape
// as its last argument.
template <typename... Args>
__device__ void __enzyme_reverse(Args...);

// A wrapper GPU kernel for calling the forward pass of collide. The wrapper code stores
// the tape generated by Enzyme into a unique location per thread
__global__ void aug_collide(float* src, float* dsrc, float* dst, float* ddst, void** tape)
{
    size_t idx = threadIdx.x;
    tape[idx] = __enzyme_augmentfwd((void*)collide_body, src, dsrc, dst, ddst);
}

// A wrapper GPU kernel for calling the reverse pass of collide. The wrapper code retrieves
// the corresponding tape per thread being executed.
__global__ void rev_collide( float* src, float* dsrc, float* dst, float* ddst, void** tape)
{
    size_t idx = threadIdx.x;
    __enzyme_reverse((void*)collide_body, src, dsrc, dst, ddst, tape[idx]);
}

// The augmented forward pass of the CPU kern call, allocating and returning
// tape memory  needed to compute the reverse pass. This calls a augmented collide
// GPU kernel, passing in a unique 8-byte location to store the tape.
void* aug_kern(float* src, float* dsrc, float* dst, float* ddst) {
    void** tape;
    cudaMalloc(&tape, sizeof(void*) * /*total number of threads*/100);
    aug_collide<<<1, 100>>>(src, dsrc, dst, ddst, tape);
    return (void*)tape;
}

// The reverse pass of the CPU kern call, using tape memory passed as the
// last argument. This calls a reverse collide GPU kernel.
void rev_kern(float* src, float* dsrc, float* dst, float* ddst, void* tape) {
    rev_collide<<<1, 100>>>(src, dsrc, dst, ddst, (void**)tape);
    cudaFree(tape);
}

// Here we register the custom forward pass aug_kern and reverse pass rev_kern
void* __enzyme_register_gradient_kern[3] = { (void*)kern, (void*)aug_kern, (void*)rev_kern };

int main() {

    double *x, *d_x, *y, *d_y; // device pointers

    cudaMalloc(&x, sizeof(*x));
    cudaMalloc(&d_x, sizeof(*d_x));
    cudaMalloc(&y, sizeof(*y));
    cudaMalloc(&d_y, sizeof(*d_y));

    double host_x = 1.4;
    double host_d_x = 0.0;
    double host_y;
    double host_d_y = 1.0;

    cudaMemcpy(x,   &host_x,   sizeof(*x),   cudaMemcpyHostToDevice);
    cudaMemcpy(d_x, &host_d_x, sizeof(*d_x), cudaMemcpyHostToDevice);
    cudaMemcpy(y,   &host_y,   sizeof(*y),   cudaMemcpyHostToDevice);
    cudaMemcpy(d_y, &host_d_y, sizeof(*d_y), cudaMemcpyHostToDevice);

    // ToDo

    cudaDeviceSynchronize(); // synchroniz

    cudaMemcpy(&host_x,   x,   sizeof(*x),   cudaMemcpyDeviceToHost);
    cudaMemcpy(&host_d_x, d_x, sizeof(*d_x), cudaMemcpyDeviceToHost);
    cudaMemcpy(&host_y,   y,   sizeof(*y),   cudaMemcpyDeviceToHost);
    cudaMemcpy(&host_d_y, d_y, sizeof(*d_y), cudaMemcpyDeviceToHost);

    printf("%f %f\n", host_x,   host_y);
    printf("%f %f\n", host_d_x, host_d_y);
}
