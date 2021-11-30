#include <stdio.h>

extern double __enzyme_autodiff(void *, double);

double square(double x) {
    return x * x;
}

double dsquare(double x) {
    // This returns the derivative of square or 2 * x
    return __enzyme_autodiff((void *) square, x);
}

int main() {
    for(double i=1; i<5; i++)
        printf("%f, square(%f)=%f, dsquare(%f)=%f \n", i, i, square(i), i, dsquare(i));
}

/*
  We can generate LLVM from this code by calling clang as follows. 
  Note that clang should be the path to whatever clang you built Enzyme against.

    clang test001.c -S -emit-llvm -o input.ll -O2 -fno-vectorize -fno-slp-vectorize -fno-unroll-loops

  where
    -S -emit-llvm   -> we want to emit LLVM bitcode rather than an executable
    -o input.ll     -> we want the output to be in a file input.ll
    -O2 -ffast-math -> runs optimizations (with fast-math) before we run Enzyme’s AD process
    -fno-vectorize -fno-slp-vectorize -fno-unroll-loops ->
         we don’t want to run vectorization or loop unrolling. 
         It is better for performance to only run these scheduling optimizations after AD.
  Perfoming AD Enzyme:
    opt input.ll -load=/path/to/Enzyme/enzyme/build/Enzyme/LLVMEnzyme-<VERSION>.so -enzyme -o output.ll -S
    
*/

/*

clang test001.c -S -emit-llvm -o input1.ll -O2 -fno-vectorize -fno-slp-vectorize -fno-unroll-loops

opt input1.ll -load=/home/leila/Enzyme/enzyme/build12DHB/Enzyme/ClangEnzyme-12.so -enzyme -o output1.ll -S

opt output1.ll -O2 -o output_opt1.ll -S

clang output_opt1.ll -o a1.exe

./a1.exe

*/

// OR EVEN EASIER WITH ALL OPTIMIZATIONS
/*

clang test001.c -Xclang -load -Xclang /home/linuxbrew/.linuxbrew/Cellar/enzyme/0.0.19/lib/ClangEnzyme-12.so -O2 -fno-vectorize -fno-unroll-loops; ./a.out

*/
