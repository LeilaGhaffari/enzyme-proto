#include <stdio.h>

int main() {
  printf("hello world\n");
  return 0;
}

/*
  $ clang hello.c -o hello
  $ ./hello

  -- lli directly executes programs in LLVM bitcode format. 
  It takes a program in LLVM bitcode format and executes it 
  using a just-in-time compiler or an interpreter.
  lli is not an emulator. It will not execute IR of different 
  architectures and it can only interpret (or JIT-compile) for 
  the host architecture.

  $ clang -O3 -emit-llvm hello.c -c -o hello.bc
  $ lli hello.bc

  $ clang -O3 -emit-llvm hello.c -S -o hello.ll
  $ lli hello.ll

  LLC - LLVM static compiler
  The llc command compiles LLVM source inputs into assembly language 
  for a specified architecture. The assembly language output can 
  then be passed through a native assembler and linker to generate 
  a native executable.

  $ llc hello.bc -o hello.s

  Assemble the native assembly language file into a program
  $ gcc hello.s -o hello.native
  $ ./hello.native
*/
