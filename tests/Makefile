
# input.ll: test034.c
# 	clang-13 test034.c -S -emit-llvm -o input.ll -O2 -fno-vectorize -fno-slp-vectorize -fno-unroll-loops
# 
# output.ll: input.ll
# 	/usr/lib/llvm-13/bin/opt input.ll -load=/home/leila/Enzyme/enzyme/buildglit/Enzyme/LLVMEnzyme-13.so -enzyme -o output.ll -S
# 
# main: output.ll
# 	clang-13 output.ll -o main

main: test034.c
	clang-13 -o main test034.c -Xclang -load -Xclang /home/leila/Enzyme/enzyme/buildglit/Enzyme/ClangEnzyme-13.so -O2 -flegacy-pass-manager 

# clang-13 -o main test034.c -Xclang -load -Xclang /home/leila/Enzyme/enzyme/buildglit/Enzyme/ClangEnzyme-13.so -O2 -flegacy-pass-manager 

main2: test034.c
	clang-13 -fuse-ld=lld -flto -o main2 test034.c -Wl,--lto-legacy-pass-manager -Wl,-mllvm=-load=/home/leila/Enzyme/enzyme/buildglit/Enzyme/LLDEnzyme-13.so
# clang-13 -std=c11 -O2 %s -S -emit-llvm -o - %loadClangEnzyme | %lli -