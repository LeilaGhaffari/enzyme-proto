// freeing memory

#include <stdio.h>
#include <stdlib.h>


int main() {

    void *tape;
    void *tape_;
    tape = malloc(8);
    tape_ = tape; 
    free(tape_);
    // free(tape); This would be double freeing

    return 0;
}

/*
clang test015.c

valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose --log-file=valgrind-out.txt ./a.out 
*/
