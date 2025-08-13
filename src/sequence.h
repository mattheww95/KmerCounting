/*
 * Create an arena of sequence data for either
 * fastq data or fasta data. This will consist 
 * of a struct of charactar arrays holding the
 * individual sequnce information e.g. array of
 * pointers.
 * */


//#include <zlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <errno.h>
#include <stdlib.h>




bool is_gzip(const char *filepath){
    size_t error_message_size = 1000;
    FILE *file = NULL;
    unsigned char magic_bytes[2] = {0, 0};
    char error_message[error_message_size];

    file = fopen(filepath, "rb");
    if(NULL == file){
        snprintf(error_message, error_message_size, "Could not open file: %s", filepath);
        perror(error_message);
        exit(EXIT_FAILURE);
    }

    fread(magic_bytes, sizeof(unsigned char), 2, file);
    if(magic_bytes[0] == 0x1f && magic_bytes[1] == 0x8b){
        return true; 
    }
    return false;

}



