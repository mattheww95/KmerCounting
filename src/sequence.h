/*
 * Create an arena of sequence data for either
 * fastq data or fasta data. This will consist 
 * of a struct of charactar arrays holding the
 * individual sequnce information e.g. array of
 * pointers.
 * */


#include "zlib.h"
#include <stdio.h>
#include <stdbool.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef struct _sequence_file_data{
   ssize_t length;
   ssize_t size;
   char* data; 
}_sequence_file_data;

_sequence_file_data* init_file_data(ssize_t size){
    _sequence_file_data* data = (_sequence_file_data*)malloc(sizeof(_sequence_file_data));
   data->data = (char*)malloc(sizeof(char) * size);
   data->size = size;
   data->length = 0;
   return data;  
}


void destroy_sequence_file_data(_sequence_file_data* data){ 
    free(data->data);
    free(data);

}


void add_data(_sequence_file_data* data, char* new_data){
    size_t new_data_length = strlen(new_data);   
    if(new_data_length + data->length >= data->size){
        size_t new_arr_size = data->size * 2;
        data->data = (char*)realloc(data->data, sizeof(char) * new_arr_size);
        data->size = new_arr_size; 
    }
    memcpy(&data->data[data->length], new_data, new_data_length); 
    data->length += new_data_length;
}

/* 
 * Test if file is gzipped compressed based on magic numbers.
 * If the file cannot be read, an error is thrown and the 
 * program exits. A side effect of this function is that
 * it checks if the file exists
 *
 * filepath: Input file to decompress
 * */
bool is_gzip(const char *filepath){
    size_t error_message_size = 1000;
    unsigned char magic_bytes[2] = {0, 0};
    char error_message[error_message_size];

    FILE* file = fopen(filepath, "rb");
    if(NULL == file){
        snprintf(error_message, error_message_size, 
                "Could not open file: %s", filepath);
        perror(error_message);
        exit(EXIT_FAILURE);
    }

    fread(magic_bytes, sizeof(unsigned char), 2, file);
    fclose(file);
    if(magic_bytes[0] == 0x1f && magic_bytes[1] == 0x8b){
        return true; 
    }
    return false;
}


bool read_file(const char* filepath){
    
    //bool is_gzipped = is_gzip(filepath);
    const size_t buffer_length = 4096;
    char* buffer = (char*)malloc(sizeof(char) * buffer_length); 
    _sequence_file_data* sequence_data = init_file_data(buffer_length);
    
    gzFile file = gzopen(filepath, "r");
    if(!file){
        fprintf(stderr, "Could not open file: %s. failed: %s \n", 
                filepath, strerror(errno)); 
        exit(EXIT_FAILURE);
    }

    while(1){
        int err;
        int bytes_read = gzread(file, buffer, buffer_length - 1);
        buffer[bytes_read] = '\0';
        add_data(sequence_data, buffer);
        if(bytes_read < buffer_length - 1){
            if(gzeof(file)){break;}
            else{
                const char* error_string = gzerror(file, &err);
               if(err){
                    fprintf(stderr, "Error: %s.\n", error_string);
                    exit(EXIT_FAILURE);
               } 
            }
        }
    } 
    fprintf(stdout, "Data: ");
    for(int i = 0; i < sequence_data->length; i++){ 
        fprintf(stdout, "%c", sequence_data->data[i]); 
    }
    fprintf(stdout, "\n");

    gzclose(file);
    return true;
}



