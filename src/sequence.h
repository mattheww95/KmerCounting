/*
 * Create an arena of sequence data for either
 * fastq data or fasta data. This will consist 
 * of a struct of charactar arrays holding the
 * individual sequnce information e.g. array of
 * pointers.
 * */


#include "zlib.h"
#include "kseq.h"
#include <stdio.h>
#include <stdbool.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>


KSEQ_INIT(gzFile, gzread);

typedef struct SeqData{
    size_t size;
    size_t length;
    char* sequence;
}SeqData;


void add_sequence_data(SeqData* data, char* new_data, size_t data_length){
    if(data_length + data->length >= data->size){
        size_t new_arr_size = data->size * 2;
        data->sequence = (char*)realloc(data->sequence, sizeof(char) * new_arr_size);
        data->size = new_arr_size; 
    }
    memcpy(&data->sequence[data->length], new_data, data_length); 
    data->length += data_length;
}

SeqData* get_data(const char* filepath){
    
    gzFile fp;
    kseq_t *seq;
    int l = 0;
    fp = gzopen(filepath, "r");
    seq = kseq_init(fp);
    size_t seq_data_start = 32768; 
    SeqData* seq_data = (SeqData*)malloc(sizeof(SeqData));
    seq_data->sequence = (char*)malloc(sizeof(char) * seq_data_start);
    seq_data->length = 0;
    seq_data->size = seq_data_start;
    
    while((l = kseq_read(seq)) >= 0){
        add_sequence_data(seq_data, seq->seq.s, seq->seq.l); 
    }
    kseq_destroy(seq);
    gzclose(fp);
    return seq_data;
}



// Abandoning the below work for now
  
typedef enum DataType{
    FASTQ,
    FASTA
}DataType;


/* 
typedef struct CaratVector{
    size_t items;
    size_t size;
    SequenceIndex* locations;

}CaratVector;
*/

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


void add_data(_sequence_file_data* data, char* new_data, size_t data_length){
    if(data_length + data->length >= data->size){
        size_t new_arr_size = data->size * 2;
        data->data = (char*)realloc(data->data, sizeof(char) * new_arr_size);
        data->size = new_arr_size; 
    }
    memcpy(&data->data[data->length], new_data, data_length); 
    data->length += data_length;
}


/*  
CaratVector* find_fasta_sequences(const char* restrict buffer){

     size_t buffer_length = strlen(buffer);  
     CaratVector* vector = (CaratVector*)malloc(sizeof(CaratVector));
     size_t initial_size = 200;
     vector->locations = (SequenceIndex*)malloc(sizeof(SequenceIndex) * initial_size); 
     vector->size = buffer_length;
     vector->items = initial_size;
     bool in_sequence = false;
     size_t header_start = 0; 
     size_t sequence_start = 0; 
     size_t sequence_end = 0; 
     for(size_t i = 0; i < buffer_length; i++){
         char test = buffer[i];
         if(test == '>'){
            header_start = i;
         }else if(!in_sequence && test == '\n'){
             i++;
             in_sequence = true;
             sequence_start = i; // Next char is sequence 
         }else if(in_sequence && test == '>'){
           sequence_end = i--;  
           vector->locations[vector->items] = SequenceIndex{header_start, 
               sequence_start, sequence_end};
           vector->items++; 
           in_sequence = false;
           if(vector->items >= vector->size){
                vector->size = vector->size * 2;
                realloc(vector->locations, sizeof(vector->size));  
           }
         } 
     }
    if(sequence_end <= header_start){
        // 0 in sequence end will indicate sequence is incomplete
        vector->locations[vector->items] = {header_start, sequence_start, 0}; 
    }
     
    return vector;

}

*/

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

char get_data_type(const char* filepath){

    int err = 0;
    char first_char = '\0';
    char eof = -1;


    gzFile file = gzopen(filepath, "r");
    if(!file){ 
        fprintf(stderr, "Error: %s.\n", gzerror(file, &err));
        exit(EXIT_FAILURE);
    }

    first_char = gzgetc(file);
    if(first_char == eof){ 
        fprintf(stderr, "File is empty: %s.\n", filepath);
        exit(EXIT_FAILURE);
    }
    
    gzclose(file);
    if(first_char == '>'){
        return FASTA;
    }else if(first_char == '@'){
        return FASTQ;
    }else{ 
        fprintf(stderr, "Cannot discern file type for: %s.\n", filepath);
        exit(EXIT_FAILURE);
    } 

}


/* 
 * Extract fasta sequences from the input buffer, 
 * and create sequences, filling a storage buffer
 * of the input sequence data. 
 *
 * Over view of this algo is that we will ignore comments in the fastas
 * as those are rare. We will get the start sequence '>' then walk forward
 * till we find a '\n'. Upon find a new line charactar we will copy in charactars
 * discarding newline charactars until we reach another '>' or '\0'.
 *
 * risk getting only partial fastas doing this however from the buffer
 * */
void extract_fasta(const _sequence_file_data* fasta_data, const char* sequence_buffer){
   
   /*   
    const char* seq_data = fasta_data->data; 
    bool seq_data_start = false;
    for(ssize_t i = 0; i < fasta_data->length; i++ ){
        if(seq_data[i] == '>') {
            while(seq_data[i] != '>' && i < fasta_data->length){
                i++;
                if(seq_data[i] == '\n'){
                    seq_data_start = true; 
                    break
                }  
            
            }  
        }else if(seq_data_start){ 
            while(seq_data[i] != '>' && i < fasta_data->length){
        
        }
    
    
    }
    */

}

/*  

bool read_file(const char* filepath){
    
    //bool is_gzipped = is_gzip(filepath);
    const size_t buffer_length = 4096;
    char* buffer = (char*)malloc(sizeof(char) * buffer_length); 
    char* temporary_buffer = (char*)malloc(sizeof(char) * buffer_length)
    _sequence_file_data* sequence_data = init_file_data(buffer_length);
    
    gzFile file = gzopen(filepath, "r");
    if(!file){
        fprintf(stderr, "Could not open file: %s. failed: %s \n", 
                filepath, strerror(errno)); 
        exit(EXIT_FAILURE);
    }

    char data_type = get_data_type(filepath);
    
    bool sequence_start = false;
    if(data_type == FASTA){
        bool sequence_started = false;
        while(1){
            int err = 0;
            int bytes_read = gzread(file, buffer, buffer_length-1); 
            buffer[bytes_read] = '\0';
            // The CaratVector provides the locations of each sequence start and
            // stop in the buffer and if the sequence is located in the next buffer
            CaratVector* buffer_seq_data = find_fasta_sequences(buffer);
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

            for(size_t i = 0; i < buffer_seq_data->items; i++){
                SequenceIndex seq_data = buffer_seq_data->locations[i]; 
                if(seq_data.sequence_end == 0 && seq_data.header_start == 0){
                    // no sequence read just skip
                    break; 
                }
                if(seq_data.sequence_start < seq_data.sequence_end){
                    // Store sequence in a temporary array 
                    // till we get the rest
                    
                }else{
                    size_t data_length = seq_data.sequence_end - seq_data.sequence_start;
                    add_data(sequence_data, &buffer[seq_data.sequence_start], data_length);  
                }
            }
            
            // Started to find sequence, keep copying it into a buffer
            //add_data(sequence_data, sequence_next) 
        }    


        bool sequence_found = false;
        size_t sequence_start_pos = 0;
        while(1){
            int err;
            int bytes_read = gzread(file, buffer, buffer_length - 1);

            buffer[bytes_read] = '\0';      
            // Get bytes and look for start token (>)
            // walk ahead until next start toke (>) or EOF is found
            // performing additional reads of the file
            // add seq to stack
            if(buffer[0] == '>') {
            // Check if more than one seqeunce is present in the buffer
            char* next_occurence = strchr(buffer[1], '>');    
            // this needs to be rethought...
            if(next_occurence){}
            
            }
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
    }
    gzclose(file);



    return true;
}


*/
