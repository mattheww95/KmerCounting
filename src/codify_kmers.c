/*
 * =====================================================================================
 *
 *       Filename:  codify_kmers.c
 *
 *    Description: Generate compression codes for the ingested k-mers.  
 *
 *        Version:  1.0
 *        Created:  08/18/25 10:29:57
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */


#include "codify_kmers.h"


static const unsigned char A = 0b00000000; // 0
static const unsigned char T = 0b00000001; // 1, 4, 16, 64
static const unsigned char C = 0b00000010; // 2, 8, 32, 128
static const unsigned char G = 0b00000011; // 3, 12, 48, 192


unsigned char set_value(char input, bool* bad_kmer){
    const unsigned char error_value = -1; 
    switch(input){
    case 'A':
       return A;
    case 'a':
      return A; 
    case 'T':
       return T;
    case 't':
      return T; 
    case 'C':
       return C;
    case 'c':
      return C; 
    case 'G':
       return G;
    case 'g':
      return G; 
    default:
      //Setting a flag on return was not working
      //so a reference to a varible is set to avoid a goto
      // although a goto may be faster...
        *bad_kmer = true; 
        return error_value; 
    }

}


//Create the code arena for storing new codes
//The number of codes stored is the starting size time the code size
//as the buffer is treated as blocks.
CodeArena* init_code_arena(unsigned short int code_size, size_t starting_size){
    
    CodeArena* code_arena = (CodeArena*)malloc(sizeof(CodeArena));
    code_arena->items = 0;
    code_arena->code_size = code_size; 
    code_arena->size = code_size * starting_size;
    code_arena->codes = (unsigned char*)malloc(sizeof(unsigned char) * starting_size * code_size);

    return code_arena;
}


void destroy_code_arena(CodeArena* code_arena){
    free(code_arena->codes);
    free(code_arena);
}


void add_code(unsigned char* buffer, CodeArena* code_arena){
    
    if((code_arena->code_size * code_arena->items) >= code_arena->size){
        code_arena->codes = realloc(code_arena->codes, code_arena->size * 2 ); 
        code_arena->size = code_arena->size * 2; 
    }
    size_t start_pos = code_arena->items * code_arena->code_size;
    
    // add in code size here
    for(size_t i = 0; i < code_arena->code_size; i++){
         code_arena->codes[i + start_pos]  = buffer[i];
    }
    code_arena->items++;
}


void print_codes(CodeArena* code_arena){

    size_t items = code_arena->items * code_arena->code_size;
    for(size_t i = 0; i < items; i=i+code_arena->code_size){
        size_t y = 0;
        for(;y < code_arena->code_size-1; y++){
            printf("%d.", code_arena->codes[i+y]);         
        }
        printf("%d\n", code_arena->codes[i+y]);         
    }
}

void compress_kmers(const SeqData* seq_data, size_t kmer_length){
    
    // Allow k-mers 
    if(kmer_length > MAX_KMER_LENGTH){
        fprintf(stdout, "K-mer size selected exceeds max k-mer length of %d", MAX_KMER_LENGTH);
        exit(EXIT_FAILURE); 
    }
    
    unsigned char buffer[code_buffer_size] = {0};
    const unsigned short int error_value = 255;
    // integer ceiling division
    unsigned short int code_size = (kmer_length + chunk_size - 1) / chunk_size;
    const size_t starting_buffer_size = 1024; 

    CodeArena* code_arena = init_code_arena(code_size, starting_buffer_size);

    for(size_t i = 0; i <= seq_data->length-kmer_length; i++){
        // compression loop here 
        bool bad_kmer = false;
        size_t counter = 0;
        size_t y = 0;
        // Issue getting all k-mers for shorter values
        for(;y <= kmer_length - chunk_size; y=y+chunk_size){
            
            //unsigned short int value = 0;
            unsigned char value = 0;
            char v1 = seq_data->sequence[i+y], 
                 v2 = seq_data->sequence[i+y+1], 
                 v3 = seq_data->sequence[i+y+2], 
                 v4= seq_data->sequence[i+y+3];
           
            
            value = value | set_value(v1, &bad_kmer);
            value = value | (set_value(v2, &bad_kmer) << 2);
            value = value | (set_value(v3, &bad_kmer) << 4); 
            value = value | (set_value(v4, &bad_kmer) << 6); 
              
            if(value > error_value){
                bad_kmer = true;
                break;  
            }
            
            buffer[counter] = value;
            counter++;
        }

        // If non-zero value is present we have some parts of the k-mer left over
        size_t leftover = kmer_length - y;
        if(leftover > 0){
            unsigned char value = 0;
            for(size_t f = 0; f < leftover; f++){
                value = value | (set_value(seq_data->sequence[i+y+f], &bad_kmer) << (f * 2));
            }
  
            if(value > error_value){
                bad_kmer = true;
            }

            buffer[counter] = value;
            counter++;
        }

        if(bad_kmer){
            continue; 
        }
        add_code(buffer, code_arena);

#ifdef DEBUG
        for(size_t f = 0; f < counter; f++){
            printf("%d.", buffer[f]);
        }
        printf("\n");
        for(int t = 0; t < kmer_length; t++){ 
            printf("%c", seq_data->sequence[i + t]);
        }
        printf("\n");
#endif
        memset(buffer, 0, code_buffer_size);
    }

#ifdef DEBUG
    print_codes(code_arena);
#endif

}



