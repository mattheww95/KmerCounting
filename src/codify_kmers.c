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


unsigned char set_value(char input){
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
      return error_value; 
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

    for(size_t i = 0; i <= seq_data->length-kmer_length; i++){
        // compression loop here 
        bool bad_kmer = false;
        size_t counter = 0;
        size_t y = 0;
        // Issue getting all k-mers for shorter values
        for(;y <= kmer_length - chunk_size; y=y+chunk_size){
            
            //unsigned short int value = 0;
            unsigned short int value = 0;
            char v1 = seq_data->sequence[i+y], 
                 v2 = seq_data->sequence[i+y+1], 
                 v3 = seq_data->sequence[i+y+2], 
                 v4= seq_data->sequence[i+y+3];

            value = value | set_value(v1);
            value = value | (set_value(v2) << 2);
            value = value | (set_value(v3) << 4); 
            value = value | (set_value(v4) << 6); 

            if(value > error_value){
                bad_kmer = true;
                break;  
            }
            buffer[counter] = value;
            counter++;
        }

        // If non-zero value is present we have some parts of the k-mer left over
        size_t leftover = kmer_length - y;
        if(leftover){
            unsigned short int value = 0;
            for(size_t f = 0; f < leftover; f++){
                value = value | (set_value(seq_data->sequence[i+y+f]) << (f * 2));
            }

            if(value > error_value){
                bad_kmer = true;
                continue;  
            }
            buffer[counter] = value;
            counter++;
        }

        if(bad_kmer){
            continue; 
        }

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

}



