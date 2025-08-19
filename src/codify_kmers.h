/*
 * =====================================================================================
 *
 *       Filename:  codify_kmers.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08/18/25 10:57:21
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */

#ifndef _CODIFY_KMERS_H_
#define _CODIFY_KMERS_H_

#include <stdlib.h>
#include "sequence.h"

static const unsigned short int MAX_KMER_LENGTH = 256;
static const unsigned char chunk_size = 4;
// max_kmer_length / chunk_size
#define code_buffer_size  64


// Buffer struct for storage of codes
typedef struct CodeArena{
    unsigned short int code_size;
    size_t items;
    size_t size;
    unsigned char* codes;
}CodeArena;

CodeArena* init_code_arena(unsigned short int code_size, size_t starting_size);

void add_code(unsigned char* buffer, CodeArena* code_arena);

void destroy_code_arena(CodeArena* code_arena);

void compress_kmers(const SeqData*, size_t kmer_length);

void print_codes(CodeArena* code_arena);

#endif
