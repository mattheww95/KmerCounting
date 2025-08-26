/*
 * Create an arena of sequence data for either
 * fastq data or fasta data. This will consist 
 * of a struct of charactar arrays holding the
 * individual sequnce information e.g. array of
 * pointers.
 * */

#ifndef _SEQUENCE_H_
#define _SEQUENCE_H_

#include "zlib.h"
#include <stdio.h>
#include <stdbool.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "kseq.h"

static const unsigned char A = 0b00000000; // 0
static const unsigned char T = 0b00000001; // 1, 4, 16, 64
static const unsigned char C = 0b00000010; // 2, 8, 32, 128
static const unsigned char G = 0b00000011; // 3, 12, 48, 192
static const size_t starting_buffer_size = 1024; // 3, 12, 48, 192


typedef struct SeqData{
    size_t size;
    size_t length;
    char* sequence;
}SeqData;

void add_sequence_data(SeqData* data, char* new_data, size_t data_length);


static const unsigned short int MAX_KMER_LENGTH = 256;
static const unsigned char chunk_size = 4;
// max_kmer_length / chunk_size
#define code_buffer_size  64


// Buffer struct for storage of codes
typedef struct CodeArena{
    unsigned char code_size;
    size_t items;
    size_t size;
    // Maybe better to store data as a min-heap, pairing heap or fibonacci heap on insertion
    unsigned char* codes;
    // can probably write a function to convert each set of codes into a scalar value 
}CodeArena;

CodeArena* get_data(const char* filepath, size_t kmer_length);

CodeArena* init_code_arena(unsigned short int code_size, size_t starting_size);

void print_code_count(unsigned char* data, size_t code_size, size_t counter);

void add_code(unsigned char* buffer, CodeArena* code_arena);

void destroy_code_arena(CodeArena* code_arena);

void compress_kmers(const SeqData*, size_t kmer_length);

bool compare_codes(const unsigned char* c1, const unsigned char* c2, size_t code_length);

void print_codes(CodeArena* code_arena);

void compress_sequence(CodeArena* code_arena, char const* sequence, size_t length, size_t kmer_length);

#endif


