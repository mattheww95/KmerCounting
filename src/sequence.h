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



typedef struct SeqData{
    size_t size;
    size_t length;
    char* sequence;
}SeqData;

void add_sequence_data(SeqData* data, char* new_data, size_t data_length);

SeqData* get_data(const char* filepath);


#endif


