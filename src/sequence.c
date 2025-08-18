/*
 * =====================================================================================
 *
 *       Filename:  sequence.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08/18/25 13:23:54
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
#include "sequence.h"

KSEQ_INIT(gzFile, gzread);


void add_sequence_data(SeqData* data, char* new_data, size_t data_length){
    size_t padding = 1;
    if(data_length + data->length + padding >= data->size){
        size_t new_arr_size = (data->size * 2);
        data->sequence = (char*)realloc(data->sequence, sizeof(char) * new_arr_size);
        data->size = new_arr_size; 
    }
    memcpy(&data->sequence[data->length], new_data, data_length); 
    data->length = data->length + data_length;
    // Hash character denotes seperation between sequences, it will
    // prevent k-mers from different records that are side by side
    // in the buffer from showing up
    data->sequence[data->length] = '#';
    data->length += padding;
  
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

#ifdef DEBUG
    printf("In debug mode. \n");
    for(size_t i = 0; i <= seq_data->length; i++){
        printf("%c", seq_data->sequence[i]); 
    }
    printf("\n");
#endif
    return seq_data;
}

