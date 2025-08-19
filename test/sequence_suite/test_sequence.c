/*
 * =====================================================================================
 *
 *       Filename:  test_sequence.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08/12/25 14:56:02
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
#include "codify_kmers.h"
#include "unity.h"


void setUp(void){}
void tearDown(void){}



/*  
void test_isgzip(void){ 
    TEST_ASSERT_TRUE(is_gzip("test.gzip"));
}

  
void test_read_file(void){
    TEST_ASSERT_TRUE(read_file("test.txt"));
    TEST_ASSERT_TRUE(read_file("test.gzip"));
    TEST_ASSERT_TRUE(read_file("test_long.txt"));
}


  
void test_get_data_type(void){
    TEST_ASSERT(get_data_type("test.fasta") == FASTA);
    TEST_ASSERT(get_data_type("test.fastq") == FASTQ);
}
*/

void test_data_reader(void){
    TEST_ASSERT(get_data("realn02.fa"));
}


  
void test_compress_kmers(void){
    const SeqData* seq_data = get_data("realn02.fa");
    compress_kmers(seq_data, 13); 
}



int main(void){
    UNITY_BEGIN();
    RUN_TEST(test_data_reader);
    RUN_TEST(test_compress_kmers);
//    RUN_TEST(test_read_file);
//    RUN_TEST(test_get_data_type);
    return UNITY_END();
}

