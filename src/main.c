#include "sequence.h"
#include <getopt.h>


typedef struct Opt{
    const char *name;
    int has_arg;
    int* flag;
    int val;
    const char* help;
}Opt;


void print_help(Opt* opts){
    printf("Welcome to our slow k-mer counting program!\n\n");
    printf("Options:\n");
    size_t counter = 0;
    while(1){
        if(opts[counter].name == 0 && opts[counter].has_arg == 0){
            exit(EXIT_SUCCESS); 

        }
        Opt option = opts[counter];
        if(opts[counter].has_arg == required_argument){
            printf("  [REQUIRED] -%c|--%s:\n\t%s\n", 
                    option.val, 
                    option.name, 
                    option.help );
        }else{
            printf("  -%c|--%s:\n\t%s\n", 
                    option.val, 
                    option.name, 
                    option.help );
        }
        counter++; 
    }

}


int main(int argc, char** argv){
    int c = 0;
    static Opt long_options[] = {
        {"input",required_argument,0,'i', "Input file."}, 
        {"size",required_argument,0,'s',"Kmer size."}, 
        {"help",no_argument,0,'h', "Show help and exit."},
        {0,0,0,0,0} 
    };
    const char* input_file = NULL;
    size_t kmer_size = 0;
    
    if(argc <= 1){
        print_help(long_options); 
    }
    while(1){
    
        int option_index = 0; 
        c = getopt_long(argc, argv, "hi:s:", long_options, &option_index);
        if(c==-1)
            break;
        switch(c){
            case 'i':
                input_file = optarg;
                break;
             case 's':
                kmer_size = strtoull(optarg, NULL, 10);
                break;
             case 'h':
                print_help(long_options);
                break;
             default:
                print_help(long_options);
                exit(EXIT_FAILURE);
            
        }
    }


    get_data(input_file, kmer_size);
    return 0;
}
