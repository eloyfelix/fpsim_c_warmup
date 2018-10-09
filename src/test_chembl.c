#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "sim_search.h"

int main(int argc, char *argv[]) {
    // 0 program name, 1 threshold, 2 input file, 3 fp query string, 4 fp size
    if(argc < 5) {
        fprintf(stderr, "\nYou forgot the parameters!\n"); 
        exit(1);
    }
    
    double threshold = atof(argv[1]);
    int fp_size = atoi(argv[4]);
    char query_string[1024];
    strcpy(query_string, argv[3]);
    char input_filename[256];
    strcpy(input_filename, argv[2]);

    // takes a while :( needs a binary format...
    MolDB *moldb = load_stuff(input_filename);

    // start the clock to time stuff
    clock_t begin = clock();

    // run the similarity search
    Result *results = similarity_search(query_string, threshold, moldb, fp_size);

    // time how long took the file load into memory
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("%fs computing similarity\n", time_spent);

    // write results to a file and free memory
    write_results_to_file(results);
    free(results->simres);
    free(results);
    return 0;
}
