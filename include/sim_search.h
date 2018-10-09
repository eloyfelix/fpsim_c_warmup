#include <stdint.h>
#include <stdio.h>

// structure to store fingerprints with id
typedef struct Molecule {
    uint64_t *fp;
    int   mol_id;
} Molecule;  

// structure to store mol db
typedef struct MolDB {
    int size;
    Molecule  *molecules;
} MolDB;  

// structure to store a similarity result
typedef struct SimRes {
    double sim;
    int   mol_id;
} SimRes;

// structure to store all results
typedef struct Result {
    int size;
    SimRes  *simres;
} Result;  

const char* getfield(char* line, int num);

int comparator(const void *p1, const void *p2);

void *aligned_malloc(size_t alignment, size_t size);

double similarity(int size, uint64_t fp1[], uint64_t fp2[]);

int count_mols_in_fp_file(char filename[]);

MolDB * load_fp_struct(int n_mols, char filename[], int fp_size);

uint64_t * load_query_mol(char query_string[], int fp_size);

Result * similarity_search(char query_string[], double threshold, MolDB * moldb, int fp_size);

void write_results_to_file(Result* results);

MolDB * load_stuff(char input_filename[]);
