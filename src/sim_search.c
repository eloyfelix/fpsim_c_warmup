#include <stdlib.h>
#include <string.h>
#include "score_functions.h"
#include "sim_search.h"


// function to get fields from csv line
const char* getfield(char* line, int num)
{
	const char* tok;
	for (tok = strtok(line, ",");
			tok && *tok;
			tok = strtok(NULL, ",\n"))
	{
		if (!--num)
			return tok;
	}
	return NULL;
}

// comparator function to run qsort
int comparator(const void *p1, const void *p2)  
{ 
    double l = ((struct SimRes *)p1)->sim; 
    double r = ((struct SimRes *)p2)->sim;
    if (l > r)
        return -1;
    else if (r > l)
        return 1;
    else
        return 0;
} 

// function to allocate memory needed for SIMD instructions
void *aligned_malloc(size_t alignment, size_t size) {
    void *mem;
    if (posix_memalign(&mem, alignment, size)) exit(1);
    return mem;
}

// similarity function
double similarity(int size, uint64_t fp1[], uint64_t fp2[]) {
    double t_score;
    tanimoto_score(fp1, fp2, size, &t_score);
    return t_score;
}

// counts mols in fp file
int count_mols_in_fp_file(char filename[]){
    FILE* input_fp;
    if ((input_fp = fopen(filename, "r")) == NULL)
    {
        fprintf(stderr, "\nError opening file\n"); 
        exit(1);   
    }
    char line[1024];
	int n_mols=0;
	// counts the number of mols in the file	
	while(fgets(line, 1024, input_fp))
	{
        n_mols++;
	}
	//close the file
	fclose(input_fp);	
    return n_mols;
}

// counts mols in fp bin file
int count_mols_in_fp_bin_file(char filename[]){
    FILE* input_fp;
    if ((input_fp = fopen(filename, "rb")) == NULL)
    {
        fprintf(stderr, "\nError opening file\n"); 
        exit(1);   
    }

    uint64_t fp_input;
    int mol_id_input;
    int fp_size;
    int n_mols = 0;

    // read fp_size
    fread(&fp_size, sizeof(int), 1, input_fp);

    int remaining_mols = 1;
    int fp_size_loop = fp_size;
    while(remaining_mols) {
        remaining_mols = fread(&mol_id_input, sizeof(int), 1, input_fp);
        if (remaining_mols)
            n_mols += 1;
        while (remaining_mols && (fp_size_loop > 0)){
            remaining_mols = fread(&fp_input, sizeof(uint64_t), 1, input_fp);
            fp_size_loop -= 1;
        }
        fp_size_loop = fp_size;
    }
	//close the file
	fclose(input_fp);	
    return n_mols;
}


// open the file to load the struct of fps from csv
MolDB * load_fp_struct(int n_mols, char filename[], int fp_size){
    FILE* input_fp;

    if ((input_fp = fopen(filename, "r")) == NULL)
    {
        fprintf(stderr, "\nError opening file\n"); 
        exit(1);   
    }

    MolDB *moldb = (MolDB *) malloc(sizeof(MolDB));
    moldb->molecules = malloc(sizeof(Molecule) * n_mols);
    moldb->size = n_mols;

    char line[1024];
    // fill the array of structs with the contents of the CSV file
    for (int i=0; i < n_mols; ++i ) {
        fgets(line, 1024, input_fp);
        // fp_size + 1 (mol_id + fp_size)
        for(int k = 0; k < fp_size + 1; ++k) {
            char* tmp = strdup(line);
            if (k == 0){
                moldb->molecules[i].mol_id = atoi(getfield(tmp, k+1));
                moldb->molecules[i].fp = aligned_malloc(32, sizeof(uint64_t) * fp_size);
            }
            else{
                moldb->molecules[i].fp[k-1] = strtoul(getfield(tmp, k+1), NULL, 10);
            }
            free(tmp);
        }
    }
    fclose(input_fp);
    return moldb;
}

// open the file to load the struct of fps from binary file
MolDB * load_fp_bin_struct(int n_mols, char filename[]) 
{ 
    FILE *infile = fopen(filename, "rb"); 
    if (infile == NULL) 
    { 
        fprintf(stderr, "\nError opening file\n"); 
        exit (1); 
    } 
    int mol_id_input;
    int fp_size;

    // read fp_size
    fread(&fp_size, sizeof(int), 1, infile);
    uint64_t fp_input;

    MolDB *moldb = (MolDB *) malloc(sizeof(MolDB));
    moldb->molecules = malloc(sizeof(Molecule) * n_mols);
    moldb->size = n_mols;

    for (int n = 0; n < n_mols; ++n){
        fread(&mol_id_input, sizeof(int), 1, infile);
        moldb->molecules[n].mol_id = mol_id_input;
        moldb->molecules[n].fp = aligned_malloc(32, sizeof(uint64_t) * fp_size);
        for (int f = 0; f < fp_size; ++f){
            fread(&fp_input, sizeof(uint64_t), 1, infile);
            moldb->molecules[n].fp[f] = fp_input;
        }
    }
    // close file 
    fclose (infile);
    return moldb; 
} 


uint64_t * load_query_mol(char query_string[], int fp_size){
    int i = 0;
    // strtok modifies original value
    char *tmp = strtok(strdup(query_string), ",");
    uint64_t *query = malloc(fp_size * sizeof(uint64_t));
    while (tmp != NULL)
    {
        query[i++] = strtoul(tmp, NULL, 10);
        tmp = strtok(NULL, ",");
    }
    free(tmp);
    return query;
}

// similarity search function, returns results ordered by similarity
Result * similarity_search(char query_string[], double threshold, MolDB * moldb, int fp_size){
    // load query mol
    uint64_t *query = load_query_mol(query_string, fp_size);

    int simres_length = 128;
    Result *results = (Result *) malloc(sizeof(Result));
    results->simres = (SimRes *) malloc(sizeof(SimRes) * simres_length);

    int total_sims = 0;
    double sim;
    for (int i=0; i < moldb->size; ++i) {
        sim = similarity(fp_size, query, moldb->molecules[i].fp);
        if (sim > threshold){
            results->simres[total_sims].mol_id = moldb->molecules[i].mol_id;
            results->simres[total_sims].sim = sim;
            total_sims++;
        }
        if (total_sims == simres_length){
            simres_length *= 2;
            // reallocating memory
            results->simres = (SimRes*)realloc(results->simres, simres_length * sizeof(SimRes));
        }
    }
    results->size = total_sims;
    free(query);

    // sort by similarity
    qsort(results->simres, total_sims, sizeof(struct SimRes), comparator);
    return results;
}

// write similarity results into a file
void write_results_to_file(Result* results){
    FILE* out_fp;
    if ( (out_fp = fopen("output.bin", "wb")) == NULL )
    {
        fprintf(stderr, "\nError opening file\n"); 
        exit(1);   
    }
    fwrite(results->simres, sizeof(struct SimRes) * results->size, 1, out_fp);
    fclose(out_fp);
}

// load fp file into memory
MolDB * load_stuff(char input_filename[]){
    // count the number of fps in the input file
    // int n_mols = count_mols_in_fp_file(input_filename);
    int n_mols = count_mols_in_fp_bin_file(input_filename);

    // load struct of fps with the file
    // MolDB *moldb = load_fp_struct(100, input_filename, fp_size);
    MolDB *moldb = load_fp_bin_struct(n_mols, input_filename);

    return moldb;
}
