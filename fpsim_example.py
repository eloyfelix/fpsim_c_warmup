from fpsim import load_fps, run_similarity
import time

sim_thres = 0.7

# aspirin
query_smiles = "CC(=O)Oc1ccccc1C(=O)O"

t0 = time.time()
moldb = load_fps("chembl.bin")
t1 = time.time()
print("Time loading the fp file", t1-t0)

print("")

t0 = t1
res = run_similarity(query_smiles, sim_thres, moldb)
t1 = time.time()
print("Similarity search against chembl", t1-t0)
print("Num results:", res.size)
print("Results:")
for i in range(res.size):
    print(res.simres[i].mol_id, res.simres[i].sim)

print("")

t0 = t1
for i in range(85):
    res = run_similarity(query_smiles, sim_thres, moldb)
t1 = time.time()
print("Similarity search against chembl * 85 (~UCHEM size)", t1-t0)

