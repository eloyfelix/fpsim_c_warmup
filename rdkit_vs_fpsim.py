from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit import DataStructs
import csv


query_smiles = "CC(=O)Oc1ccccc1C(=O)O"

############################################ RDKIT SIMILARITY


def gen_fp(smiles):
    rdmol = Chem.MolFromSmiles(smiles)
    if rdmol:
        # rdkit "canonicalise"
        smiles = Chem.MolToSmiles(rdmol, canonical=True, isomericSmiles=True)
        rdmol = Chem.MolFromSmiles(smiles)
        if rdmol:
            fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(rdmol, 2, nBits=1024, useChirality=True, useBondTypes=True, useFeatures=True)
            return fp
    return None


def calc_similarity(query_fp, fps):
    sims = map(lambda x: [x[0], DataStructs.FingerprintSimilarity(query_fp, x[1])], fps)
    return list(sims)

with open('100mols.smi') as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')
    mols = []
    for row in readCSV:
        mols.append([int(row[0]), row[1]])

query_fp = gen_fp(query_smiles)
fps = list(filter(lambda x: x is not None, map(lambda x: [x[0], gen_fp(x[1])], mols)))
rdsims = calc_similarity(query_fp, fps)
# sort by id, easier comparison
rdsims.sort(key = lambda x:x[0], reverse=True)

######################################## FPSIM SIMILARITY


from fpsim import load_fps, run_similarity

# keep all sims
sim_thres = -0.01

# load the DB
moldb = load_fps("100mols.bin")

# same query_smiles
res = run_similarity(query_smiles, sim_thres, moldb)

fpsims = []
for i in range(res.size):
    fpsims.append([res.simres[i].mol_id, res.simres[i].sim])

# sort by id, easier comparison
fpsims.sort(key = lambda x:x[0], reverse=True)

######################################## COMPARE ALL SIMILARITIES

print("Same results for all mols?", fpsims == rdsims)
print("RDKit first 10 sims(ordered by mol_id):\n", rdsims[0:10])
print("FPSim first 10 sims(ordered by mol_id):\n", fpsims[0:10])
