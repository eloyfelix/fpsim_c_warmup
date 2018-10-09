from ctypes import POINTER, c_uint64, c_int, c_double, c_char, Structure, cdll, byref, c_char_p
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import textwrap
import time
import os

sslib = cdll.LoadLibrary(os.path.join(os.path.dirname(__file__), 'sim_search.so'))

# FP_SIZE = 16

class SimRes(Structure):
    _fields_ = [('sim', c_double),
                ('mol_id', c_int)]

class Result(Structure):
    _fields_ = [('size', c_int),
                ('simres', POINTER(SimRes))]

class Molecule(Structure):
    _fields_ = [('fp', POINTER(c_uint64)),
                ('mol_id', c_int)]

class MolDB(Structure):
    _fields_ = [('size', c_int),
                ('molecules', POINTER(Molecule))]

############################################################################# 

# def load_query_mol(c_mol):
#     load_query_mol = sslib.load_query_mol
#     load_query_mol.restype = POINTER(c_uint64)
#     c_query = load_query_mol(c_mol, c_int(FP_SIZE))
#     py_mol = ','.join([str(c_query[x]) for x in range(FP_SIZE)])
#     return c_query

#############################################################################


def load_fps(filename):
    filename = filename.encode()
    c_filename = ((c_char * len(filename))(*filename))
    load_stuff = sslib.load_stuff
    load_stuff.restype = POINTER(MolDB)
    moldb = load_stuff(c_filename)#, c_int(fp_size))
    return moldb

#############################################################################


def run_similarity(smiles, threshold, moldb):
    rdmol = Chem.MolFromSmiles(smiles)
    if rdmol:
        # rdkit canonicalise
        smiles = Chem.MolToSmiles(rdmol, canonical=True, isomericSmiles=True)
        rdmol = Chem.MolFromSmiles(smiles)
        fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(rdmol, 2, nBits=1024, useChirality=True, useBondTypes=True, useFeatures=True)
        splited = textwrap.wrap(fp.ToBitString(), 64)
        fp_size = len(splited)
        molfp = ','.join([str(int(x, 2)) for x in splited])

        molfp_b = molfp.encode()
        c_mol = ((c_char * len(molfp_b))(*molfp_b))

        similarity_search = sslib.similarity_search
        similarity_search.restype = POINTER(Result)
        results = similarity_search(c_mol, c_double(threshold), moldb, c_int(fp_size))
        return results.contents
    else:
        return None
