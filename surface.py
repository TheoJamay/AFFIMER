from Bio.PDB import Select, PDBIO, Selection
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from os import listdir
from os.path import isfile, join
import os


class ChainSelect(Select):
    def __init__(self, list_chain):
        self.chain = list_chain

    def accept_chain(self, chain):
        if chain.get_id() in self.chain :
            return 1
        else:          
            return 0


class AtomSelect(Select):
    def __init__(self,a_tot):
        self.atom = a_tot

    def accept_atom(self, atom):
        if atom in self.atom :
            return 1
        else:
            return 0


def decoup(structure,file) :
    list_chain = []
    for chain in structure.get_chains() :
        if chain.get_full_id()[2] != "L" and chain.get_full_id()[2] != "H" :
            list_chain.append(chain.get_full_id()[2])
    pdb_chain_file = "spilt_{}.pdb".format(file[:-4])
    io_w_no_h = PDBIO()
    io_w_no_h.set_structure(structure)
    io_w_no_h.save("{}".format(pdb_chain_file), ChainSelect(list_chain))


if __name__ == '__main__':
    w=True
    files = [f for f in listdir("flavivirus_AbDb") if isfile(join("flavivirus_AbDb", f))]
    for file in files :
        print(file)
        parser = PDBParser()
        structure = parser.get_structure(file[:-4], "flavivirus_AbDb"+file)
        sr = ShrakeRupley()
        chains = list(structure.get_chains())
        a_tot = []
        if len(chains) == 3 :
            print(len(chains))
            for chain in structure.get_chains() :
                if chain.get_full_id()[2] != "L" and chain.get_full_id()[2] != "H" :
                    sr.compute(chain, level="A")
                    for residue in chain :
                        if residue.get_full_id()[3][0][0] != "H" :
                            for a in residue :
                                if round(a.sasa,2) != 0 :
                                    a_tot.append(a)
                    pdb_chain_file = "{}{}.pdb".format(file[:-4],chain.get_full_id()[2])
                    io_w_no_h = PDBIO()
                    io_w_no_h.set_structure(structure)
                    io_w_no_h.save("surface/{}".format(pdb_chain_file), AtomSelect(a_tot))
        else :
            print(len(chains))
            decoup(structure,file)
            struc_H = parser.get_structure("HPDB","spilt_{}.pdb".format(file[:-4]))
            os.remove("spilt_{}.pdb".format(file[:-4]))
            chains = list(structure.get_chains())
            sr.compute(struc_H, level="A")
            for chain in struc_H.get_chains() :
                for residue in chain :
                    if residue.get_full_id()[3][0][0] != "H" :
                        for a in residue :
                            if round(a.sasa,2) != 0 :
                                a_tot.append(a)
            for chain in struc_H.get_chains() :
                pdb_chain_file = "{}{}.pdb".format(file[:-4],chain.get_full_id()[2])
                io_w_no_h = PDBIO()
                io_w_no_h.set_structure(struc_H)
                io_w_no_h.save("surface/{}".format(pdb_chain_file), AtomSelect(a_tot))