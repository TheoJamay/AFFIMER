from Bio.PDB import Select, PDBIO, Selection
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from os import listdir
from os.path import isfile, join
from  joblib import Parallel, delayed
import os

class ResiduSelect(Select):
    def __init__(self, res_tot):
        self.residue = res_tot

    def accept_residue(self, residue):
        if residue in self.residue :
            if residue.get_full_id()[3][0][0] != "H" :
                return 1
            else :
                return 0
        else:
            return 0


class AtomSelect(Select):
    def __init__(self,list_asa):
        self.atom = list_asa

    def accept_atom(self, atom):
        if atom.get_serial_number() in self.atom :
            return 1
        else:
            return 0


def decoup(structure,res_list,ab,file) :
    i=2
    for chain in structure.get_chains() :
        if chain.get_full_id()[2] != "L" and chain.get_full_id()[2] != "H" :
            res_tot = res_list.copy()
            for residue in chain :
                res_tot.append(residue)
            pdb_chain_file = "spilt_tempo{}_{}_{}.pdb".format(ab,i,file[:-4])
            io_w_no_h = PDBIO()
            io_w_no_h.set_structure(structure)
            io_w_no_h.save("{}".format(pdb_chain_file), ResiduSelect(res_tot))
            i += 1


def seq_bind_site(structure,file,ab) :
    sr = ShrakeRupley()
    for chain in structure.get_chains() :
        if chain.get_full_id()[2] != "L" and chain.get_full_id()[2] != "H" :
            list_ref = []
            sr.compute(structure, level="A")
            for residue in chain :
                if residue.get_full_id()[3][0][0] != "H" :
                    for atom in residue :
                        list_ref.append(round(atom.sasa,2))
            sr.compute(chain, level="A")
            i=0
            list_asa = []
            for residue in chain :
                if residue.get_full_id()[3][0][0] != "H" :
                    for atom in residue :
                        if list_ref[i] != round(atom.sasa,2) :
                            list_asa.append(atom.get_serial_number())
                        i+=1
            if list_asa : #Il peut y avoir des cas où aucun atomes n'est en intéraction
                pdb_chain_file = "{}{}:{}.pdb".format(file[:-4],chain.get_full_id()[2],ab)
                io_w_no_h = PDBIO()
                io_w_no_h.set_structure(structure)
                io_w_no_h.save("binding_site_AbDb_Martin/{}".format(pdb_chain_file), AtomSelect(list_asa))


def work(file):
    print(file)
    parser = PDBParser()
    structure = parser.get_structure(file[:-4], "../AbDb/LH_Protein_Martin/"+file)
    res_L = Selection.unfold_entities(structure[0]["L"], "R")
    res_H = Selection.unfold_entities(structure[0]["H"], "R")
    chains = list(structure.get_chains())
    decoup(structure,res_L,"L",file)
    decoup(structure,res_H,"H",file)
    if len(chains) == 3 :
        struc_H = parser.get_structure("HPDB","spilt_tempoH_2_{}.pdb".format(file[:-4]))
        seq_bind_site(struc_H,file,"H")
        os.remove("spilt_tempoH_2_{}.pdb".format(file[:-4]))
        struc_L = parser.get_structure("LPDB","spilt_tempoL_2_{}.pdb".format(file[:-4]))
        seq_bind_site(struc_L,file,"L")
        os.remove("spilt_tempoL_2_{}.pdb".format(file[:-4]))
    else : #Si plus de 3 chaînes dans le pdb c'est probablement que l'antigène est multimérique
        for i in range(2,len(chains)) :
            struc_H = parser.get_structure("HPDB","spilt_tempoH_{}_{}.pdb".format(i,file[:-4]))
            seq_bind_site(struc_H,file,"H")
            os.remove("spilt_tempoH_{}_{}.pdb".format(i,file[:-4]))
            struc_L = parser.get_structure("LPDB","spilt_tempoL_{}_{}.pdb".format(i,file[:-4]))
            seq_bind_site(struc_L,file,"L")
            os.remove("spilt_tempoL_{}_{}.pdb".format(i,file[:-4]))
    

if __name__ == '__main__' :
    Parallel(n_jobs = 3, prefer = "processes")(delayed(work)(file) for file in listdir("../AbDb/LH_Protein_Martin") 
    if isfile(join("../AbDb/LH_Protein_Martin",file)))