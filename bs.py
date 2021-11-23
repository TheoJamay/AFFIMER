from Bio.PDB import Select, PDBIO, Selection
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.SASA import ShrakeRupley


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
    def __init__(self, ab,list_asa):
        self.chain = ab
        self.atom = list_asa

    def accept_atom(self, atom):
        if atom.get_full_id()[2] == self.chain :
            return 1
        else : 
            if atom.get_serial_number() in self.atom :
                return 1
            else:
                return 0


def decoup(structure,res_list,ab) :
    i=2
    for chain in structure.get_chains() :
        if chain.get_full_id()[2] != "L" and chain.get_full_id()[2] != "H" :
            res_tot = res_list.copy()
            for residue in chain :
                res_tot.append(residue)
            pdb_chain_file = "spilt_tempo{}_{}.pdb".format(ab,i)
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
            pdb_chain_file = "{}{}:{}.pdb".format(file[:-4],ab,chain.get_full_id()[2])
            io_w_no_h = PDBIO()
            io_w_no_h.set_structure(structure)
            io_w_no_h.save("bank_bs/{}".format(pdb_chain_file), AtomSelect(ab, list_asa))


if __name__ == '__main__':
    with open("../test_bdd/bdd_flavivirus.txt","r") as bdd_file :
        for line in bdd_file :
            words = line.split()
            file = words[0].strip()
            print(file)
            virus = words[1].strip()
            parser = PDBParser()
            structure = parser.get_structure(file[:-4], "../test_bdd/"+file)
            res_L = Selection.unfold_entities(structure[0]["L"], "R")
            res_H = Selection.unfold_entities(structure[0]["H"], "R")
            chains = list(structure.get_chains())
            decoup(structure,res_L,"L")
            decoup(structure,res_H,"H")
            if len(chains) == 3 :
                struc_H = parser.get_structure("HPDB","spilt_tempoH_2.pdb")
                seq_bind_site(struc_H,file,"H")
                struc_L = parser.get_structure("LPDB","spilt_tempoL_2.pdb")
                seq_bind_site(struc_L,file,"L")
            else :
                for i in range(2,len(chains)) :
                    struc_H = parser.get_structure("HPDB","spilt_tempoH_{}.pdb".format(i))
                    seq_bind_site(struc_H,file,"H")
                    struc_L = parser.get_structure("LPDB","spilt_tempoL_{}.pdb".format(i))
                    seq_bind_site(struc_L,file,"L")