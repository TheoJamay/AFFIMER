from Bio.PDB import *
from Bio.PDB.SASA import ShrakeRupley
import os


class ChainSelect(Select):
    """
    class permettant de ne conserver que nos chaîne d'intéret dans le nouveau target pdb créer
    Parameters
    ----------
    list_chain : liste contenant nos 2 chaîne d'intéret
    return
    -------
   	boolean
    """
    def __init__(self, list_chain):
        self.chain = list_chain

    def accept_chain(self, chain):
        if chain.get_id() in self.chain :
            return 1
        else:          
            return 0


def consecutive(list_pos,list_aa) :
    """
    Générateur qui retourne les séquences consécutive (+/- 5) de l'anticorps (H ou L) d'une taille sup à 4 acides aminés en contact avec l'antigène
    Parameters
    ----------
    list_pos : liste de la position des aa de l'anticorps en contact avec l'antigène
    list_aa : liste code 1 lettre des aa de l'anticorps en contact avec l'antigène
    Yield
    -------
   	seq : séquence consécutive (+/- 5) d'une taille sup à 4 acides aminés
    """
    pred = -10
    lgt = 1
    seq = ""
    for i in range(len(list_pos)) :
        now = list_pos[i]
        if now < pred+5 :
            lgt += 1
        if now > pred+5 or i == len(list_pos)-1 :
            if lgt >= 4 :
                if i == len(list_pos)-1 and now < pred+5 :
                    seq += list_aa[i]
                    yield seq
                else :
                    yield seq
            lgt = 1
            seq = ""
        seq += list_aa[i]
        pred = now


def seq_bind_site(structure,antigene,dict_aa) :
    """
    Calcul de l'accessibilité au solvant (SAS)du complexe puis de chaque chaîne séparée. 
    La zone de contact comprend les atomes dont l'accessibilité a changé.
    Parameters
    ----------
    structure : structure PDB parser
    antigene : nom chaîne antigène
    dict_aa : dictionnaire servant à la conversion du code 3 lettre aa en code 1 lettre
    Returns
    -------
   	list_pos : liste de la position des aa de l'anticorps en contact avec l'antigène
    list_aa : liste code 1 lettre des aa de l'anticorps en contact avec l'antigène
    """
    sr = ShrakeRupley()
    list_pos = []
    list_aa = []
    for chain in structure.get_chains() :
        if chain.get_full_id()[2] != antigene :
            list_ref = []
            sr.compute(structure, level="A") #SAS du complexe au niveau atomique
            for residue in chain :
                if residue.get_full_id()[3][0][0] != "H" : #On exclu les HETATM
                    for atom in residue :
                        list_ref.append(round(atom.sasa,2)) #list de référence contenant la SAS pour chaque atome du complexe
            sr.compute(chain, level="A") #SAS de la chaîne au niveau atomique
            i=0
            res=0        
            for residue in chain :
                if residue.get_full_id()[3][0][0] != "H" :
                    res += 1
                    flag = True
                    for a in residue :
                        if list_ref[i] != round(a.sasa,2) and flag :
                            list_pos.append(res)
                            list_aa.append(dict_aa[residue.get_resname()])
                            flag = False #Un atome du résidu en intéraction suffit pour considérer le résidu comme en contact avec l'antigène                     
                        i+=1
    return list_pos, list_aa


if __name__ == '__main__':
    dict_aa = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    """
    Le fichier bs_13_NS1.tx contient les binding sites identifiés par PatchSearch comme similaire à la surface NS1 : 5GS6
    Le fichier best_boucle_NS1.txt contient les séquences consécutive (+/- 5) d'une taille sup à 4 acides aminés identifier dans 
    le complexe d'origine des binding sites
    """
    with open("NS1_5GS6/dscore_sup_13.txt","r") as bs_file,  open("NS1_5GS6/best_boucle.txt","a") as out_file :
        for line in bs_file : 
            file = "{}.pdb".format(line[:6])
            antigene = line[6] #chaîne antigène à sélectionner
            anticorps = line[8] #chaîne anticorps (lourde ou lègère) à sélectionner
            parser = PDBParser()
            structure = parser.get_structure(file[:-4], "pdb_files/"+file) #structure complexe entier
            list_chain = [antigene,anticorps] # On souhaite conserver uniquement ces 2 chaîne d'intéret sur tout le complexe
            pdb_chain_file = "target_{}.pdb".format(file[:-4])
            io_w_no_h = PDBIO()
            io_w_no_h.set_structure(structure)
            io_w_no_h.save("{}".format(pdb_chain_file), ChainSelect(list_chain)) #Création d'un fichier pdb avec uniquement nos 2 chaînes d'intéret
            target_structure = parser.get_structure("TPDB","target_{}.pdb".format(file[:-4])) #structure du nouveau pdb
            os.remove("target_{}.pdb".format(file[:-4]))
            list_pos, list_aa = seq_bind_site(target_structure,antigene,dict_aa)
            out_file.write(line.strip()) 
            for seq in consecutive(list_pos,list_aa) :
                out_file.write(" {}".format(seq))
            out_file.write("\n")