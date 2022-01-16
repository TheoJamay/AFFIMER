from os import listdir
from os.path import isfile, join


def sort_bs(file) :
    with open("ps_out/"+file,"r") as bs_file :
        dic_bs = {}
        first = True
        for line in bs_file :
            if first :
                first = False
                continue
            words = line.split()
            pdb = words[0].split("/")
            if pdb[4] not in dic_bs :
                dic_bs[pdb[4]] = float(words[8])
            else :
                print("error")
            
    dic_sort = dict(sorted(dic_bs.items(), key=lambda item: item[1], reverse=True))

    with open("ps_sort/"+file[:-3]+"sort","w") as sort_file :
        for key in dic_sort :
            sort_file.write("{}\t{}\n".format(dic_sort[key],key))

# files = [f for f in listdir("ps_out") if f[-3:] == "out"]
# for file in files :
#     sort_bs(file)
sort_bs("5GS6_1A_out")
sort_bs("5GS6_1B_out")