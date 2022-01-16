from os import listdir
from os.path import isfile, join


def select_pdb(file,dic_selec) :
    list_doublon = []
    with open("../new_app/ps_sort/"+file,"r") as sort_file :
        for line in sort_file :
            words = line.split()
            score = float(words[0])
            if score > 10 :
                pdb = words[1][:4]
                ch = "{}_{}".format(pdb,words[1][-6])
                #print(ch)
                AC = words[1]
                if pdb != file[:-8] and ch not in list_doublon :
                    if not AC in dic_selec :
                        dic_selec[AC] = [(file[:-5],score)]
                    else :
                        dic_selec[AC].append((file[:-5],score))
                    list_doublon.append(ch)
            else :
                break


# with open("ps_sort_files.txt","r") as tot_file :
#     dic_selec = {}
#     for file in tot_file :
#         #print(file)
dic_selec = {}

select_pdb("5GS6_1A_sort",dic_selec)
select_pdb("5GS6_1B_sort",dic_selec)
# list_select = list(set(list_select))
# print(list_select)
with open("NS1/5GS6_10.txt","w") as out_file :
    for key in dic_selec :
        out_file.write("{} {}\n".format(key,dic_selec[key]))