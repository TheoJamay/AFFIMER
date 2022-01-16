# Projet Long Affimer

best_boucle.py retourne les séquences d’acides aminés consécutive (+/- 5) de l'anticorps (chaîne H ou L) d'une taille supérieure à 4 acides aminés en contact avec l'antigène pour les meilleurs binding sites.

surface.py convertit un complexe Ac/Ag en une surface protéique

binding_site_AbDb.py convertit les complexes Ac/Ag de l’AbDb en binding sites côté antigène.

Au sein du dossier codes_output_patchsearch on retrouve AC_select.py qui fait selectionne les binding sites les plus similaire. Ensuite sort_result.py ordonne les binding sites en  fonction du D-score.

Le dossier flavivirus contient des complexes Ac/Ag de flavivirus de l'AbDb.

Le dossier NS1_5GS6 contient best_boucle.txt qui le fichier de sortie de best_boucle.py. Le fichier 5GS6_score_sup10.txt contient les binding dites des chaîne A et B de 5GS6.pdb avec un score supérieur à 10. Le fichier dscore_sup_13.txt est utilisé par best_boucle.py.

Le dossier patchsearch_output contient les outputs généraux de différent antigène de flavivirus dont la NS1 sans anticorps 5GS6. 

Le dosser patchsearch_sort contient la version trié des outputs généraux contenu dans le dossier patchsearch_output.

Le dossier patchsearch_al contient les fichier output d'alignement fournit par patchsearch

Le dossier pdb_files est utile pour best_boucles.py
