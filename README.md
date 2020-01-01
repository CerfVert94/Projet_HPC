# Projet_HPC

Roqyun Ko:
Progres
1. L'optimisation normale / scalaire d'une operation de morphologie est finie.

#TODO 
1. L'optimisation normale / scalaire des operations sequentielles de morpho.
2. Le traitement du bord pour minimiser les erreurs qui se produisent. 
  a.  Les valeurs de bord ne sont ni initialisees lors de la creation ni ecrasees lors d'un operation de morpho.
  b. Les valeurs aleatoires sont le facteur principal d'un mauvais resultat lors d'un traitement enchaine des morpho (Erosion - Dilatation - Dilatation - Erosion). 
3. Restructurer les entetes des fonction de morpho pour 
  - potentiellement enlever le parametre "struct struct_element_dim s". Il est peu utile. 
  - ajouter un parametre pour un buffer intermediaire. 
    a. la fonction avec la meilleure optimisation perd le temps avec creation d'un buffer intermediaire.
    b. le traitement enchaine (E-D-D-E) en a besoin. 
4. L'optimisation normale / scalaire de Sigma Delta
5. Trouver un moyen pour optimiser le memoire. Il y a trop de defauts de cache.
