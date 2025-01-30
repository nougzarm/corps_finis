# Définition, rappels 

  Un corps fini, souvent noté F_q (où q désigne son nombre d'éléments) est un anneau dont tous les éléments non nuls sont inversibles pour la multiplication. Tout corps à q élément peut être noté F_q sans ambiguité puisqu'il est unique à isomorphisme près. Etant donné un corps fini F_q, on peut montrer que q est forcément une puissance d'un nombre premier p. Il existe donc un entier naturel e tel que : q = p^e. Réciproquement, pour toute telle paire (p,e), il existe un corps fini de cardinal p^e.

Une façon de construire des corps finis :
  Etant donné un polynome f de degré e et à coefficients dans Z/pZ, si f est irréducible alors l'anneau quotient des polynômes (Z/pZ[X])/(f) est un corps car l'idéal (f) engendré par f est maximal. De plus, ce corps (qui peut aussi être vu comme un espace vectoriel) est de cardinal p^e. C'est pourquoi, les éléments de corps finis peuvent être vus comme des (images de) polynômes.

  La librairie corps_finis contient les outils de base pour travailler dans les espaces des polynômes Z[X] et Z/pZ[X] ou encore dans les corps finis F_q.

  Le programme compilé via le fichier test propose quelques exemples d'opérations dans ces espaces.
  