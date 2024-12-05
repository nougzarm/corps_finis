# corps-finis

  Un corps fini, souvent noté F_q (où q désigne son nombre d'éléments) est un anneau dont tous les éléments, sauf 0 sont inversibles pour la multiplication. N'importe quel corps de cardinal q peut être désigné par F_q sans ambiguité puisqu'un tel corps est unique à isomorphisme près. Il se trouve que q est forcément une puissance d'un nombre premier, on se donnera donc le premier p et le nombre naturel e tels que : q = p^e. À l'inverse, pour toute telle paire (p,e), il existe un corps fini de cardinal p^e.

Une façon de construire des corps finis :
  Etant donné un polynome f de degré e et à coefficients dans Z/pZ, si f est irréducible alors l'anneau quotient des polynômes (Z/pZ[X])/(f) est un corps car l'idéal (f) engendré par f est maximal. De plus, ce corps (qui peut aussi être vu comme un espace vectoriel) est de cardinal p^e. C'est pourquoi, les éléments de corps finis peuvent être vus comme des (images de) polynômes.

  La librairie corps-finis contient les outils de base pour travailler dans les espaces des polynômes Z[X] et Z/pZ[X] ou encore dans les corps finis F_q.

  Le programme compilé via le fichier test propose quelques exemples d'opérations dans ces espaces.
