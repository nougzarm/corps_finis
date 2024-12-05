#ifndef CORPSFINI_INCLUDED
#define CORPSFINI_INCLUDED

#include <stdio.h>
#include <stdlib.h>

// Structure ---------------------------------------------------
typedef struct{
    int* coeff;
    int degre;
} polynome;


/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                                     1. MATHS                                                   |   
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|                                                     
 */
int max(int a, int b);
int min(int a, int b);
int modulo(int a, int b);
int puissance(int a, int e);
int puissance_modulo(int m, int e, int p);
int inverse_mod(int a, int p);


/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                                  2. OUTILS                                                     |   
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|                                                     
 */
int conversion_scalaire(polynome* P);                       // Si P est constant alors renvoie son coeff
void scalaire(int n, polynome* P);                          // P -> nP
void scalaire_mod(int n, polynome* P, int p);               // P -> nP  (dans F_p[X])
void scalaire_Fq(int n, polynome* P, int p, polynome* f);   // P -> nP  (dans F_q)
void modulo_transfo(polynome* P, int p);                    // P -> P [p]
void unitaire(polynome* P, int p);                          // P [p] -> P/CoeffDom(P) [p]
polynome difference_etendu(polynome* A, polynome* Q, polynome* B);              // retourne A-QB   (dans Z[X])
polynome difference_etendu_mod(polynome* A, polynome* Q, polynome* B, int p);   // retourne A-QB   (dans F_p[X])
polynome surjection(polynome*, int p, polynome* f);                // Surjection Z[X] ->> F_p[X] ->> F_p[X]/(f)   
polynome demi_surjection(polynome*, int p, polynome* f);           // Surjection F_p[X] ->> F_p[X]/(f)


/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                      3. INITIALISATION DE POLYNOMES                                            |   
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|                                                     
 */
polynome polynull();                                    // Renvoie le polynome nul
polynome copie(polynome* P);                            // Renvoie une copie de P
polynome monome(int coefficient, int exposant);         // Renvoie le monome de degré et coeff voulu
polynome polynome_init(int* coefficient, int degre);    // Initialise un polynôme [A programmer]


/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                                   4. GESTION                                                   |   
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|                                                     
 */
void afficher(polynome*);     // Affiche le polynome
void vider(polynome*);        // Libère la mémoire allouée pour la déf du polynome


/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                           5. OPERATIONS DANS Z[X]                                              |   
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|                                                     
 */
polynome addition(polynome*, polynome*);                 // Addition dans Z[X]
polynome oppose(polynome*);                              // Opposé dans Z[X]
polynome soustraction(polynome*, polynome*);             // Soustraction dans Z[X]
polynome multiplication(polynome*, polynome*);           // Multiplication dans Z[X]
polynome puissance_polynome(polynome*, int exposant);    // Puissance dans Z[X]
polynome reduction_modulo(polynome*, int p);             // Réduit un polynome de Z[X] modulo p


/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                          6. OPERATIONS DANS F_p[X]                                             |   
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|                                                     
 */
polynome addition_mod(polynome*, polynome*, int p);                    // Addition dans F_p[X]
polynome soustraction_mod(polynome*, polynome*, int p);                // Soustraction dans F_p[X]
polynome multiplication_mod(polynome*, polynome*, int p);              // Multiplication dans F_p[X]
polynome puissance_mod(polynome*, int exposant, int p);                // Puissance dans F_p[X]
polynome division_euclid(polynome*, polynome*, int p, int i);          // Div. euclidienne dans F_p[X] (i=0: le quotient, i=1: le reste)
polynome algo_euclide(polynome*, polynome*, int p);                    // PGCD dans F_p[X] (Algorithme d'Euclide)
polynome algo_euclide_etendu(polynome* P, polynome* Q, int p, int i);  // Euclide étendu (si uP+vQ=pgcd alors i=0: v, sinon: u)


/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                        7. OPERATIONS DANS F_q = F_p[X]/(f)                                     |   
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|                                                     
 */
int cardinal(int p, polynome* f);                                      // Renvoie q = card(F_q)
polynome addition_Fq(polynome*, polynome*, int p, polynome* f);        // Addition dans F_q
polynome multiplication_Fq(polynome*, polynome*, int p, polynome* f);  // Multiplication dans F_q
polynome puissance_Fq(polynome*, int exposant, int p, polynome* f);    // Puissance dans F_q
polynome inverse(polynome*, int p, polynome* f);                       // Renvoie l'inverse d'un elt dans F_q
polynome division(polynome*, polynome*, int p, polynome* f);           // Division dans F_q
int ordre(polynome* P, int p, polynome* f);                            // Ordre de P dans F_q
int verif_generateur(polynome* P, int p, polynome* f);                 // Verifie si P est générateur de F_q*


#endif // CORPSFINI_INCLUDED