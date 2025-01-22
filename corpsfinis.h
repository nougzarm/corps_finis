#ifndef CORPSFINI_INCLUDED
#define CORPSFINI_INCLUDED

#include <stdio.h>
#include <stdlib.h>

/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                                 0. STRUCTURES                                                  |
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
Structures utilisées:
    polynome : permet de définir un polynome (à coefficients entiers)
    corps_fini : permet de définir un corps fini F_q = (Z/pZ[X])/(f) où p = corps_fini.car et f = corps_fini.relation
        f doit donc être irréductible. 
    element : permet de définir un élément d'un corps fini via sa forme polynomiale element.representation, dans
        un corps donné element.corps
 */
typedef struct{
    int* coeff;
    int degre;
} polynome;

typedef struct{
    int car; 
    polynome* relation;
} corps_fini;

typedef struct{
    polynome* representation;
    corps_fini* corps;
} element;

/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                              1. OUTILS MATHS                                                   |
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
 */
int max(int a, int b);          // max(a, b)
int min(int a, int b);          // min(a, b)
int modulo(int a, int b);       // a mod b
int puissance(int a, int e);    // a^e
int puissance_modulo(int m, int e, int p);  // m^e mod p
int inverse_mod(int a, int p);


/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                 2. INITIALISATION(/DEFINITION) DE POLYNOMES                                    |
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
 */
/*  Initialisation d'un polynôme P à 0  */
void initp_polynull(polynome* P);

/*  Initialisation d'un polynôme P à partir de Q   
    Retourne:   1 si l'initalisation a échoué (Q n'est pas correctement défini)
                0 si la copie est réussie     */
int initp_copie(polynome* P);

/*  Initialisation d'un polynome P par le monomme coeff*X^exp
    Retourne:   1 si l'initalisation a échoué (exp n'est pas positif)
                0 si l'initialisation   */
int initp_monome(polynome* P, int coeff, int exp);

/*  Initialisation d'un polynome P à partir d'une liste de coefficients
    - La liste coefficient contient les coefficients (en commençant par le coefficient du degré constant)
    - Si degre = -1 alors P est initialisé à 0   */
void initp_polynome(polynome* P, int* coeff, int degre);


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