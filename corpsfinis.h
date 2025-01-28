#ifndef CORPSFINI_INCLUDED
#define CORPSFINI_INCLUDED

#include <stdio.h>
#include <stdlib.h>

/*  Lib de travail polynomes/corps finis
    @nougzarm
    
    Sommaire:
        1. Structures
        2. Outils nombres entiers
        3. Initialisation de polynomes
        4. Gestion de polynomes
        5. Opérations dans Z[X]
        6. Opérations dans Z/pZ[X]
        7. Initialisation/gestion de corps finis
        8. Initialisation/gestion d'éléments de corps finis
        9. Opérations dans un corps fini
        ...     */

/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                                 1. STRUCTURES                                                  |
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|  
*/    
/*  Structure de polynôme dans Z[X], représenté par ses coefficients ainsi que son degré    */   
typedef struct{
    int* coeff;
    int degre;
} polynome;

/*  Structure de corps fini défini via sa caractéristique car = p ainsi qu'un polynome relation = f irréductible 
    dans F_p[X], afin que F_q = (Z/pZ[X])/(f) soit un corps fini de cardinal q = p^deg(f)    */
typedef struct{
    int car;
    polynome relation;
} corpsfini;

/*  Structure d'élément de F_q, représenté par un polynome de F_p[X] modulo f    */
typedef struct{
    polynome representation;
    corpsfini* corps;
} element;

/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                         2. OUTILS NOMBRES ENTIERS                                              |
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
*/
/*  Fonctions max et min entre deux entiers  */
int max(int a, int b);
int min(int a, int b);

/*  Retourne la réduction POSITIVE de a modulo b     */
int modulo(int a, int b);

/*  Fonctions exponentiations "optimisée" */
int puissance(int a, int e);
int puissance_modulo(int m, int e, int p);

/*  Retourne l'inverse de a dans Z/pZ    */
int inverse_mod(int a, int p);


/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                 3. INITIALISATION(/DEFINITION) DE POLYNOMES                                    |
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
*/
/*  Remarque importante: suivant le cas où le polynome a deja été (manuellement) initialisé ou non, il est
        nécessaire de bien choisir entre initp et setp
        initp : à utiliser pour initialiser un polynome (directement après une déclaration)
        setp : à utiliser pour attribuer une valeur à un polynôme UNIQUEMENT lorsque le polynome a deja
            été initialisé via initp     */

/*  Initialisation d'un polynôme P à 0  */
int cf_initp_polynull(polynome* P);
int cf_setp_polynull(polynome* P);

/*  Initialisation d'un polynôme P à partir d'une copie de Q     */
int cf_initp_copie(polynome* P, polynome* Q);
int cf_setp_copie(polynome* P, polynome* Q);

/*  Initialisation d'un polynome P par le monomme coeff*X^exp   */
int cf_initp_monome(polynome* P, int coeff, int exp);
int cf_setp_monome(polynome* P, int coeff, int exp);

/*  Initialisation d'un polynome P à partir d'une liste de coefficients
        - La liste coefficient contient les coefficients (en commençant par le coefficient du degré constant)
        - Si degre = -1 alors P est initialisé à 0   */
int cf_initp_liste(polynome* P, int* coeff, int degre);
int cf_setp_liste(polynome* P, int* coeff, int degre);


/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                           4. GESTION DE POLYNOMES                                              |
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
*/
/*  Affichage d'un polynome P   */
int cf_afficherp(polynome* P);

/*  Libération de la mémoire occupée par un polynome P (libére la liste des coefficients)
    Remarque: Utiliser cette fonction uniquement si P a été initialisé (p.ex avec une fonction initp)  */
int cf_viderp(polynome* P);

/*  Echange les polynomes P et Q    */
void cf_swapp(polynome* P, polynome* Q);

/*  Stocke: r0 <- r1 et r1 <- r2 (et vide r2)   */
void cf_flipp(polynome* r0, polynome* r1, polynome* r2);

/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                           5. OPERATIONS DANS Z[X]                                              |
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
*/
/*  Compare deux polynomes P et A (leurs coefficients) 
    Retourne 1 si P=A et 0 sinon     */
int cf_comparp(polynome* P, polynome* A);

/*  Compare P au monome X (Retourne 1 si P=X et 0 sinon)    */
int cf_comparp_X(polynome* P);

/*  Retourne le coefficient dominant d'un polynome P    */
int cf_cdp(polynome* P);

/*  Stocker le produit nA dans P (pour le cas A=P voir cf_mulp_int_tr)    */
int cf_mulp_int(polynome* P, polynome* A, int n);

/*  Tranformation P <- nP    */
int cf_mulp_int_tr(polynome* P, int n);

/*  Stocker le polynome somme A+B dans P (pour le cas A=P voir cf_addp_tr)    */
int cf_addp(polynome* P, polynome* A, polynome* B);

/*  Transforamtion P <- P*A     */
int cf_addp_tr(polynome* P, polynome* A);

/*  Stocker le polynome opposé -A dans P    */
int cf_opposep(polynome* P, polynome* A);

/*  Transformation P <- -P    */
int cf_opposep_tr(polynome* P);

/*  Stocker le polynome différence A+B dans P    */
int cf_subp(polynome* P, polynome* A, polynome* B);

/*  Transformation P <- P-A    */
int cf_subp_tr(polynome* P, polynome* A);

/*  Stocker le polynome produit A*B dans P    */
int cf_mulp(polynome* P, polynome* A, polynome* B);

/*  Transformation P <- P*A    */
int cf_mulp_tr(polynome* P, polynome* A);

/*  Stocker le polynome A-Q*B dans P    */
int cf_diffetnd(polynome* P, polynome* A, polynome* Q, polynome* B);

/*  Stocker le polynome A^exp dans P    */
int cf_puissancep(polynome* P, polynome* A, int exp);

/*  Stocker le polynome A mod p dans P    */
int cf_redp_int(polynome* P, polynome* A, int p);

/*  Transormation P <- P mod p    */
int cf_redp_int_tr(polynome* P, int p);


/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                          6. OPERATIONS DANS F_p[X]                                             |
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
*/
/*  Stocker la somme A+B mod p dans P   */
int cf_addp_mod(polynome* P, polynome* A, polynome* B, int p);

/*  Transformation P <- P+A mod p   */
int cf_addp_mod_tr(polynome* P, polynome* A, int p);

/*  Stocker l'opposé -A mod p dans P   */
int cf_opposep_mod(polynome* P, polynome* A, int p);

/*  Transformation P <- -P mod p   */
int cf_opposep_mod_tr(polynome* P, int p);

/*  Stocker la différence A-B mod p dans P   */
int cf_subp_mod(polynome* P, polynome* A, polynome* B, int p);

/*  Transformation P <- P-A mod p   */
int cf_subp_mod_tr(polynome* P, polynome* A, int p);

/*  Stocker le produit nA mod p dans P   */
int cf_mulp_int_mod(polynome* P, polynome* A, int n, int p);

/*  Transformation P <- nP mod p   */
int cf_mulp_int_mod_tr(polynome* P, int n, int p);

/*  Transformation P <- P*(1/cd(P) mod p),  où cd = coeff dominant  */
int cf_unitp_mod(polynome* P, int p);

/*  Stocker le produit A*B mod p dans P   */
int cf_mulp_mod(polynome* P, polynome* A, polynome* B, int p);

/*  Transformation P <- P*A mod p */
int cf_mulp_mod_tr(polynome* P, polynome* A, int p);

/*  Stocker le polynome A-Q*B mod p dans P    */
int cf_diffetnd_mod(polynome* P, polynome* A, polynome* Q, polynome* B, int p);

/*  Stocker la puissance A^exp mod p dans P   */
int cf_puissancep_mod(polynome*P, polynome* A, int exp, int p);

/*  Division euclidienne de A par B dans F_p[X] (i=0: quotient, i=1: reste)
    Stocke dans:    Quotient de la div. si i=0
                    Reste de la div. si i=1     */
int cf_divp_mod(polynome* P, polynome* A, polynome* B, int p, int i);

/*  Stocke PGCD(A, B) mod p dans P  */
int cf_pgcdp_mod(polynome* P, polynome* A, polynome* B, int p);

/*  Algorithme d'Euclide étendu dans F_p[X] (si uA+vB=pgcd(A, B) alors i=0: P <- v, sinon: P <- u)  */
int cf_bezoutp_mod(polynome* P, polynome* A, polynome* B, int p, int i);

/*  Stocke A mod f dans P (dans F_p[X])  */
int cf_redp_pol(polynome* P, polynome* A, polynome* f, int p); 

/*  Transformation P <- P mod f (dans F_p[X])  */
int cf_redp_pol_tr(polynome* P, polynome* f, int p);


/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                 7. INITIALISATION/GESTION CORPS FINIS                                          |
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
*/
/*  Libération de la mémoire occupée par un corps fini F    */
int cf_vidercf(corpsfini* F);

/*  Initialisation d'un corps fini F_q = F_p[X]/(f)    
    Remarques:  - p doit être un nombre premier
                - f doit être un polynome irréductible de F_p[X]
                - Pour définir F_q = F_p = Z/pZ, on peut choisir f = X (voir fonction init suivante)
                - Penser à vider F AINSI que f après avoir fini (cf_vidercf et cf_viderp)
    Retourne:   1 si echec
                0 si réussi    */
int cf_initcf_pol(corpsfini* F, int p, polynome* f);

/*  Initialisation d'un corps fini F_q = F_p = Z/pZ    */
int cf_initcf_int(corpsfini* F, int p);

/*  Retourne le cardinal d'un corps fini F  */
int cf_cardinalcf(corpsfini* F);

/*  Compare les corps finis F et K (retourne 1 si F=K, 0 sinon) */
int cf_comparcf(corpsfini* F, corpsfini* K);

/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                            8. INITIALISATION/GESTION D'ELEMENTS DE CORPS FINIS                                 |
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
*/
/*  Libération de la mémoire occupée par un element x d'un corps fini    */
int cf_viderEl(element* x);

/*  Initialisation d'un élément x = 0 dans F   */
int cf_initEl_null(element* x, corpsfini* F);
int cf_setEl_null(element* x, corpsfini* F);

/*  Initialisation d'un élément de F via sa forme polynomiale P
    Remarque: penser à vider x AINSI que P après avoir fini   */
int cf_initEl_pol(element* x, corpsfini* F, polynome* P);
int cf_setEl_pol(element* x, corpsfini* F, polynome* P);

/*  Initialisation de x=n (utile lorsque F est de la forme Z/pZ)   */
int cf_initEl_int(element* x, corpsfini* F, int n);
int cf_setEl_int(element* x, corpsfini* F, int n);

/*  Initialisation d'un élément de x en y stockant une copie de y   */
int cf_initEl_copie(element* x, element* y);
int cf_setEl_copie(element* x, element* y);

/*  Initialisation de x à 1 dans F  */
int cf_initEl_unite(element* x, corpsfini* F);
int cf_setEl_unite(element* x, corpsfini* F);


/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                   9. OPERATIONS DANS UN CORPS FINI F                                           |
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
*/
/*  Stocker la somme y+z dans x 
    Retourne:   0 si réussi
                1 si echec (p.ex pas définis sur le même corps)*/
int cf_addEl(element* x, element* y, element* z);

/*  Transformation x <- y dans F 
    Retourne:   0 si réussi
                1 si echec (p.ex pas définis sur le même corps) */
int cf_addEl_tr(element* x, element* y);

/*  Stocker l'opposé -y dans x  */
int cf_opposeEl(element* x, element* y);

/*  Transformation x <- -x  */
int cf_opposeEl_tr(element* x);

/*  Stocker la différence y-z dans x 
    Retourne:   0 si réussi
                1 si echec (p.ex pas définis sur le même corps)*/
int cf_subEl(element* x, element* y, element* z);

/*  Transformation x <- x-y
    Retourne:   0 si réussi
                1 si echec (p.ex pas définis sur le même corps)  */
int cf_subEl_tr(element* x, element* y);

/*  Stocker le produit y*z dans x  
    Retourne:   0 si réussi
                1 si echec (p.ex pas définis sur le même corps)     */
int cf_mulEl(element* x, element* y, element* z);

/*  Transformation x <- x*y dans F
    Retourne:   0 si réussi
                1 si echec (p.ex pas définis sur le même corps)  */
int cf_mulEl_tr(element* x, element* y);

/*  Stocker y^exp dans x */
int cf_puissanceEl(element* x, element* y, int exp);

/*  Stocker y^-1 dans x
    Retourne:   0 si réussi
                1 si echec (p.ex y n'est pas inversible)     */
int cf_invEl(element* x, element* y);

/*  Stocker le quotient y/z dans x  
    Retourne:   0 si réussi
                1 si echec (p.ex z non inversible ou != corps) */
int cf_divEl(element* x, element* y, element* z);

/*  Transformation x <- x/y 
    Retourn:    0 si réussi
                1 si echec (p.ex y non inversible ou != corps)  */
int cf_divEl_tr(element* x, element* y);

/*  Retourne l'ordre d'un élément x dans F  */
int cf_ordreEl(element* x);

/*  Vérifie si x est générateur de F* 
    Retourne:   1 si x est générateur
                0 sinon  */
int cf_verifgenEl(element* x);


#endif // CORPSFINI_INCLUDED