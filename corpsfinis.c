#include "corpsfinis.h"

/*  Lib de travail polynomes/corps finis
    @nougzarm

    Sommaire:
        (1. Structures - voir .h)
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
    |                                         2. OUTILS NOMBRES ENTIERS                                              |
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
 */
int max(int a, int b){
    if (a>b)
        return a;
    else
        return b;
}

int min(int a, int b){
    if (a>b)
        return b;
    else
        return a;
}

int modulo(int a, int b){
    return a%b + (b * (a%b<0));
}

int puissance(int m, int e){
    if ( e==0 ){
        return 1;
    }
    else if (( e%2==0 )&( e!=0 )){
        return puissance(m*m, e/2);
    }
    else {
        return m*puissance(m*m, (e-1)/2);
    }
}

int puissance_modulo(int m, int e, int p){
    int result = 1;
    for (int i = 1; i <= e; i++){
        result = modulo(result, p)*m;
    }
    return modulo(result, p);
}

int cf_inv_mod(int a, int p){
    int r0 = p, r1 = a, r2 = r0%r1;
    int v0 = 0, v1 = 1, v2 = v0 - (r0/r1)*v1;
    while ( r2 != 0 ){
        r0 = r1;
        r1 = r2;
        r2 = modulo(r0, r1);
        v0 = v1;
        v1 = v2;
        v2 = v0 - (r0/r1)*v1;
    }
    v1 = modulo(v1, p);
    return v1;
}


/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                 3. INITIALISATION(/DEFINITION) DE POLYNOMES                                    |
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
*/
int cf_initp_polynull(polynome* P){
    P->coeff = NULL;
    P->degre = -1;
    return 0;
}

int cf_setp_polynull(polynome* P){
    if (P->coeff != NULL){
        free(P->coeff);
    }
    return cf_initp_polynull(P);
}

int cf_initp_copie(polynome* P, polynome* Q){
    if(Q->degre == -1){
        return cf_initp_polynull(P);
    }
    P->degre = Q->degre;
    P->coeff = (int*)calloc(Q->degre+1, sizeof(int));
    for(int i=0; i < Q->degre+1; i++){
        P->coeff[i] = Q->coeff[i];
    }
    return 0;   // P <- Q
}

int cf_setp_copie(polynome* P, polynome* Q){
    if(P->coeff != NULL){
        free(P->coeff);
    }
    return cf_initp_copie(P, Q);
}

int cf_initp_monome(polynome* P, int coeff, int exp){
    if(coeff == 0){
        cf_initp_polynull(P);
        return 0;
    }
    P->degre = exp;
    P->coeff = calloc(exp + 1, sizeof(int));
    for (int i = 0; i < exp; i++){
        P->coeff[i] = 0;
    }
    P->coeff[exp] = coeff;
    return 0;
}

int cf_setp_monome(polynome* P, int coeff, int exp){
    if(P->coeff != NULL){
        free(P->coeff);
    }
    return cf_initp_monome(P, coeff, exp);
}

int cf_initp_liste(polynome* P, int* coeff, int degre){
    P->degre = degre;
    P->coeff = calloc(degre+1, sizeof(int));
    for(int i = 0; i <= degre; i++){
        P->coeff[i] = coeff[i];
    }
    return 0;
}

int cf_setp_liste(polynome* P, int* coeff, int degre){
    if(P->coeff != NULL){
        free(P->coeff);
    }
    return cf_initp_liste(P, coeff, degre);
}


/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                           4. GESTION DE POLYNOMES                                              |
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
*/

int cf_afficherp(polynome* P){
    if (P->degre == -1){
        printf("0");
        return 0;
    }
    else {
        int coeff = P->coeff[P->degre];     // Coefficient dominant de P
        // Affichage du monôme dominant
        if(coeff == 1){
            printf("X^%d", P->degre);
        }
        else{
            printf("%dX^%d", coeff, P->degre);
        }
        // Affichage des autres monômes
        for (int i = 1; i <= P->degre; i++){
            coeff = P->coeff[P->degre - i]; 
            if (coeff != 0 && coeff == 1){
                printf(" + X^%d", P->degre - i);
            }
            else if(coeff != 0 && coeff != 1){
                printf(" + %dX^%d", coeff, P->degre - i);
            }
        }
        return 0;
    }
}

int cf_viderp(polynome* P){
    return cf_setp_polynull(P);
}

void cf_swapp(polynome* P, polynome* Q){
    polynome T;
    cf_initp_copie(&T, P);
    cf_setp_copie(P, Q);
    cf_setp_copie(Q, &T);
    cf_viderp(&T);
    return;
}

void cf_flipp(polynome* r0, polynome* r1, polynome* r2){
    cf_setp_copie(r0, r1);
    cf_setp_copie(r1, r2);
    cf_viderp(r2);
    return;
}


/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                           5. OPERATIONS DANS Z[X]                                              |   
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|                                                     
 */
int cf_comparp(polynome* P, polynome* A){
    if (P->degre != A->degre){
        return 0;
    }
    else {
        for (int i = 0; i <= P->degre; i++){
            if (P->coeff[i] != A->coeff[i]){
                return 0;
            }
        }
        return 1;
    }
}

int cf_comparp_X(polynome* P){
    if (P->degre != 1){
        return 0;
    }
    else if (P->coeff[0] == 0 && P->coeff[1] == 1){
        return 1;
    }
    else {
        return 0;
    }
}

int cf_cdp(polynome* P){
    if (P->degre == -1){
        return 0;
    }
    else {
        return P->coeff[P->degre];
    }
}

int cf_mulp_int(polynome* P, polynome* A, int n){
    cf_setp_polynull(P);
    if (n == 0 || A->degre == -1){
        return 0;
    }
    else {
        P->degre = A->degre;
        P->coeff = (int*)calloc(P->degre + 1, sizeof(int));
        for (int i = 0; i <= A->degre; i++){
            P->coeff[i] = n * A->coeff[i];
        }
        return 0;
    }
}

int cf_mulp_int_tr(polynome* P, int n){
    if (P->degre == -1)
        return 0;   // P est nul donc pas besoin
    else if (n == 0){
        cf_setp_polynull(P);
        return 0;
    }
    else {
        for (int i = 0; i<= P->degre; i++){
            P->coeff[i] = n * P->coeff[i];
        }
        return 0;
    }
}

int cf_addp(polynome* P, polynome* A, polynome* B){
    cf_setp_polynull(P);
    if (A->degre != B->degre){
        P->degre = max(A->degre, B->degre);
        P->coeff = calloc(P->degre + 1, sizeof(int));
        int mini = min(A->degre, B->degre);
        for (int i = 0; i <= mini; i++){
            P->coeff[i] = A->coeff[i] + B->coeff[i];
        }
        if (A->degre > B->degre){
            for (int i = mini + 1; i <= P->degre; i++){
                P->coeff[i] = A->coeff[i];
            }
        }
        else {
            for (int i = mini + 1; i <= P->degre; i++){
                P->coeff[i] = B->coeff[i];
            }
        }
        return 0;
    }
    else {
        int n = A->degre;
        P->degre = -1;
        for (int i = n; i >= 0; i--){
            if (A->coeff[i] + B->coeff[i] != 0) {
                P->degre = i;
                break;
            }
        }
        if (P->degre == -1){
            P->coeff = NULL;
            return 0;
        }
        P->coeff = calloc(P->degre + 1, sizeof(int));
        for (int i = 0; i <= P->degre; i++){
            P->coeff[i] = A->coeff[i] + B->coeff[i];
        }
        return 0;
    }
}

int cf_addp_tr(polynome* P, polynome* A){
    if (P->degre != A->degre){
        int degP = max(P->degre, A->degre);
        int min_ = min(P->degre, A->degre);
        int* coeffsP = calloc(degP + 1, sizeof(int));
        for (int i = 0; i <= min_; i++){
            coeffsP[i] = P->coeff[i] + A->coeff[i];
        }
        if (P->degre > A->degre){
            for (int i = min_ + 1; i <= degP; i++){
                coeffsP[i] = P->coeff[i];
            }
        }
        else {
            for (int i = min_ + 1; i <= degP; i++){
                coeffsP[i] = A->coeff[i];
            }
        }
        free(P->coeff);
        P->coeff = coeffsP;
        P->degre = degP;
        return 0;
    }
    else{
        int degP = -1;
        for (int i = P->degre; i >= 0; i--){
            if (P->coeff[i] + A->coeff[i] != 0){
                degP = i;
                break;
            }
        }
        if (degP == -1){
            if(P->coeff != NULL){
                free(P->coeff);
            }
            P->coeff = NULL;
            P->degre = -1;
            return 0;
        }
        else {
            int* coeffsP = calloc(degP + 1, sizeof(int));
            for (int i = 0; i <= degP; i++){
                coeffsP[i] = P->coeff[i] + A->coeff[i];
            }
            if (P->coeff != NULL){
                free(P->coeff);
            }
            P->coeff = coeffsP;
            P->degre = degP;
            return 0;
        }
    }
}

int cf_opposep(polynome* P, polynome* A){
    cf_setp_polynull(P);
    if(A->degre == -1){
        return 0;
    }
    P->degre = A->degre;
    P->coeff = calloc(P->degre + 1, sizeof(int));
    for (int i = 0; i <= P->degre; i++){
        P->coeff[i] = - A->coeff[i];
    }
    return 0;
}

int cf_opposep_tr(polynome* P){
    if (P->degre == -1){
        return 0;
    }
    for (int i = 0; i <= P->degre; i++){
        P->coeff[i] = - P->coeff[i];
    }
    return 0;
}

int cf_subp(polynome* P, polynome* A, polynome* B){
    cf_setp_polynull(P);
    polynome T;
    cf_initp_polynull(&T);
    cf_opposep(&T, B);
    cf_addp(P, A, &T);
    cf_viderp(&T);
    return 0;
}

int cf_subp_tr(polynome* P, polynome* A){
    polynome T;
    cf_initp_polynull(&T);
    cf_opposep(&T, A);
    cf_addp_tr(P, &T);
    cf_viderp(&T);
    return 0;
}

int cf_mulp_int_mod(polynome* P, polynome* A, int n, int p){
    cf_mulp_int(P, A, n);
    cf_redp_int_tr(P, p);
    return 0;
}

int cf_mulp_int_mod_tr(polynome* P, int n, int p){
    cf_mulp_int_tr(P, n);
    cf_redp_int_tr(P, p);
    return 0;
}

int cf_unitp_mod(polynome* P, int p){
    int C = P->coeff[P->degre];     // coeff dominant de P
    int C_inv = cf_inv_mod(C, p);   // son inverse
    cf_mulp_int_mod_tr(P, C_inv, p);
    return 0;
}

int cf_mulp(polynome* P, polynome* A, polynome* B){
    cf_setp_polynull(P);
    if ( (A->degre == -1) || (B->degre == -1) ){
        return 0;
    }
    else {
        P->degre = A->degre + B->degre;
        P->coeff = calloc(P->degre + 1, sizeof(int));
        int C;
        for (int k = 0; k <= P->degre; k++){
            C = 0;
            for (int i = 0; i <= k; i++){
                if ( (i <= A->degre) && (k-i <= B->degre) ){
                    C = C + (A->coeff[i])*(B->coeff[k-i]);
                }
            }
            P->coeff[k] = C;
        }
        return 0;
    }
}

int cf_mulp_tr(polynome* P, polynome* A){
    if ( (P->degre == -1) || (A->degre == -1) ){
        cf_setp_polynull(P);
        return 0;
    }
    else {
        polynome T;
        cf_initp_copie(&T, P);
        cf_mulp(P, &T, A);
        cf_viderp(&T);
        return 0;
    }
}

int cf_diffetnd(polynome* P, polynome* A, polynome* Q, polynome* B){
    cf_initp_polynull(P);
    polynome T;
    cf_initp_polynull(&T);
    cf_mulp(&T, Q, B);
    cf_subp(P, A, &T);
    cf_viderp(&T);
    return 0;
}

int cf_puissancep(polynome* P, polynome* A, int exp){
    cf_setp_polynull(P);
    if (exp == 0){
        cf_setp_monome(P, 1, 0);
        return 0;
    }
    else if (exp == 1){
        cf_setp_copie(P, A);
        return 0;
    }
    else {
        cf_setp_copie(P, A);
        for (int i = 2; i<=exp; i++){
            cf_mulp_tr(P, A);
        }
        return 0;
    }
}

int cf_redp_int(polynome* P, polynome* A, int p){
    cf_setp_polynull(P);
    for (int i = A->degre; i >= 0; i--){
        if ( (A->degre)%p != 0 ){
            P->degre = i;
            break;
        }
    }
    if ( P->degre == -1 ){
        return 0;
    }
    else {
        P->coeff = calloc(P->degre + 1, sizeof(int));
        for (int i = 0; i <= P->degre; i++){
            P->coeff[i] = modulo(A->coeff[i], p);
        }
        return 0;
    }
}

int cf_redp_int_tr(polynome* P, int p){
    if (P->degre == -1){
        return 0; 
    }
    int d = -1;
    for (int i = P->degre; i >= 0; i--){
        if ( (P->coeff[i])%p != 0 ){
            d = i;
            break;
        }
    }
    if ( d == -1 ){
        cf_setp_polynull(P);
        return 0;
    }
    P->degre = d;
    P->coeff = realloc(P->coeff, (d+1)*sizeof(int));
    for (int i = 0; i <= P->degre; i++){
        P->coeff[i] = modulo(P->coeff[i], p);
    }
    return 0;
}


/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                          6. OPERATIONS DANS F_p[X]                                             |   
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|                                                     
 */
int cf_addp_mod(polynome* P, polynome* A, polynome* B, int p){
    int result = cf_addp(P, A, B);
    cf_redp_int_tr(P, p);
    return result;
};

int cf_addp_mod_tr(polynome* P, polynome* A, int p){
    int result = cf_addp_tr(P, A);
    cf_redp_int_tr(P, p);
    return result;
};

int cf_opposep_mod(polynome* P, polynome* A, int p){
    int result = cf_opposep(P, A);
    cf_redp_int_tr(P, p);
    return result;
};

int cf_opposep_mod_tr(polynome* P, int p){
    int result = cf_opposep_tr(P);
    cf_redp_int_tr(P, p);
    return result;
};

int cf_subp_mod(polynome* P, polynome* A, polynome* B, int p){
    int result = cf_subp(P, A, B);
    cf_redp_int_tr(P, p);
    return result;
};

int cf_subp_mod_tr(polynome* P, polynome* A, int p){
    int result = cf_subp_tr(P, A);
    cf_redp_int_tr(P, p);
    return result;
};

int cf_mulp_mod(polynome* P, polynome* A, polynome* B, int p){
    int result = cf_mulp(P, A, B);
    cf_redp_int_tr(P, p);
    return result;
}

int cf_mulp_mod_tr(polynome* P, polynome* A, int p){
    int result = cf_mulp_tr(P, A);
    cf_redp_int_tr(P, p);
    return result;
}

int cf_diffetnd_mod(polynome* P, polynome* A, polynome* Q, polynome* B, int p){
    cf_diffetnd(P, A, Q, B);
    cf_redp_int_tr(P, p);
    return 0;
}

int cf_puissancep_mod(polynome* P, polynome* A, int exp, int p){
    if(exp == 0){
        cf_initp_monome(P, 1, 0);
        return 0;
    }
    cf_setp_copie(P, A);
    if(P->degre == -1){
        return 0;
    }
    int result = 0;
    for (int i = 2; i<=exp; i++){
        result = cf_mulp_mod_tr(P, A, p);
    }
    return result;
}

int cf_divp_mod(polynome* P, polynome* A, polynome* B, int p, int i){
    cf_redp_int_tr(A, p);
    cf_redp_int_tr(B, p);
    polynome Q, M, R;
    cf_initp_polynull(&Q);
    cf_initp_polynull(&M);
    cf_initp_copie(&R, A);
    while ( R.degre >= B->degre ) {
        cf_setp_monome(&M, R.coeff[R.degre] * cf_inv_mod(B->coeff[B->degre], p), R.degre - B->degre );
        cf_addp_tr(&Q, &M);
        cf_redp_int_tr(&Q, p);
        cf_diffetnd_mod(&R, A, &Q, B, p);
    }
    if (i == 0){
        cf_viderp(&R);
        cf_viderp(&M);
        cf_initp_copie(P, &Q);
        cf_viderp(&Q);
        return 0;
    }
    else {
        cf_viderp(&Q);
        cf_viderp(&M);
        cf_initp_copie(P, &R);
        cf_viderp(&R);
        return 0;
    }
}

int cf_pgcdp_mod(polynome* P, polynome* A, polynome* B, int p){
    if (B->degre > A->degre){
        cf_swapp(A, B);
    }
    polynome r0;
    polynome r1;
    polynome r2;
    cf_initp_copie(&r0, A);
    cf_initp_copie(&r1, B);
    int result = cf_divp_mod(&r2, &r0, &r1, p, 1);
    while ( r2.degre != -1 ){
        cf_setp_copie(&r0, &r1);
        cf_setp_copie(&r1, &r2);
        result = cf_divp_mod(&r2, &r0, &r1, p, 1);
    }
    cf_viderp(&r0); 
    cf_viderp(&r2);
    cf_unitp_mod(&r1, p);
    cf_setp_copie(P, &r1);
    cf_viderp(&r1);
    return result;
}

// Ici le degré de A est supposé supérieur à celui de B
int cf_bezoutp_mod(polynome* P, polynome* A, polynome* B, int p, int i){
    polynome r0, r1, r2;
    cf_initp_copie(&r0, A);
    cf_initp_copie(&r1, B);
    cf_initp_polynull(&r2);
    
    polynome u0, u1, u2;
    polynome v0, v1, v2;
    cf_initp_monome(&u0, 1, 0);
    cf_initp_monome(&u1, 0, 0);
    cf_initp_monome(&v0, 0, 0);
    cf_initp_monome(&v1, 1, 0);
    cf_divp_mod(&r2, &r0, &r1, p, 1);   // r2 = r0 - r1*q
    polynome q;
    cf_initp_polynull(&q);
    cf_divp_mod(&q, &r0, &r1, p, 0);
    cf_initp_polynull(&u2);
    cf_initp_polynull(&v2);
    cf_diffetnd_mod(&u2, &u0, &q, &u1, p);
    cf_diffetnd_mod(&v2, &v0, &q, &v1, p);
    while ( r2.degre != -1 ){
        cf_flipp(&r0, &r1, &r2);
        cf_divp_mod(&r2, &r0, &r1, p, 1);
        cf_divp_mod(&q, &r0, &r1, p, 0);
        cf_flipp(&u0, &u1, &u2);
        cf_diffetnd_mod(&u2, &u0, &q, &u1, p);
        cf_flipp(&v0, &v1, &v2);
        cf_diffetnd_mod(&v2, &v0, &q, &v1, p);
    }
    int C = r1.coeff[r1.degre];
    int C_inv = cf_inv_mod(C, p);
    cf_viderp(&r0); cf_viderp(&r1); cf_viderp(&r2);
    cf_viderp(&u0); cf_viderp(&u2);
    cf_viderp(&v0); cf_viderp(&v2);
    cf_viderp(&q);
    if ( i == 0 ){
        cf_viderp(&v1);
        cf_mulp_int_tr(&u1, C_inv);
        cf_redp_int_tr(&u1, p);
        cf_setp_copie(P, &u1);
        cf_viderp(&u1);
        return 0;
    }
    else {
        cf_viderp(&u1);
        cf_mulp_int_tr(&v1, C_inv);
        cf_redp_int_tr(&v1, p);
        cf_setp_copie(P, &v1);
        cf_viderp(&v1);
        return 0;
    }
}

int cf_redp_pol(polynome* P, polynome* A, polynome* f, int p){
    cf_redp_int_tr(A, p);
    cf_divp_mod(P, A, f, p, 1);
    return 0;
}

int cf_redp_pol_tr(polynome* P, polynome* f, int p){
    cf_redp_int_tr(P, p);
    polynome T;
    cf_divp_mod(&T, P, f, p, 1);
    cf_setp_copie(P, &T);
    cf_viderp(&T);
    return 0;
}


/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                 7. INITIALISATION/GESTION CORPS FINIS                                          |
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
*/
int cf_vidercf(corpsfini* F){
    cf_viderp(&F->relation);
    F->car = -1;
    return 0;
}

int cf_initcf_pol(corpsfini* F, int p, polynome* f){
    if (p < 2){
        return 1;   // p doit être un nombre premier
    }
    if (f->degre < 0){
        return 1;   // f ne doit pas être nul
    }
    F->car = p;
    cf_initp_copie(&F->relation, f);
    return 0;   // Initialisation réussie
}

int cf_initcf_int(corpsfini* F, int p){
    if (p < 2){
        return 1;   // p doit être un nombre premier
    }
    F->car = p;
    cf_initp_monome(&F->relation, 1, 1);
    return 0;   // Initialisation réussie
}

int cf_cardinalcf(corpsfini* F){
    return puissance(F->car, F->relation.degre);
}

int cf_comparcf(corpsfini* F, corpsfini* K){
    if (F->car != K->car){
        return 0;
    }
    if (cf_comparp(&F->relation, &F->relation) == 0){
        return 0;
    }
    else{
        return 1;
    }
}


/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                            8. INITIALISATION/GESTION D'ELEMENTS DE CORPS FINIS                                 |
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
*/
int cf_viderEl(element* x){
    cf_viderp(&x->representation);
    x->corps = NULL;
    return 0;
};

int cf_initEl_null(element* x, corpsfini* F){
    x->corps = F;
    cf_initp_polynull(&x->representation);
    return 0;
}

int cf_setEl_null(element* x, corpsfini* F){
    cf_viderEl(x);
    return cf_initEl_null(x, F);
}

int cf_initEl_pol(element* x, corpsfini* F, polynome* P){
    x->corps = F;
    cf_initp_polynull(&x->representation);
    cf_redp_pol(&x->representation, P, &F->relation, F->car);
    return 0;   // Initialisation réussie
}

int cf_setEl_pol(element* x, corpsfini* F, polynome* P){   
    cf_viderEl(x);
    return cf_initEl_pol(x, F, P);
}

int cf_initEl_int(element* x, corpsfini* F, int n){
    x->corps = F;
    cf_initp_polynull(&x->representation);
    cf_initp_monome(&x->representation, n, 1);
    return 0;
}

int cf_setEl_int(element* x, corpsfini* F, int n){
    cf_viderEl(x);
    return cf_initEl_int(x, F, n);
}

int cf_initEl_copie(element* x, element* y){
    x->corps = y->corps;
    cf_initp_polynull(&x->representation);
    cf_initp_copie(&x->representation, &y->representation);
    return 0;
}

int cf_setEl_copie(element* x, element* y){
    cf_viderEl(x);
    return cf_initEl_copie(x, y);
}

int cf_initEl_unite(element* x, corpsfini* F){
    x->corps = F;
    cf_initp_polynull(&x->representation);
    cf_initp_monome(&x->representation, 1, 0);
    return 0;
}

int cf_setEl_unite(element* x, corpsfini* F){
    cf_viderEl(x);
    return cf_initEl_unite(x, F);
}


/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                      8. OPERATIONS DANS UN CORPS FINI F                                        |
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
*/
int cf_addEl(element* x, element* y, element* z){
    if (cf_comparcf(y->corps, z->corps) == 0){
        return 1;   // y et z ne sont pas définis sur le même corps
    }
    cf_viderEl(x);
    // Données du corps fini
    corpsfini* corps = y->corps;
    int p = corps->car;
    polynome* f = &corps->relation;
    // Définition de x
    x->corps = corps;
    cf_addp_mod(&x->representation, &y->representation, &z->representation, p);
    cf_redp_pol_tr(&x->representation, f, p);
    return 0;
}

int cf_addEl_tr(element* x, element* y){
    if (cf_comparcf(x->corps, y->corps) == 0){
        return 1;   // x et y ne sont pas définis sur le même corps
    }
    int p = x->corps->car;
    polynome* f = &x->corps->relation;
    cf_addp_mod_tr(&x->representation, &y->representation, p);
    cf_redp_pol_tr(&x->representation, f, p);
    return 0;
}

int cf_opposeEl(element* x, element* y){
    cf_viderEl(x);
    x->corps = y->corps;
    cf_opposep_mod(&x->representation, &y->representation, x->corps->car);
    return 0;
}

int cf_opposeEl_tr(element* x){
    cf_opposep_mod_tr(&x->representation, x->corps->car);
    return 0;
}

int cf_subEl(element* x, element* y, element* z){
    if (cf_comparcf(y->corps, z->corps) == 0){
        return 1;   // y et z ne sont pas définis sur le même corps
    }
    cf_viderEl(x);
    element t;
    cf_initEl_copie(&t, z);
    cf_opposeEl(&t, z);
    cf_addEl(x, y, &t);
    cf_viderEl(&t);
    return 0;
}

int cf_subel_tr(element* x, element* y){
    if (cf_comparcf(x->corps, y->corps) == 0){
        return 1;   // x et y ne sont pas définis sur le même corps
    }
    element t;
    cf_initEl_copie(&t, y);
    cf_addEl_tr(x, &t);
    cf_viderEl(&t);
    return 0;
}

int cf_mulEl(element* x, element* y, element* z){
    if (cf_comparcf(y->corps, z->corps) == 0){
        return 1;   // y et z ne sont pas définis sur le même corps
    }
    cf_viderEl(x);
    x->corps = y->corps;
    cf_mulp_mod(&x->representation, &y->representation, &z->representation, x->corps->car);
    cf_redp_pol_tr(&x->representation, &(x->corps->relation), x->corps->car);
    return 0;
}

int cf_mulEl_tr(element* x, element* y){
    if (cf_comparcf(x->corps, y->corps) == 0){
        return 1;   // x et y ne sont pas définis sur le même corps
    }
    cf_mulp_mod_tr(&x->representation, &y->representation, x->corps->car);
    cf_redp_pol_tr(&x->representation, &(x->corps->relation), x->corps->car);
    return 0;
}

int cf_puissanceEl(element* x, element* y, int exp){
    cf_viderEl(x);
    cf_initEl_unite(x, y->corps);
    for (int i = 1; i <= exp; i++){
        cf_mulEl_tr(x, y);
    }
    return 0;
}

int cf_invEl(element* x, element* y){
    if (y->representation.degre == -1){
        return 1;   // y est nul donc non inversible
    }
    cf_viderEl(x);
    x->corps = y->corps;
    cf_bezoutp_mod(&x->representation, &x->corps->relation, &y->representation, x->corps->car, 1);
    return 0;
}

int cf_divEl(element* x, element* y, element* z){
    if (z->representation.degre == -1){
        return 1;
    }
    if (cf_comparcf(y->corps, z->corps) == 0){
        return 1;
    }
    cf_invEl(x, z);
    cf_mulEl_tr(x, y);
    return 0;
}

int cf_divEl_tr(element* x, element* y){
    if (y->representation.degre == -1){
        return 1;
    }
    if (cf_comparcf(x->corps, y->corps) == 0){
        return 1;
    }
    element t;
    cf_initEl_null(&t, x->corps);
    cf_invEl(&t, y);
    cf_mulEl_tr(x, &t);
    cf_viderEl(&t);
    return 0;
}

int cf_ordreEl(element* x){
    if (x->representation.degre == -1){
        return 0;
    }
    else if (x->representation.degre == 0 && x->representation.coeff[0] == 1){
        return 1;
    }
    int ord = 2;
    element pdt;
    cf_initEl_null(&pdt, x->corps);
    cf_mulEl(&pdt, x, x);
    while(pdt.representation.degre != 0 || pdt.representation.coeff[0] != 1){
        cf_mulEl_tr(&pdt, x);
        ord++;
    }
    cf_viderEl(&pdt);
    return ord;
}

int cf_verifgenEl(element* x){
    int ord = cf_ordreEl(x);
    if (ord == cf_cardinalcf(x->corps) - 1){
        return 1;
    }
    else {
        return 0;
    }
}