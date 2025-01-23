#include "corpsfinis.h"

/*  SOMMAIRE :
    1. Maths
    2. Outils
    3. Initialisation de polynômes
    4. Gestion
    5. Opérations dans Z[X]
    6. Opérations dans F_p[X]
    7. Opérations dans F_q = F_p[X]/(f)
 */

/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                              1. OUTILS MATHS                                                   |   
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
    |                                      3. INITIALISATION DE POLYNOMES                                            |   
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|                                                     
 */
void cf_initp_polynull(polynome* P){
    if (P->coeff != NULL){
        free(P->coeff);
    }
    P->coeff = NULL;
    P->degre = -1;
    return;
}

int cf_initp_copie(polynome* P, polynome* Q){
    // Si P est deja initialisé, on le vide
    if(P->degre != NULL){
        free(P->coeff);
    }
    if(Q->degre < -1){
        return 1;   // Q n'est pas correctement défini
    }
    else if(Q->degre == -1){
        cf_init_polynull(P);
        return 0;   // P <- 0
    }
    else{
        P->degre = Q->degre;
        P->coeff = calloc(Q->degre+1, sizeof(int));
        for(int i=0; i < Q->degre+1; i++){
            P->coeff[i] = Q->coeff[i];
        }
        return 0;   // P <- Q
    }
}

int cf_initp_monome(polynome* P, int coeff, int exp){
    if(exp < 0){
        return 1;   // Choisir un exposant positif
    }
    if(P->coeff != NULL){
        free(P->coeff);
    }
    if(coeff == 0){
        cf_initp_polynull(P);
        return 0;   // P <- 0
    }
    else{
        P->degre = exp;
        P->coeff = calloc(exp + 1, sizeof(int));
        for (int i = 0; i < exp; i++){
            P.coeff[i] = 0;
            }
        P.coeff[exp] = coeff;
        return 0;   // P <- coeff*X^exp
    }
}

// Initialisation d'un polynôme à partir d'une liste contenant les coefficients souhaités
void cf_initp_polynome(polynome* P, int* coeff, int degre){
    if(degre < 0){
        return 1;   // Choisir un degré positif (ou utiliser initp_polynull pour degre = -1)
    }
    // Vider P si il est déjà initialisé
    if(P->coeff != NULL){
        free(P->coeff);
    }
    // Début de l'initialisation
    P->degre = degre;
    P->coeff = calloc(degre+1, sizeof(int));
    for(int i = 0; i <= degre; i++){
        P->coeff[i] = coeff[i];
    }
    return 0;   // Initialisation réussie
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

void cf_viderp(polynome* P){
    if (P->coeff != NULL){
        free(P->coeff);
        P->degre = -1;
        P->coeff = NULL;
    }
}


/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                        4. OUTILS POUR LES POLYNOMES                                            |   
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|                                                     
 */
void cf_swapp(polynome* P, polynome* Q){
    polynome T;
    cf_initp_copie(&T, P);
    cf_initp_copie(P, Q);
    cf_initp_copie(Q, &T);
    cf_viderp(&T);
    return;
}

void cf_flipp(polynome* r0, polynome* r1, polynome* r2){
    cf_initp_copie(r0, r1);
    cf_initp_copie(r1, r2);
    cf_viderp(r2);
    return;
}

// Surjection Z[X] ->> F_p[X] ->> F_p[X]/(f) 
polynome surjection(polynome* P, int p, polynome* f){
    cf_redp_int_tr(P, p);
    polynome P_mod = division_euclid(P, f, p, 1);
    return P_mod;
}

// Surjection F_p[X] ->> F_p[X]/(f)  
polynome demi_surjection(polynome* P, int p, polynome* f){
    polynome P_mod = division_euclid(P, f, p, 1);
    return P_mod;
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

int cf_cdp(polynome* P){
    if (P->degre == -1){
        return 0;
    }
    else {
        return P->coeff[P->degre];
    }
}

void cf_mulp_int(polynome* P, polynome* A, int n){
    cf_initp_polynull(P);
    if (n == 0 || A->degre == -1){
        return 0;
    }
    else {
        P->degre = A->degre;
        P->coeff = calloc(P->degre + 1, sizeof(int));
        for (int i = 0; i <= A->degre; i++){
            P->coeff[i] = n * A->coeff[i];
        }
        return 0;
    }
}

void cf_mulp_int_tr(polynome* P, int n){
    if (P->degre == -1)
        return ; // P est nul
    if (n == 0){
        cf_initp_polynull(P);
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
    cf_initp_polynull(P);
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
    cf_initp_polynull(P);
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
    cf_initp_polynull(P);
    polynome T;
    cf_opposep(&T, B);
    cf_addp(P, A, &T);
    cf_viderp(&T);
    return 0;
}

int cf_subp_tr(polynome* P, polynome* A){
    polynome T;
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
    cf_mulp_int_mod_tr(P, C_inv);
    return 0;
}

int cf_mulp(polynome* P, polynome* A, polynome* B){
    cf_initp_polynull(P);
    if ( (A->degre == -1) || (B->degre == -1) ){
        return 0;
    }
    else {
        P->degre = A->degre + B->degre;
        P->coeff = calloc(P.degre + 1, sizeof(int));
        int C;
        for (int k = 0; k <= P.degre; k++){
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
        cf_initp_polynull(P);
        return 0;
    }
    else {
        polynome T;
        cf_copie(&T, P);
        cf_mulp(P, &T, A);
        cf_viderp(&T);
        return 0;
    }
}

int cf_diffetnd(polynome* P, polynome* A, polynome* Q, polynome* B){
    polynome T;
    cf_mulp(&T, Q, B);
    cf_subp(P, A, &T);
    cf_viderp(&T);
    return 0;
}

int cf_puissancep(polynome* P, polynome* A, int exp){
    cf_initp_polynull(P);
    if (exp == 0){
        cf_initp_monome(P, 1, 0);
        return 0;
    }
    else if (exp == 1){
        cf_initp_copie(P, A);
        return 0;
    }
    else {
        cf_initp_copie(P, A);
        polynome T;
        for (int i = 2; i<=exp; i++){
            cf_mulp(&T, P, A);
            cf_initp_copie(P, &T);
        }
        cf_viderp(&T);
        return 0;
    }
}

int cf_redp_int(polynome* P, polynome* A, int p){
    cf_initp_polynull(P);
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
        cf_initp_polynull(P);
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
        initp_monome(P, 1, 0);
        return 0;
    }
    initp_copie(P, A);
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
    cf_initp_copie(&R, A);
    while ( R.degre >= B->degre ) {
        cf_initp_monome(&M, R.coeff[R.degre] * cf_inv_mod(B->coeff[B->degre], p), R.degre - B->degre );
        cf_addp_tr(&Q, &M);
        cf_redp_int_tr(&Q, p);
        cf_viderp(&M);
        cf_viderp(&R);
        cf_diffetnd_mod(&R, A, &Q, B, p);
    }
    if (i == 0){
        cf_viderp(&R);
        cf_initp_copie(P, &Q);
        cf_viderp(&Q);
        return 0;
    }
    else {
        cf_viderp(&Q);
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
    initp_copie(&r0, A);
    initp_copie(&r1, B);
    int result = cf_divp_mod(&r2, &r0, &r1, p, 1);
    while ( r2.degre != -1 ){
        cf_initp_copie(&r0, &r1);
        cf_initp_copie(&r1, &r2);
        result = cf_divp_mod(&r2, &r0, &r1, p, 1);
    }
    viderp(&r0); 
    viderp(&r2);
    unitaire(&r1, p);
    cf_initp_copie(P, &r1);
    viderp(&r1);
    return result;
}

// Ici le degré de A est supposé supérieur à celui de B
int cf_bezoutp_mod(polynome* P, polynome* A, polynome* B, int p, int i){
    polynome r0, r1, r2;
    initp_copie(&r0, A);
    initp_copie(&r1, B);

    polynome u0, u1, u2;
    polynome v0, v1, v2;
    initp_monome(&u0, 1, 0);
    initp_monome(&u1, 0, 0);
    initp_monome(&v0, 0, 0);
    initp_monome(&v1, 1, 0);

    cf_divp_mod(&r2, &r0, &r1, p, 1);   // r2 = r0 - r1*q
    polynome q;
    cf_divp_mod(&q, &r0, &r1, p, 0);
    
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
        cf_initp_copie(P, &u1);
        cf_viderp(&u1);
        return 0;
    }
    else {
        cf_viderp(&u1);
        cf_mulp_int_tr(&v1, C_inv);
        cf_redp_int_tr(&v1, p);
        cf_initp_copie(P, &v1);
        cf_viderp(&v1);
        return 0;
    }
}


/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                        7. OPERATIONS DANS F_q = F_p[X]/(f)                                     |   
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|                                                     
 */
int cardinal(int p, polynome* f){
    return puissance(p, f->degre);
}

polynome addition_Fq(polynome* A, polynome* B, int p, polynome* f){
    polynome S_inter = addition_mod(A, B, p);
    polynome S = demi_surjection(&S_inter, p, f);
    vider(&S_inter);
    return S;
}

polynome multiplication_Fq(polynome* A, polynome* B, int p, polynome* f){
    polynome P_inter = multiplication_mod(A, B, p);
    polynome P = demi_surjection(&P_inter, p, f);
    vider(&P_inter);
    return P;
}

polynome puissance_Fq(polynome* A, int exposant, int p, polynome* f){
    polynome result;
    if (exposant == 0){
        result = monome(1, 0);
        return result;
    }
    else if (exposant == 1){
        result = copie(A);
        return result;
    }
    else {
        result = copie(A);
        polynome inter;
        for (int i = 2; i<=exposant; i++){
            inter = multiplication_Fq(A, &result, p, f);
            vider(&result);
            result = copie(&inter);
            vider(&inter);
        }
        return result;
    }
}

polynome inverse(polynome* P, int p, polynome* f){
    polynome P_inter = surjection(P, p, f);
    polynome u = algo_euclide_etendu(f, &P_inter, p, 1);
    vider(&P_inter);
    return u;
}

polynome division(polynome* A, polynome* B, int p, polynome* f){
    polynome B_inv = inverse(B, p, f);
    polynome D = multiplication_Fq(A, &B_inv, p, f);
    vider(&B_inv);
    return D;
}

int ordre(polynome* P, int p, polynome* f){
    if (P->degre == -1) { return 0; }
    if (P->degre == 0 && P->coeff[0] == 1) { return 1; }
    int i = 2;
    polynome PDT = copie(P);
    polynome PDTT = multiplication_Fq(&PDT, P, p, f);
    while (PDTT.degre != 0 || PDTT.coeff[0] !=1){
        vider(&PDT);
        PDT = copie(&PDTT);
        vider(&PDTT);
        PDTT = multiplication_Fq(&PDT, P, p, f);
        i++;
    }
    return i;
}

int verif_generateur(polynome* P, int p, polynome* f){
    int i = ordre(P, p, f);
    if (i == puissance(p, f->degre) - 1) {
        return 1;
    }
    else {
        return 0;
    }
}











