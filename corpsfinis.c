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

int inverse_mod(int a, int p){
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
void initp_polynull(polynome* P){
    if (P->coeff != NULL){
        free(P->coeff);
    }
    P->coeff = NULL;
    P->degre = -1;
    return;
}

int initp_copie(polynome* P, polynome* Q){
    // Si P est deja initialisé, on le vide
    if(P->degre != NULL){
        free(P->coeff);
    }
    if(Q->degre < -1){
        return 1;   // Q n'est pas correctement défini
    }
    else if(Q->degre == -1){
        init_polynull(P);
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

int initp_monome(polynome* P, int coeff, int exp){
    if(exp < 0){
        return 1;   // Choisir un exposant positif
    }
    if(P->coeff != NULL){
        free(P->coeff);
    }
    if(coeff == 0){
        initp_polynull(P);
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
void initp_polynome(polynome* P, int* coeff, int degre){
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

void modulo_transfo(polynome* P, int p){
    if ( P->degre == -1 ){
        return; 
    }
    int d = -1;
    for (int i = P->degre; i >= 0; i--){
        if ( (P->coeff[i])%p != 0 ){
            d = i;
            break;
        }
    }
    if ( d == -1 ){
        free(P->coeff);
        P->degre = -1;
        P->coeff = NULL;
        return;
    }
    P->degre = d;
    P->coeff = realloc(P->coeff, (d+1)*sizeof(int));
    for (int i = 0; i <= P->degre; i++){
        P->coeff[i] = modulo(P->coeff[i], p);
    }
    return;
}

void scalaire_mod(int n, polynome* P, int p){
    scalaire(n, P);
    modulo_transfo(P, p);
    return;
}

void unitaire(polynome* P, int p){
    int C = P->coeff[P->degre];     // coeff dominant de P
    int C_inv = inverse_mod(C, p);  // son inverse
    scalaire(C_inv, P);
    modulo_transfo(P, p);
    return;
}

polynome difference_etendu(polynome* A, polynome* Q, polynome* B){
    polynome S, S_inter;
    S = multiplication(Q, B);
    S_inter = soustraction(A, &S);
    vider(&S);
    S = copie(&S_inter);
    vider(&S_inter);
    return S;
}

polynome difference_etendu_mod(polynome* A, polynome* Q, polynome* B, int p){
    polynome S, S_inter;
    S = multiplication(Q, B);
    S_inter = soustraction(A, &S);
    vider(&S);
    S = copie(&S_inter);
    vider(&S_inter);
    modulo_transfo(&S, p);
    return S;
}

void ajout(polynome* A, polynome* B){
    polynome S = addition(A, B);
    vider(A);
    *A = copie(&S);
    vider(&S);
    return;
}

static void swap(polynome* P, polynome* Q){
    polynome R = copie(P);
    vider(P);
    *P = copie(Q);
    vider(Q);
    *Q = copie(&R);
    vider(&R);
}

// Surjection Z[X] ->> F_p[X] ->> F_p[X]/(f) 
polynome surjection(polynome* P, int p, polynome* f){
    modulo_transfo(P, p);
    polynome P_mod = division_euclid(P, f, p, 1);
    return P_mod;
}

// Surjection F_p[X] ->> F_p[X]/(f)  
polynome demi_surjection(polynome* P, int p, polynome* f){
    polynome P_mod = division_euclid(P, f, p, 1);
    return P_mod;
}

void scalaire_Fq(int n, polynome* P, int p, polynome* f){
    scalaire(n, P);
    polynome P_inter = surjection(P, p, f);
    vider(P);
    *P = copie(&P_inter);
    vider(&P_inter);
    return;
}

static void flip(polynome* r0, polynome* r1, polynome* r2){
    vider(r0);
    *r0 = copie(r1);
    vider(r1);
    *r1 = copie(r2);
    vider(r2);
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

int cf_cdp(polynome* P){
    if (P->degre == -1){
        return 0;
    }
    else {
        return P->coeff[P->degre];
    }
}

void cf_mulp_int(polynome* P, polynome* A, int n){
    initp_polynull(P);
    if (n == 0 || A->degre == -1){
        return;
    }
    else {
        P->degre = A->degre;
        P->coeff = calloc(P->degre + 1, sizeof(int));
        for (int i = 0; i <= A->degre; i++){
            P->coeff[i] = n * A->coeff[i];
        }
        return;
    }
}

void cf_mulp_int_tr(polynome* P, int n){
    if (P->degre == -1)
        return ; // P est nul
    if (n == 0){
        initp_polynull(P);
        return;
    }
    else {
        for (int i = 0; i<= P->degre; i++){
            P->coeff[i] = n * P->coeff[i];
        }
        return;
    }
}

int cf_addp(polynome* P, polynome* A, polynome* B){
    initp_polynull(P);
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

int cf_opposep(polynome* P, polynome* A){
    initp_polynull(P);
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

int cf_subp(polynome* P, polynome* A, polynome* B){
    initp_polynull(P);
    polynome T;
    cf_opposep(&T, B);
    cf_addp(P, A, &T);
    vider(&T);
    return 0;
}

int cf_mulp(polynome* P, polynome* A, polynome* B){
    initp_polynull(P);
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

int cf_puissancep(polynome* P, polynome* A, int exp){
    initp_polynull(P);
    if (exp == 0){
        initp_monome(P, 1, 0);
        return 0;
    }
    else if (exp == 1){
        initp_copie(P, A);
        return 0;
    }
    else {
        initp_copie(P, A);
        polynome T;
        for (int i = 2; i<=exp; i++){
            cf_mulp(&T, P, A);
            initp_copie(P, &T);
        }
        viderp(&T);
        return 0;
    }
}

int cf_redp_int(polynome* P, polynome* A, int p){
    initp_polynull(P);
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

/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                          6. OPERATIONS DANS F_p[X]                                             |   
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|                                                     
 */
polynome addition_mod(polynome* A, polynome* B, int p){
    polynome S = addition(A, B);
    modulo_transfo(&S, p);
    return S;
};

polynome soustraction_mod(polynome* A, polynome* B, int p){
    polynome S = soustraction(A, B);
    modulo_transfo(&S, p);
    return S;
};

polynome multiplication_mod(polynome* A, polynome* B, int p){
    polynome P = multiplication(A, B);
    modulo_transfo(&P, p);
    return P;
}

polynome puissance_mod(polynome* A, int exposant, int p){
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
            inter = multiplication_mod(A, &result, p);
            vider(&result);
            result = copie(&inter);
            vider(&inter);
        }
        return result;
    }
}

polynome division_euclid(polynome* A, polynome* B, int p, int i){
    modulo_transfo(A, p); modulo_transfo(B, p);
    polynome Q, M;
    Q = polynull();
    polynome R = copie(A);
    while ( R.degre >= B->degre ) {
        M = monome( R.coeff[R.degre] * inverse_mod(B->coeff[B->degre], p), R.degre - B->degre );
        ajout(&Q, &M);
        modulo_transfo(&Q, p);
        vider(&M);
        vider(&R);
        R = difference_etendu_mod(A, &Q, B, p);
    }
    if (i == 0){return Q;}
    else {return R;}
}

polynome algo_euclide(polynome* P, polynome* Q, int p){
    if (Q->degre > P->degre){
        swap(P, Q);
    }
    polynome r0 = copie(P);
    polynome r1 = copie(Q);
    polynome r2 = division_euclid(&r0, &r1, p, 1);
    while ( r2.degre != -1 ){
        vider(&r0);
        r0 = copie(&r1);
        vider(&r1);
        r1 = copie(&r2);
        vider(&r2);
        r2 = division_euclid(&r0, &r1, p, 1);
    }
    vider(&r0); vider(&r2);
    unitaire(&r1, p);
    return r1;
}

// Ici le degré de P est supposé supérieur à celui de Q
polynome algo_euclide_etendu(polynome* P, polynome* Q, int p, int i){
    polynome r0 = copie(P);
    polynome r1 = copie(Q);
    polynome u0 = monome(1, 0); polynome u1 = monome(0, 0);
    polynome v0 = monome(0, 0); polynome v1 = monome(1, 0);
    polynome r2 = division_euclid(&r0, &r1, p, 1);
    polynome q = division_euclid(&r0, &r1, p, 0);
    polynome u2 = difference_etendu_mod(&u0, &q, &u1, p);
    polynome v2 = difference_etendu_mod(&v0, &q, &v1, p);
    while ( r2.degre != -1 ){
        flip(&r0, &r1, &r2);
        r2 = division_euclid(&r0, &r1, p, 1);
        vider(&q);
        q = division_euclid(&r0, &r1, p, 0);
        flip(&u0, &u1, &u2);
        u2 = difference_etendu_mod(&u0, &q, &u1, p);
        flip(&v0, &v1, &v2);
        v2 = difference_etendu_mod(&v0, &q, &v1, p);
    }
    int C = r1.coeff[r1.degre];
    int C_inv = inverse_mod(C, p);
    vider(&r0); vider(&r1); vider(&r2);
    vider(&u0); vider(&u2);
    vider(&v0); vider(&v2);
    vider(&q);
    if ( i == 0 ){
        vider(&v1);
        scalaire(C_inv, &u1);
        modulo_transfo(&u1, p);
        return u1;
    }
    else {
        vider(&u1);
        scalaire(C_inv, &v1);
        modulo_transfo(&v1, p);
        return v1;
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











