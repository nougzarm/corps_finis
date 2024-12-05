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
    |                                                     1. MATHS                                                   |   
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
    |                                                  2. OUTILS                                                     |   
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|                                                     
 */
int conversion_scalaire(polynome* P){
    if (P->degre == 0){
        return P->coeff[0];
    }
    else {
        return 0;
    }
}

void scalaire(int n, polynome* P){
    if (P->degre == -1)
        return;
    else {
        for (int i = 0; i<= P->degre; i++){
            P->coeff[i] = n * P->coeff[i];
        }
        return;
    }
}

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
    int C = P->coeff[P->degre]; // coeff dominant de P
    int C_inv = inverse_mod(C, p); // son inverse
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
    |                                      3. INITIALISATION DE POLYNOMES                                            |   
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|                                                     
 */
polynome polynull(){
    int* coeff = NULL;
    polynome polynome_nul = {coeff, -1};
    return polynome_nul;
}

polynome copie(polynome* P){
    polynome Q;
    Q.degre = P->degre;
    if (Q.degre == -1){
        Q.coeff = NULL;
    }
    else {
        Q.coeff = calloc(Q.degre + 1, sizeof(int));
        for(int i = 0; i <= Q.degre; i++) {
            Q.coeff[i] = P->coeff[i];
        }
    }
    return Q;
}

polynome monome(int coefficient, int exposant){
    polynome P;
    if (coefficient == 0){
        P = polynull();
        return P;
    }
    else {
        P.degre = exposant;
        P.coeff = calloc(exposant + 1, sizeof(int));
        for (int i = 0; i < exposant; i++) { P.coeff[i] = 0; }
        P.coeff[exposant] = coefficient;
        return P;
    }
}

// Initialisation d'un polynôme à partir d'une liste contenant les coefficients souhaités
polynome polynome_init(int* liste, int taille_liste){
    int* coeff = NULL;
    polynome nouveau_polynome = {coeff, -1};
    if(taille_liste < 1){
        return nouveau_polynome;
    }
    else{
        nouveau_polynome.degre = taille_liste - 1;
        nouveau_polynome.coeff = calloc(taille_liste, sizeof(int));
        for(int i = 0; i < taille_liste; i++){
            nouveau_polynome.coeff[nouveau_polynome.degre - i] = liste[i];
        }
        return nouveau_polynome;
    }
}


/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                                   4. GESTION                                                   |   
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|                                                     
 */
void vider(polynome* P){
    if (P->coeff != NULL){
        free(P->coeff);
        P->degre = -1;
        P->coeff = NULL;
    }
}

void afficher(polynome* P){
    if (P->degre == -1){
        printf("0");
        return;
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
        return;
    }
}

/*  |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|
    |                                           5. OPERATIONS DANS Z[X]                                              |   
    |----------------------------------------------------------------------------------------------------------------|
    |----------------------------------------------------------------------------------------------------------------|                                                     
 */
polynome addition(polynome* P, polynome* Q){
    polynome S;
    if (P->degre != Q->degre){
        S.degre = max(P->degre, Q->degre);
        S.coeff = calloc(S.degre + 1, sizeof(int));
        int mini = min(P->degre, Q->degre); 
        for (int i = 0; i <= mini; i++){
            S.coeff[i] = P->coeff[i] + Q->coeff[i];
        }
        if (P->degre > Q->degre){
            for (int i = mini + 1; i <= P->degre; i++){
                S.coeff[i] = P->coeff[i];
            }
        }
        else {
            for (int i = mini + 1; i <= Q->degre; i++){
                S.coeff[i] = Q->coeff[i];
            }
        }
        return S;
    }
    else {
        int n = P->degre;
        S.degre = -1;
        for (int i = n; i >= 0; i--){
            if (P->coeff[i] + Q->coeff[i] != 0) {
                S.degre = i;
                break;
            }
        }
        if (S.degre == -1){
            S.coeff = NULL;
            return S;
        }
        S.coeff = calloc(S.degre + 1, sizeof(int));
        for (int i = 0; i <= S.degre; i++){
            S.coeff[i] = P->coeff[i] + Q->coeff[i];
        }
        return S;
    }
}


polynome oppose(polynome* P){
    polynome Q;
    Q.degre = P->degre;
    if (Q.degre == -1) {
        Q.coeff = NULL;
    }
    else {
        Q.coeff = calloc(Q.degre + 1, sizeof(int));
        for (int i = 0; i <= Q.degre; i++){
            Q.coeff[i] = - P->coeff[i];
        }
    }
    return Q;
}


polynome soustraction(polynome* P, polynome* Q){
    polynome T = oppose(Q);
    return addition(P, &T);
}

polynome multiplication(polynome* A, polynome* B){
    polynome P;
    if ( (A->degre == -1) || (B->degre == -1) ){
        P.degre = -1;
        P.coeff = NULL;
        return P;
    }
    else {
        P.degre = A->degre + B->degre;
        P.coeff = calloc(P.degre + 1, sizeof(int));
        int C;
        for (int k = 0; k <= P.degre; k++){
            C = 0;
            for (int i = 0; i <= k; i++){
                if ( (i <= A->degre) && (k-i <= B->degre) ){
                    C = C + (A->coeff[i])*(B->coeff[k-i]);
                }
            }
            P.coeff[k] = C;
        }
        return P;
    }
}

polynome puissance_polynome(polynome* A, int exposant){
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
            inter = multiplication(A, &result);
            vider(&result);
            result = copie(&inter);
            vider(&inter);
        }
        return result;
    }
}

polynome reduction_modulo(polynome* P, int p){
    polynome Q; 
    Q.degre = -1;
    for (int i = P->degre; i >= 0; i--){
        if ( (P->degre)%p != 0 ){
            Q.degre = i;
            break;
        }
    }
    if ( Q.degre == -1 ){
        Q.coeff = NULL;
        return Q;
    }
    else {
        Q.coeff = calloc(Q.degre + 1, sizeof(int));
        for (int i = 0; i <= Q.degre; i++){
            Q.coeff[i] = modulo(P->coeff[i], p);
        }
        return Q;
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











