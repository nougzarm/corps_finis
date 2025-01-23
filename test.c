#include "test.h"

/* En fonction de la valeur de choix_test : 

    0 : Affiche la somme de P et Q dans Z[X]
    1 : Affiche la différence P-Q dans Z[X]
    2 : Affiche la somme de P et Q dans F_p[X]
    3 : Affiche la réduction modulo p de P.   i.e  Z[X] ->> F_p[X] 
    4 : Affiche le produit de P et Q dans Z[X]
    5 : Affiche le produit de P et Q dans F_p[X]
    6 : Affiche la division euclidienne de P par Q dans F_p[X]
    7 : Affiche le PGCD de P et Q (dans F_p[X])
    8 : Affiche la formule de Bezout entre P et Q dans F_p[X]
    9 : Affiche la réduction de P et Q dans F_q
    10 : Affiche le produit de P et Q dans F_q
    11 : Affiche l'inverse de P dans F_q
    12 : Affiche le quotient de P par Q dans F_q
    13 : Affiche l'ordre de P dans F_q
    14 : Affiche si P est générateur de F_q
*/

int main() {
    //  CHOIX DU TEST À EFFECTUER ----------------------------------------------------
    int choix_test = 9;

    //  CONFIGURATION ----------------------------------------------------------------
    int p = 3;
    int A[] = {0, 1, 5, 6, 0, 1};   //  polynome P   (Placer le coeff dominant en fin de liste)
    int B[] = {1, 0, 1, 2, 2};      //  polynome Q   
    int C[] = {1, 0, 1};            //  polynome f   (Irréductible dans F_p[X])
    int a = sizeof(A)/sizeof(int)-1; 
    int b = sizeof(B)/sizeof(int)-1; 
    int c = sizeof(C)/sizeof(int)-1;

    //  Initialisation des polynomes -------------------------------------------------
    polynome P, Q, f;

    cf_initp_liste(&P, A, a);
    cf_initp_liste(&Q, B, b);
    cf_initp_liste(&f, C, c);

    //  Initialisation du corps F_p[X]/(f)
    corpsfini F;
    int r0 = cf_initcf_pol(&F, p, &f);
    //  Initialisation des éléments de F


    //  Affichage et déroulement du test ---------------------------------------------
    printf("----------------------------------------------------------------- \n");
    printf("Définitions :\n - Premier p = %d \n", p);
    printf(" - Polynome irréductible dans F_%d[X] : f = ", p); cf_afficherp(&f); 
    printf("\n - P = "); cf_afficherp(&P);
    printf("\n - Q = "); cf_afficherp(&Q); printf("\n\n");

    test(choix_test, &P, &Q, p, &f, &F); 
    printf("----------------------------------------------------------------- \n");
    //  Libération de la mémoire -----------------------------------------------------
    cf_viderp(&P);
    cf_viderp(&Q);
    cf_viderp(&f);

    if(r0 == 0){
        cf_vidercf(&F); // Cas où F a été initialisé
    }
}



void test(int choix_test, polynome* P, polynome* Q, int p, polynome* f, corpsfini* F){

    if(choix_test == 0){
        polynome S;
        cf_initp_polynull(&S);
        cf_addp(&S, P, Q);
        printf("Dans Z[X],  P + Q = "); cf_afficherp(&S); printf("\n");
        cf_viderp(&S);
    }

    else if(choix_test == 1){
        polynome S;
        cf_initp_polynull(&S);
        cf_subp(&S, P, Q);
        printf("Dans Z[X],  P - Q = "); cf_afficherp(&S); printf("\n");
        cf_viderp(&S);
    }

    else if(choix_test == 2){
        polynome S;
        cf_initp_polynull(&S);
        cf_addp_mod(&S, P, Q, p);
        printf("Dans F_p[X],  P + Q = "); cf_afficherp(&S); printf("\n");
        cf_viderp(&S);
    }

    else if(choix_test == 3){
        polynome T;
        cf_initp_polynull(&T);
        cf_redp_int(&T, P, p);
        printf("La réduction de P modulo %d vaut : ", p); cf_afficherp(&T); printf("\n");
    }
 
    else if(choix_test == 4){
        polynome T;
        cf_initp_polynull(&T);  
        cf_mulp(&T, P, Q);
        printf("Dans Z[X],  P*Q = "); cf_afficherp(&T); printf("\n");
        cf_viderp(&T);
    }

    else if(choix_test == 5){
        polynome T;
        cf_initp_polynull(&T);   
        cf_mulp_mod(&T, P, Q, p);
        printf("Dans F_%d[X],  P*Q = ", p); cf_afficherp(&T); printf("\n");
        cf_viderp(&T);
    }

    else if(choix_test == 6){
        polynome q, r;
        cf_initp_polynull(&q); cf_initp_polynull(&r);
        cf_divp_mod(&q, P, Q, p, 0);
        cf_divp_mod(&r, P, Q, p, 1);
        printf("Résultat de la division euclidienne de P par Q (dans Z/%dZ[X]) : \n", p);
        printf(" - Quotient : "); cf_afficherp(&q);
        printf("\n - Reste : "); cf_afficherp(&r); printf("\n");
        cf_viderp(&q); cf_viderp(&r);
    }

    else if(choix_test == 7){
        polynome T;
        cf_initp_polynull(&T);
        cf_pgcdp_mod(&T, P, Q, p);
        printf("Résultat : \n");
        printf("Le PGCD de P et Q (dans Z/%dZ[X]) est : ", p); cf_afficherp(&T); printf("\n");
        cf_viderp(&T);
    }

    else if(choix_test == 8){
        polynome D, U, V;
        cf_initp_polynull(&D); cf_initp_polynull(&U); cf_initp_polynull(&V);
        cf_pgcdp_mod(&D, P, Q, p);
        cf_bezoutp_mod(&U, P, Q, p, 0);
        cf_bezoutp_mod(&V, P, Q, p, 1);
        printf("Résultat : \n");
        printf("La relation de Bezout entre P et Q est : uP + vQ = PGCD(P,Q), où :\n");
        printf(" - u = "); cf_afficherp(&U); printf("\n");
        printf(" - v = "); cf_afficherp(&V); printf("\n");
        printf(" - PGCD(P,Q) = "); cf_afficherp(&D); printf("\n");
        cf_viderp(&D); cf_viderp(&U); cf_viderp(&V); 
    }

    else if(choix_test == 9){
        polynome P_red, Q_red;
        cf_initp_polynull(&P_red); cf_initp_polynull(&Q_red);
        cf_redp_pol(&P_red, P, f, p);
        cf_redp_pol(&Q_red, Q, f, p);
        printf("L'image de P dans F_%d est : ", puissance(p, f->degre)); cf_afficherp(&P_red); printf("\n");
        printf("L'mage de Q dans F_%d est : ", puissance(p, f->degre)); cf_afficherp(&Q_red); printf("\n");
        cf_viderp(&P_red);
        cf_viderp(&Q_red);
    }
/*
    else if(choix_test == 10) {
        polynome T;
        cf_initp_polynull(&T);
        polynome PDT = multiplication_Fq(P, Q, p, f);
        printf("Dans F_%d,  P*Q = ", q); afficher(&PDT); printf("\n");
        vider(&PDT);
    }

    else if(choix_test == 11)
    {
        polynome P_inv = inverse(P, p, f);
        printf("L'inverse de P dans F_%d est : ", q); afficher(&P_inv); printf("\n");
        vider(&P_inv);
    }

    else if(choix_test == 12){
        polynome D = division(P, Q, p, f);
        printf("Dans F_%d,  P/Q = ", q); afficher(&D); printf("\n");
        vider(&D);
    }

    else if(choix_test == 13){
        int ord = ordre(P, p, f);
        printf("Dans F_%d, l'ordre de P est : %d \n", q, ord);
    }

    else if(choix_test == 14){
        int v = verif_generateur(P, p, f);
        if(v == 0){
            printf("P n'est pas générateur de F_%d \n", q);
        }
        else {
            printf("P est un élément générateur de F_%d \n", q);
        }
    } */
}
