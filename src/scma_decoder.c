#include "scma.h" 

// Function prototypes
double log_sum_exp(double *array, int size);

void scmadec(double complex y[K][N], double complex CB[K][M][V], double complex h[K][V][N], double N0, int Nit, double LLR[log2(M)*V][N]) {
    int K = sizeof(CB) / sizeof(CB[0]);  // Number of orthogonal resources
    int M = sizeof(CB[0]) / sizeof(CB[0][0]);  // Number of codewords in each codebook
    int V = sizeof(CB[0][0]) / sizeof(CB[0][0][0]);  // Number of users (layers)
    
    double Noise = 1.0 / N0;
    double F[K][V]; // Factor graph
    int N = sizeof(y[0]) / sizeof(y[0][0]);  // Frame size

    double complex f[M][M][M][K];


    double Ap = log(1.0/4.0);
    // printf("%f", Ap);

    for(int jj=0 ; jj<N ; jj++) {
        
        for (int m1 = 0; m1 < M; m1++) {
            for (int m2 = 0; m2 < M; m2++) {
                for(int m3=0 ; m3<M ; m3++) {
                    for(int k = 0 ; k<K ; k++) {
                        f[m1][m2][m3][k] = 0 + 0 * I;  // Complex zero
                    }
                }
            }
        }

        // for( int v=0 ; v<V ; v++) {
        //     printf("%d", x[jj][v]);
        //     printf("\n");
        // }
        
        // printf("\n");

        for (int k=0 ; k<K ; k++) {
            int a = 0;
            int ind[3] = {};
            for(int j=0 ; j<V ; j++) {
                if(F[k][j] == 1) {
                    ind[a] = j;
                    // printf("%d", ind[a]);
                    a++;
                }
            }

            // ind -- generated
            for(int m1=0 ; m1<M ; m1++) {
                for(int m2=0 ; m2<M ; m2++) {
                    for(int m3=0 ; m3<M ; m3++) {
                        d = y[k][jj] - ( 
                            (CB[ind[0]][k][m1] * h[k][ind[0]]) + 
                            (CB[ind[1]][k][m2] * h[k][ind[1]]) + 
                            (CB[ind[2]][k][m3] * h[k][ind[2]]) ) ;

                        f[m1][m2][m3][k] = -noise_pow * (
                            (creal(d) * creal(d)) + (cimag(d) * cimag(d)));
                        // printf("%f", f[m1][m2][m3][k]); // some values nan -------------------
                    }
                }
            }
        }

        double complex Igv[K][V][M];
        double complex Ivg[K][V][M];    

        for(int k=0 ; k<K ; k++) {
            for(int v=0 ; v<V ; v++) {
                for(int m=0 ; m<M ; m++) {
                    Igv[k][v][m] = 0 + 0 * I;  // Complex zero
                    Ivg[k][v][m] = Ap + 0 * I;

                }
            }
        }

        for(int iter=0 ; iter<Niter ; iter++) {
            // printf("%d", iter);
            for (int k=0 ; k<K ; k++) {
                int a = 0;
                int ind[3] = {};

                for(int j=0 ; j<V ; j++) {
                    if(F[k][j] == 1) {
                        ind[a] = j;
                        // printf("%d", ind[a]);
                        a++;
                    }
                }

                double complex sIgv[M * M];

                for(int m1=0 ; m1<M ; m1++) {
                    for(int m=0 ; m<(M*M) ; m++) {
                        sIgv[m] = 0 + 0 * I;
                    }

                    for(int m2=0 ; m2<M ; m2++) {
                        for(int m3 = 0 ; m3<M ; m3++) {
                            sIgv[(((m2+1)-1)*(M)) + m3] = f[m1][m2][m3][k] + 
                            Ivg[k][ind[1]][m2] + Ivg[k][ind[2]][m3] ;

                            Igv[k][ind[0]][m1] = log_sum_exp(sIgv, M*M);
                            // printf("%f", Igv[k][ind[0]][m1]); // some values nan -------------
                        }
                    }
                }

                // ============================== m2 ==================

                for(int m2=0 ; m2<M ; m2++) {
                    for(int m=0 ; m<(M*M) ; m++) {
                        sIgv[m] = 0 + 0 * I;
                    }

                    for(int m1=0 ; m1<M ; m1++) {
                        for(int m3 = 0 ; m3<M ; m3++) {
                            sIgv[(((m1+1)-1)*(M)) + m3] = f[m1][m2][m3][k] + 
                            Ivg[k][ind[0]][m1] + Ivg[k][ind[2]][m3] ;

                            Igv[k][ind[1]][m2] = log_sum_exp(sIgv, M*M);
                            // printf("%f", Igv[k][ind[0]][m1]); // all nan --------------------
                        }
                    }
                }

                // ========================= m3 ==============================

                for(int m3=0 ; m3<M ; m3++) {
                    for(int m=0 ; m<(M*M) ; m++) {
                        sIgv[m] = 0 + 0 * I;
                    }

                    for(int m1=0 ; m1<M ; m1++) {
                        for(int m2 = 0 ; m2<M ; m2++) {
                            sIgv[(((m1+1)-1)*(M)) + m2] = f[m1][m2][m3][k] + 
                            Ivg[k][ind[0]][m1] + Ivg[k][ind[1]][m2] ;

                            Igv[k][ind[2]][m3] = log_sum_exp(sIgv, M*M);
                            // printf("%f", Igv[k][ind[0]][m1]);
                        }
                    }
                }
            }

            // for(int k=0 ; k<K ; k++){
            //     for(int v=0 ; v<V ; v++){
            //         for(int m=0 ; m<M ; m++){
                        // printf("%f", Igv[k][v][m]);
            //         }
            //     }
            // }
        
            //=========================== Ivg ====================================================

            double complex s1;
            double complex s2;

            for (int k=0 ; k<V ; k++) {
                int a = 0;
                int ind[3] = {};

                for(int j=0 ; j<K ; j++) {
                    if(F[j][k] == 1) {
                        ind[a] = j;
                        // printf("%d", ind[a]);
                        a++;
                    }
                }

                // printf("\n");

                // s1 = logsumexp(Igv, ind[0], k);
                double complex sum_exp1 = 0.0 + 0.0 * I;  // Initialize to complex 0
                double complex sum_exp2 = 0.0 + 0.0 * I;  // Initialize to complex 0

                // Compute sum of exponentials over the slice Igv[ind][k][:]
                for (int m = 0; m < M; m++) {
                    sum_exp1 += cexp(Igv[ind[0]][k][m]);  // Use cexp() for complex exponential
                                    // printf("%f", sum_exp1);

                }

                s1 = clog(sum_exp1);
                // printf("%f", s1);

                                // Compute sum of exponentials over the slice Igv[ind][k][:]
                for (int m = 0; m < M; m++) {
                    sum_exp2 += cexp(Igv[ind[1]][k][m]);  // Use cexp() for complex exponential
                }

                s2 = clog(sum_exp2);
                // printf("%f", s2);

                for(int n=0 ; n<M ; n++) {
                    Ivg[ind[0]][k][n] = Igv[ind[1]][k][n] - s2;
                    Ivg[ind[1]][k][n] = Igv[ind[0]][k][n] - s1;
                    // printf("%.2f %+.2fi", creal(Ivg[ind[1]][k][n]), cimag(Ivg[ind[1]][k][n]));
                }
            }

            // Igv =============================================

            // for(int k=0; k<K ; k++) {
            //     for(int v=0 ; v<V ; v++) {
            //         printf("%.2f%+.2fi ", creal(Igv[k][v][1]), cimag(Igv[k][v][1]));
            //     }
            //     printf("\n");
            // }

            // global function ===============================

            // for(int m1=0 ; m1<M ; m1++){
            //     for(int m2=0 ; m2<M ; m2++){
            //         printf("%.2f%+.2fi ", creal(f[m1][m2][0][0]), cimag(f[m1][m2][0][0]));
            //     }
            //     printf("\n");
            // }

            // fading coeff. ===============================

            // for(int m1=0 ; m1<K ; m1++){
            //     for(int m2=0 ; m2<V ; m2++){
            //         printf("%.2f%+.2fi ", creal(h[m1][m2][0]), cimag(h[m1][m2][0]));
            //     }
            //     printf("\n");
            // }

            // for(int k=0 ; k<K ; k++){
            //     for(int n=0 ; n<N ; n++){
            //         printf("%.4f%+.4fi ", creal(y[k][n]), cimag(y[k][n]));
            //     }
            //     printf("\n");                
            // }
        }

        // ============================= LLR =====================================
        // =======================================================================
        
        double complex Q[M][V];

        for (int m = 0; m < M; m++) {
            for (int v = 0; v < V; v++) {
                Q[m][v] = 0 + 0 * I;  // Complex zero
            }
        }

        for(int k=0 ; k<V ; k++) {
            
            int a = 0;
            int ind[3] = {};

            for(int j=0 ; j<K ; j++) {
                if(F[j][k] == 1) {
                    ind[a] = j;
                    // printf("%d", ind[a]);
                    a++;
                }
                // printf("\t");
            }

            for(int m=0 ; m<M ; m++) {
                Q[m][k] = Ap + Igv[ind[0]][k][m] + Igv[ind[1]][k][m];
                // Ap+Igv(ind(1),k,m)+Igv(ind(2),k,m); 
                // printf("%f", Q[m][k]);
            }
        }

        // for(int m=0; m<M ; m++) {
        //         for(int k=0 ; k<K ; k++) {
        //             printf("%.2f%+.2fi ", creal(Q[m][k]), cimag(Q[m][k]));
        //         }
        //         printf("\n");
        //     }

        //     LLR_tmp = zeros(log2(M)*V, 1); % temp variable for parallel work
        
        double LLR_temp[12];

        for(int i=0 ; i<12 ; i++) {
                LLR_temp[i] = 0;
            }
        
        for(int k=0 ; k<V ; k++) {
            
            LLR_temp[2*k] = clog(
                (exp(Q[0][k]) + exp(Q[1][k])) 
                / 
                (exp(Q[2][k]) + exp(Q[3][k]))
                );
            LLR_temp[(2*k) + 1] = clog((exp(Q[0][k]) + exp(Q[2][k])) / (exp(Q[1][k]) + exp(Q[3][k])));            
        }

        double LLR[12][N];

        for (int k = 0; k < 12; k++) {
            for (int n = 0; n < N; n++) {
                LLR[k][n] = 0;  // Complex zero
            }
        }        

        for(int i=0 ; i<12 ; i++) {
            LLR[jj][i] = LLR_temp[i];
            // printf("%f", LLR[jj][i]);
            // printf("\n");
        }

        // printf("\n");
    }
}

// Helper function for log-sum-exp
double log_sum_exp(double *array, int size) {
    double max_val = array[0];
    for (int i = 1; i < size; i++) {
        if (array[i] > max_val) {
            max_val = array[i];
        }
    }

    double sum_exp = 0;
    for (int i = 0; i < size; i++) {
        sum_exp += exp(array[i] - max_val);
    }

    return max_val + log(sum_exp);
}

void scmaenc(int x[V][N], double complex CB[V][K][M], double complex h[K][V][N], double complex y[K][N]) {
    // Initialize y to zero
    for (int k = 0; k < K; k++) {
        for (int n = 0; n < N; n++) {
            y[k][n] = 0 + 0 * I;  // Complex zero
        }
    }

    // Main encoding loop
    for (int n = 0; n < N; n++) {
        for (int v = 0; v < V; v++) {
            int symbol_index = x[v][n];  // Symbol index (ranges from 0 to M-1)
            for (int k = 0; k < K; k++) {
                y[k][n] += CB[v][k][symbol_index] * h[k][v][n]; 
            }
        }
    }
}

int main() {



    double complex CB[V][K][M];
    double complex h[K][V][N];
    double complex s[K][N];

    for (int v = 0; v < V; v++) {
        for (int k = 0; k < K; k++) {
            for (int m = 0; m < M; m++) {
                CB[v][k][m] = CB_real[v][k][m] + CB_imag[v][k][m] * I; 
            }
        }
    }

    for (int k = 0; k < K; k++) {
        for (int v = 0; v < V; v++) {
            for (int n = 0; n < N; n++) {
                
                h[k][v][n] = h_real[k][v][n] + h_imag[k][v][n] * I; 

            }
        }
    }


            scmaenc(x, CB, h, s);

    for (int i = 0; i < K; i++) {
        for (int j = 0; j < N; j++) {
            printf("(%.6f + %.6fi)\t", creal(s[i][j]), cimag(s[i][j]));
        }
        printf("\n"); // New line after each row
    }


return 0;
}