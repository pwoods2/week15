#include <stdio.h>    
#include <stdlib.h>   
#include <math.h>     
#include <string.h>   
#include "err_est.h"  

#define Malloc(type,n) (type *)malloc((n)*sizeof(type)) // Macro for typed malloc

int main(int argc, char *argv[]) {
    // Check if exact number of arguments is provided
    if (argc != 4) {
        printf("Usage: %s <input_file> <output_file> <number_of_features>\n", argv[0]);
        return 1;
    }

    // Get arguments from command line
    char *input_file = argv[1];   // Path to input file
    char *output_file = argv[2];  // Path to output file
    int p = atoi(argv[3]);        // Number of features (1 to 4)

    // Ensure valid feature count
    if (p < 1 || p > 4) {
        printf("Feature count must be between 1 and 4.\n");
        return 1;
    }

    // Variable declarations
    int size_class1, size_class2, D, D_total; // Sizes of each class and feature dimension
    double *gene_class1, *gene_class2;        // Flattened arrays of gene expression data
    float temp_float;                         // Temporary float for reading from file
    int i, j;                                  // Loop counters
    double **data_c1, **data_c2;              // 2D views into the gene data arrays
    int index, index2;                        // Index counters
    double *X1, *X2, *a, m;                   // X1/X2: class data in current feature space, a: LDA weights, m: mean threshold
    FILE *fp, *fpout;                         // File pointers
    double resub, bresub;                     // Error estimates

    a = Malloc(double, p); // Allocate memory for LDA coefficients

    // Open and read input file
    fp = fopen(input_file, "r");
    if (fp == NULL) {
        perror("Failed to open input file");
        return 1;
    }

    // Read number of samples and features
    fscanf(fp, "%d, %d\r\n%d\r\n", &size_class1, &size_class2, &D);
    D_total = D;

    // Allocate memory for gene expression data
    if ((gene_class1 = Malloc(double, size_class1 * D_total)) == NULL) exit(1);
    if ((gene_class2 = Malloc(double, size_class2 * D_total)) == NULL) exit(1);

    // Read class 1 data
    for (i = 0; i < size_class1; i++) {
        for (j = 0; j < D_total; j++) {
            fscanf(fp, "%f\t", &temp_float);
            gene_class1[i * D_total + j] = temp_float;
        }
        fscanf(fp, "\r\n");
    }

    // Read class 2 data
    for (i = 0; i < size_class2; i++) {
        for (j = 0; j < D_total; j++) {
            fscanf(fp, "%f\t", &temp_float);
            gene_class2[i * D_total + j] = temp_float;
        }
        fscanf(fp, "\r\n");
    }

    fclose(fp); // Close input file

    // Allocate memory for row pointers to each sample (2D view of flattened data)
    if ((data_c1 = Malloc(double*, size_class1)) == NULL) exit(1);
    if ((data_c2 = Malloc(double*, size_class2)) == NULL) exit(1);

    // Set row pointers
    for (i = 0; i < size_class1; i++) data_c1[i] = &gene_class1[i * D_total];
    for (i = 0; i < size_class2; i++) data_c2[i] = &gene_class2[i * D_total];

    // Allocate memory for reduced feature data for current p-combination
    if ((X1 = Malloc(double, size_class1 * p)) == NULL) exit(1);
    if ((X2 = Malloc(double, size_class2 * p)) == NULL) exit(1);

    // Open output file
    fpout = fopen(output_file, "w");
    if (fpout == NULL) {
        perror("Failed to open output file");
        return 1;
    }

    // Write output file header
    fprintf(fpout, "gene index\ta_n\tm\tresub\tbresub\n");

    // Iterate through all combinations of p features (p = 1–4)
    int ii, jj, pp, qq;
    for (ii = 0; ii < D_total; ii++) {
        if (p == 1) {
            int order[] = { ii };

            // Extract feature values for both classes based on selected feature
            for (index = 0; index < size_class1; index++)
                for (index2 = 0; index2 < p; index2++)
                    X1[index * p + index2] = *(data_c1[index] + order[index2]);

            for (index = 0; index < size_class2; index++)
                for (index2 = 0; index2 < p; index2++)
                    X2[index * p + index2] = *(data_c2[index] + order[index2]);

            // Run LDA and calculate error estimates
            lda(X1, X2, size_class1, size_class2, p, a, &m);
            resub = resub_lda(X1, X2, size_class1, size_class2, p, a, m);
            bresub = bresub_lda(X1, X2, size_class1, size_class2, p, a, m);

            // Write result
            fprintf(fpout, "[%i]\t[%.4f]\t%.4f\t%.4f\t%.4f\n", ii, a[0], m, resub, bresub);
        }

        for (jj = ii + 1; jj < D_total && p >= 2; jj++) {
            if (p == 2) {
                int order[] = { ii, jj };

                // Extract features
                for (index = 0; index < size_class1; index++)
                    for (index2 = 0; index2 < p; index2++)
                        X1[index * p + index2] = *(data_c1[index] + order[index2]);

                for (index = 0; index < size_class2; index++)
                    for (index2 = 0; index2 < p; index2++)
                        X2[index * p + index2] = *(data_c2[index] + order[index2]);

                // LDA and error calculations
                lda(X1, X2, size_class1, size_class2, p, a, &m);
                resub = resub_lda(X1, X2, size_class1, size_class2, p, a, m);
                bresub = bresub_lda(X1, X2, size_class1, size_class2, p, a, m);

                // Write result
                fprintf(fpout, "[%i,%i]\t[%.4f,%.4f]\t%.4f\t%.4f\t%.4f\n", ii, jj, a[0], a[1], m, resub, bresub);
            }

            for (pp = jj + 1; pp < D_total && p >= 3; pp++) {
                if (p == 3) {
                    int order[] = { ii, jj, pp };

                    for (index = 0; index < size_class1; index++)
                        for (index2 = 0; index2 < p; index2++)
                            X1[index * p + index2] = *(data_c1[index] + order[index2]);

                    for (index = 0; index < size_class2; index++)
                        for (index2 = 0; index2 < p; index2++)
                            X2[index * p + index2] = *(data_c2[index] + order[index2]);                lda(X1, X2, size_class1, size_class2, p, a, &m);
                resub = resub_lda(X1, X2, size_class1, size_class2, p, a, m);
                bresub = bresub_lda(X1, X2, size_class1, size_class2, p, a, m);

                fprintf(fpout, "[%i,%i,%i]\t[%.4f,%.4f,%.4f]\t%.4f\t%.4f\t%.4f\n",
                        ii, jj, pp, a[0], a[1], a[2], m, resub, bresub);
            }

            for (qq = pp + 1; qq < D_total && p == 4; qq++) {
                int order[] = { ii, jj, pp, qq };

                for (index = 0; index < size_class1; index++)
                    for (index2 = 0; index2 < p; index2++)
                        X1[index * p + index2] = *(data_c1[index] + order[index2]);

                for (index = 0; index < size_class2; index++)
                    for (index2 = 0; index2 < p; index2++)
                        X2[index * p + index2] = *(data_c2[index] + order[index2]);

                lda(X1, X2, size_class1, size_class2, p, a, &m);
                resub = resub_lda(X1, X2, size_class1, size_class2, p, a, m);
                bresub = bresub_lda(X1, X2, size_class1, size_class2, p, a, m);

                fprintf(fpout, "[%i,%i,%i,%i]\t[%.4f,%.4f,%.4f,%.4f]\t%.4f\t%.4f\t%.4f\n",
                        ii, jj, pp, qq, a[0], a[1], a[2], a[3], m, resub, bresub);
            }
        }
    }
}

// Free all allocated memory
free(gene_class1);
free(gene_class2);
free(data_c1);
free(data_c2);
free(X1);
free(X2);
free(a);

fclose(fpout);

return 0;
}
