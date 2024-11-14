#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "err_est.h"
#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

int main(void) {

	int size_class1;
	int size_class2;
	int D;
	int D_total;
	double* gene_class1; //raw data of class 0
	double* gene_class2; //raw data of class 1
	float temp_float;
	int i, j;
	double** data_c1; //data vectors of class 0
	double** data_c2; //data vectors of class 1
	int index, index2;
	double *X1, *X2, *a, m;
	FILE *fp;
	
//    int D_special; 
	
	//-------------------------------------------------
	int p = 2; //number of features used in the classification
	//-------------------------------------------------
	
	dmem(a,p); // need to allocate memory for a, o.w. memcpy will fail

	//open the data file
	FILE *fpout; // output file containing LDA params and Error Estimations

	double resub, bresub;
	
	//-------------------------------------------------
	fp = fopen("/home/pwoods2/week15August2024grace/week15/grace/WT_ACK-TestOct2024grace.txt", "r");
	//-------------------------------------------------
	
	if (fp==NULL)
		return -1;
		
	

	//obtain the size information of two classes
	fscanf(fp, "%d, %d\r\n", &size_class1, &size_class2);

	//the number of total features
	fscanf(fp, "%d\r\n", &D);
	D_total = D;

    //generates file that is a portion of the original input file of size: D_special columns 
    //D_special = 200;
    //fpout = fopen("/Users/pwoods/Xcode-Files/lda-code-from-dr-ivanov/data_correct_format.txt", "w");
	//fprintf(fpout,"%i, %i\n%i\n", size_class1, size_class2, D_special);
	//fclose(fpout);
	////////////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////
	//allocate memory for gene data
	if ( (gene_class1=Malloc(double, size_class1*D_total))==NULL)
		exit(1);
	if ( (gene_class2=Malloc(double, size_class2*D_total))==NULL)
		exit(1);

	//obtain the gene data from data file
	for (i=0; i<size_class1; i++) {
		for (j=0; j<D_total; j++) {
			fscanf(fp, "%f\t", &temp_float);
			//			gene_class1[i*D_total+j]=log10(temp_float);
			gene_class1[i*D_total+j]=temp_float;
			//if ( j < D_special) fprintf(fpout, "%f\t", temp_float);
			//			printf("%i,%f\n", i*D_total+j, gene_class1[i*D_total+j]);
		}
		fscanf(fp, "\r\n");
		//fprintf(fpout, "\n");
	}

	//
	for (i=0; i<size_class2; i++) {
		for (j=0; j<D_total; j++) {
			fscanf(fp, "%f\t", &temp_float);
			//			gene_class2[i*D_total+j]=log10(temp_float);
			gene_class2[i*D_total+j]=temp_float;
			//if ( j < D_special) fprintf(fpout, "%f\t", temp_float);
			//			printf("%i,%f\n", i*D_total+j, gene_class2[i*D_total+j]);
		}
		fscanf(fp, "\r\n");
		//fprintf(fpout, "\n");
	}

	//close file
	fclose(fp);
    //fclose(fpout);

	//	//////////////////////////////////////////////////////
	//allocate memory for data vector
	if ( (data_c1 = Malloc(double*, size_class1))==NULL)
		exit(1);
	if ( (data_c2 = Malloc(double*, size_class2))==NULL)
		exit(1);

	// set the pointer
	for (i = 0; i < size_class1; i++) {
		data_c1[i] = &gene_class1[i*D_total];

		//		printf("%i,%f\n", i, *data_c1[i]);
	}
	for (i = 0; i < size_class2; i++) {
		data_c2[i] = &gene_class2[i*D_total];
		//		printf("%i,%f\n", i, *data_c2[i]);

	}

	//	//allocate memory for X1,X2
	if ( (X1=Malloc(double, size_class1*p))==NULL) {
		printf("error");
		exit(1);
	}
	if ( (X2=Malloc(double, size_class2*p))==NULL) {
		printf("error");
		exit(1);
	}
	//-------------------------------------------------
	int ii, jj, pp, qq;
	fpout = fopen("/home/pwoods2/week15August2024grace/week15/grace/2-Feature-WT-ACK-Test-Output-Grace.txt", "w");
	fprintf(fpout, "gene index\ta_n\tm\tresub\tbresub\n");
	//-------------------------------------------------
	// the following for loops are for exhaustive search
	for (ii=0; ii<D_total; ii++) 
		for (jj=ii+1; jj<D_total; jj++)
//			for (pp = jj+1; pp<D_total; pp++)
//				for (qq = pp+1; qq<D_total; qq++) 
			{
//					int order[] = { ii };
                  int order[] = { ii, jj };
//  			    int order[] = { ii, jj, pp };
//					int order[] = { ii, jj, pp, qq};

					for (index = 0; index<size_class1; index++) {
						for (index2 = 0; index2<p; index2++) {
							X1[index*p+index2] = *(data_c1[index]+order[index2]);
						//	printf("X1:%i,%f\n", index*p+index2, X1[index*p+index2]);
						}
					}
								printf("\n\n");

					//	// set up X2	
					for (index = 0; index<size_class2; index++) {
						for (index2 = 0; index2<p; index2++) {
							X2[index*p+index2] = *(data_c2[index]+order[index2]);
						//	printf("X2:%i,%f\n", index*p+index2, X2[index*p+index2]);
						}
					}

					// LDA and bolstered EE	


					lda(X1, X2, size_class1,size_class2, p, a, &m);

					//			//			printf("LDA parameters:\n");
					//			//			printf("a = [");
					//			//			for (i=0; i<p; i++)
					//			//				printf(" %.4f", a[i]);
					//			//			printf(" ]\nm = %.4f\n\n", m);
					//
					//			//			 Compute various error estimates 
					resub = resub_lda(X1, X2, size_class1, size_class2, p, a, m);
					//							printf("Resubstitution estimate = %.4f\n\n", resub);
					//
					bresub = bresub_lda(X1, X2, size_class1, size_class2, p, a, m);
					//							printf("Bolstered resubstitution estimate = %.4f\n\n", bresub);
					//			  printf("X1[0]=%f\n",X1[0]);
					//			  
					//			  printf("X2[0]=%f\n",X2[0]);
					//			  printf("X1[1]=%f\n",X1[1]);
					//			  printf("X2[1]=%f\n\n",X2[1]);

//					fprintf(fpout,"[%i]\t[%.4f]\t%.4f\t%.4f\t%.4f\n",ii,a[0],m,resub,bresub);
//					fprintf(fpout, "[%i,%i]\t[%.4f  %.4f]\t%.4f\t%.4f\t%.4f\n", ii, jj,a[0], a[1], m, resub, bresub);
					if (bresub<0.10){
//  				fprintf(fpout,"[%i,%i,%i]\t[%.4f  %.4f  %.4f]\t%.4f\t%.4f\t%.4f\n", ii, jj, pp, a[0], a[1], a[2], m, resub, bresub);
                  fprintf(fpout, "[%i,%i]\t[%.4f  %.4f]\t%.4f\t%.4f\t%.4f\n", ii, jj,a[0], a[1], m, resub, bresub);
//                  fprintf(fpout,"[%i]\t[%.4f]\t%.4f\t%.4f\t%.4f\n",ii,a[0],m,resub,bresub);
					}
					//			if (bresub<=0.25) {
					//				fprintf(
					//						fpout,
					//						"[%i,%i,%i,%i]\t[%.4f  %.4f  %.4f %.4f]\t%.4f\t%.4f\t%.4f\n",
					//						ii, jj, pp, qq, a[0], a[1], a[2], a[3], m, resub,
					//						bresub);
					//			}
				
				}
            


fclose(fpout);
free(a);
free(X1);
free(X2);
free(gene_class1);
free(gene_class2);
free(data_c1);
free(data_c2);

system("PAUSE");

return 0;

}
