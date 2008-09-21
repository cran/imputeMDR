#include <stdlib.h>
#include <math.h>
#include <R.h>


void make_table_em(int* c_loc, int* c_num, int lc, int* dataset, int rd, int cd, double *table, 
					int* vec1, int* vec2, int* vec3, int len_vec1)
{
	int i, l;
	int loc;
	int nct=2*pow(4,lc);
	double* current_t= (double*) calloc(nct,sizeof(double));
	double* update_t= (double*) calloc(nct,sizeof(double));
	double* scale= (double*) calloc(len_vec1*2,sizeof(double));
	int iteration=0;
	int vec2_start=0;
	double* c_prob= (double*) calloc(2,sizeof(double));
	double update_total=1000;

	for (i=0; i<rd; i++)
	{
		loc = dataset[i];
		for (l=0; l<lc; l++) 
                    loc += c_num[l] * dataset[c_loc[l]+i];
		table[loc]++;
	}

	memset(scale, 0, len_vec1*2*sizeof(double));
    memcpy(update_t, table, nct*sizeof(double));
    memcpy(current_t, table, nct*sizeof(double));
	for (i=0; i<2; i++) c_prob[i]=0;

	// iteration & comparison
	while (iteration<30 && fabs(update_total)>0.001)
	{
		iteration++;
		vec2_start=0;
	    update_total=0;
		memset(scale, 0, len_vec1*2*sizeof(double));

	    memcpy(update_t, table, nct*sizeof(double));

		// update
		for (i=0; i<len_vec1; i++)
		{	
			for (l=0; l<vec3[i]; l++) {
				scale[i*2]=scale[i*2]+current_t[vec2[vec2_start+l]*2];			 //controls
				scale[i*2+1]=scale[i*2+1]+current_t[vec2[vec2_start+l]*2+1]; //cases
			}

				c_prob[0]=(double) 1/vec3[i];
				c_prob[1]=(double) 1/vec3[i];

			for (l=0; l<vec3[i]; l++) {
               // M-step : obtain new cell probabilities
				if (scale[i*2]>0)
						c_prob[0]=current_t[vec2[vec2_start+l]*2] / scale[i*2];
				
				if (scale[i*2+1]>0)
						c_prob[1]=current_t[vec2[vec2_start+l]*2+1] / scale[i*2+1];
				
				//E-step : impute missing values by expectation
				update_t[vec2[vec2_start+l]*2]= (double)(update_t[vec2[vec2_start+l]*2] 
						+ table[vec1[i]*2] *  c_prob[0]);//controls
				update_t[vec2[vec2_start+l]*2+1]= (double)(update_t[vec2[vec2_start+l]*2+1] 
						+ table[vec1[i]*2+1] * c_prob[1] );//cases
			}
			vec2_start=vec2_start+vec3[i];
		}

		for (i=0; i<nct; i++) 
				update_total = update_total + pow((update_t[i]-current_t[i]),2);

		memcpy(current_t, update_t, nct*sizeof(double));

	}

	for (i=0; i<nct; i++) 
            table[i]=update_t[i];
	for (i=0;i<len_vec1 ;i++ ){
		   table[vec1[i]*2]=0;
		   table[vec1[i]*2+1]=0;
	}

	free(current_t);
	free(update_t);
	free(scale);
}


/* missing: none, complete, missing category */
void make_table_none(int* c_loc, int* c_num, int lc, 
                int* dataset, int rd, int cd, double* table)
{
	int i, j;
	int loc;
    /* initial location: case or control */
    /* sample_index must start from 0 */
	/* suppose there is no missing data */
	for (i = 0; i < rd; i++) {   // for each individual data
                // case or control value sets the start position
		loc = dataset[i];		 
		for (j = 0; j < lc; j++) 
			loc += c_num[j] * dataset[c_loc[j] + i];
		table[loc]++; 
	}
}


/* use available complete data*/
void make_table_available(int* c_loc, int* c_num, int lc, 
                int* dataset, int rd, int cd, double* table, int* newcase, int* newcontrol)
{
    int i, j;
    int loc;
    int geno_val;
   *newcase = *newcontrol = 0;
	
    for (i = 0; i < rd; i++) {
        /* initial location: case or control */
        /* sample_index must start from 0 */
        loc = dataset[i];
        for (j = 0; j < lc; j++) {
            geno_val = dataset[c_loc[j] + i];
            if (geno_val == 3)  /* if missing */
				goto missing_val;
            loc += c_num[j] * geno_val; 
        }

		if (dataset[i] ==1 ){
            (*newcase)++;
			
     	} else{
           (*newcontrol)++;
		}

        table[loc]++;
		missing_val:
			;
    }
}


void err_rate(int* comb, int* rcomb, int* ccomb,
              int* train, int* rt, int* ct, int* test, int *rtest,
              double* threshold,
              double* err_train, double* err_test, int* na_method,
			  int* vec1, int* vec2, int* vec3, int* n_vec1)
{

    int level = 3;
	int len_vec1=*n_vec1;
	/* levels is 4 for missing category and EM method */
    if (*na_method == 1||*na_method == 3)
        level = 4;

    int j, k, l;
    /* for calculating error rate */
    double err1, err2;
    int nct = 2 * pow(level, *rcomb); // depending on data size
    double ratio;

    int* c_num = (int*) calloc(*rcomb, sizeof(int));
    int* c_loc = (int*) calloc(*rcomb, sizeof(int));
    /* counting table */
    //int* t_train = (int*) calloc(nct, sizeof(int));
    //int* t_test = (int*) calloc(nct, sizeof(int));

    /* counting table in case of missing */
    double* missing_t_train = (double*) calloc(nct, sizeof(double));
    double* missing_t_test = (double*) calloc(nct, sizeof(double));
	int case_train=0,control_train=0, case_test=0, control_test=0;
	int train_total=0, test_total=0;		//denominator in error computation

    /* base location of table */
    c_num[0] = 2;
    if (*rcomb > 1)
            for (l = 1; l < *rcomb; l++) 
                    c_num[l] = c_num[l - 1] * level;			// depending on data size
    
    for (j = 0; j < *ccomb; j++) {
			/* initialization */
            for (k = 0; k < nct; k++)
                    missing_t_train[k] = missing_t_test[k] = 0;


            /* count the cell from train dataset */ 
            for (l = 0; l < *rcomb; l++)
                    c_loc[l] = (comb[*rcomb * j + l] - 1) * *rt;			// get position of each SNP
            if( (*na_method==0) || (*na_method==1) ){					// complete, missing category
                    make_table_none(c_loc, c_num, *rcomb, train, *rt, *ct, missing_t_train);
            } else if (*na_method==2)	{
                    make_table_available(c_loc, c_num, *rcomb, train, *rt, *ct, missing_t_train,&case_train,&control_train);
            } else if (*na_method==3){
				make_table_em(c_loc, c_num, *rcomb, train, *rt, *ct, missing_t_train, vec1, vec2, vec3, len_vec1);
            }

            
            /* count the cell from test dataset */ 
            for (l = 0; l < *rcomb; l++)
                    c_loc[l] = (comb[*rcomb * j + l] - 1) * *rtest;
            if(*na_method==0||*na_method==1	){					 // complete, missing category
                    make_table_none(c_loc, c_num, *rcomb, test, *rtest, *ct, missing_t_test);
            } else if (*na_method==2)	{
                    make_table_available(c_loc, c_num, *rcomb, test, *rtest, *ct, missing_t_test,&case_test,&control_test);
            } else if (*na_method==3){
                    make_table_em(c_loc, c_num, *rcomb, test, *rtest, *ct, missing_t_test, vec1, vec2, vec3, len_vec1);
            }


            /* calculate the train and test error rate */
            err1 = err2 = 0;
			if (*na_method==2)	{		//available MDR uses variable threshold for each SNP combination model
				train_total= case_train+ control_train;
				test_total= case_test+ control_test;
			} else {
				train_total=*rt;
				test_total=*rtest;
			}

            for (k = 0; k < nct / 2; k++) {
                    /* calculate ratio, declaration is somewhat ambigous */
                    if (missing_t_train[k * 2] == 0)
                            ratio = *threshold;
                    else
                            ratio = (double) missing_t_train[k * 2 + 1] / missing_t_train[k * 2]; 
                    if (ratio >= *threshold) {
                            err1 += missing_t_train[k * 2];
                            err2 += missing_t_test[k * 2];
                    } else {
                            err1 += missing_t_train[k * 2 + 1];
                            err2 += missing_t_test[k * 2 + 1];
                    }
            }

			err_train[j] = (double) err1 / train_total;
            err_test[j] = (double) err2 / test_total; 


}

    free(c_num);
    free(c_loc);
}

