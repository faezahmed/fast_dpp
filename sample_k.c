/* Sampling from k-DPP inspired from
https://github.com/mehdidc/dpp

Determinantal point process sampling procedure based
on  (Fast Determinantal Point Process Sampling with
     Application to Clustering, Byungkon Kang, NIPS 2013)
*/

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

/*
   Recursive definition of determinate using expansion by minors.
*/
double Determinant(double **a,int n)
{
   int i,j,j1,j2;
   double det = 0;
   double **m = NULL;

   if (n < 1) { /* Error */

   } else if (n == 1) { 
      det = a[0][0];
   } else if (n == 2) {
      det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
   } else {
      det = 0;
      for (j1=0;j1<n;j1++) {
         m = malloc((n-1)*sizeof(double *));
         for (i=0;i<n-1;i++)
            m[i] = malloc((n-1)*sizeof(double));
         for (i=1;i<n;i++) {
            j2 = 0;
            for (j=0;j<n;j++) {
               if (j == j1)
                  continue;
               m[i-1][j2] = a[i][j];
               j2++;
            }
         }
         det += pow(-1.0,j1+2.0) * a[0][j1] * Determinant(m,n-1);
         for (i=0;i<n-1;i++)
            free(m[i]);
         free(m);
      }
   }
   //printf("\n %f ", det);
   return(det);
}

/*
   Find the cofactor matrix of a square matrix
*/
double ** CoFactor(double **a,int n)
{
   int i,j,ii,jj,i1,j1;
   double det;
   double **c;
   double **b;

   c = malloc((n-1)*sizeof(double *));
   b = malloc(n * sizeof(double *));

	

   for (i=0;i<n-1;i++)
     {
	 c[i] = (double*)malloc((n-1)*sizeof(double));
	}

   for (j=0;j<n;j++) {
	b[j] = (double*)malloc( sizeof(double) * n);
      for (i=0;i<n;i++) {

         /* Form the adjoint a_ij */
         i1 = 0;
         for (ii=0;ii<n;ii++) {
            if (ii == i)
               continue;
            j1 = 0;
            for (jj=0;jj<n;jj++) {
               if (jj == j)
                  continue;
               c[i1][j1] = a[ii][jj];
               j1++;
            }
            i1++;
         }

         /* Calculate the determinate */
         det = Determinant(c,n-1);
         /* Fill in the elements of the cofactor */
         b[j][i] = pow(-1.0,i+j+2.0) * det;
      }
   }
	//printMatrix (b,n);
   for (i=0;i<n-1;i++)
      free(c[i]);
   free(c);

return b;
}

/*
   Transpose of a square matrix, do it in place
*/
void inv_Transpose(double **a,int n,double detm)
{
   int i,j;
   double tmp;

   for (i=1;i<n;i++) {
      for (j=0;j<i;j++) {
         tmp = a[i][j]/detm;
         a[i][j] = a[j][i]/detm;
         a[j][i] = tmp;
      }
   }
}

/*
   A utility function to swap to integers
*/
void swap (int *a, int *b)
{
    int temp = *a;
    *a = *b;
    *b = temp;
}
 
/*
    A utility function to print integer array
*/
void printArray (int arr[], int n)
{
    int i;
    for (i = 0; i < n; i++)
        printf("%d ", arr[i]);
    printf("\n");
}

/*
    A utility function to print double array
*/
void printdArray (double arr[], int n)
{
    int i;
    for (i = 0; i < n; i++)
        printf("%lf ", arr[i]);
    printf("\n");
}

/*
 	A utility function to print an array
*/
void printMatrix (double **arr, int n)
{
    int i,j;
    for (i = 0; i<n; i++)
	{
		for(j=0;j<n;j++)
			printf("%lf ", arr[i][j]);
		printf("\n");
	}
}

/* 
	A function to generate a random permutation of arr[]
*/
void randomize ( int arr[], int n )
{
    // Use a different seed value so that we don't get same
    // result each time we run this program
    srand ( time(NULL) );
    int i;
    for (i = n-1; i > 0; i--)
    {
        // Pick a random index from 0 to i
        int j = rand() % (i+1);
 
        // Swap arr[i] with the element at random index
        swap(&arr[i], &arr[j]);
    }
}



/*
    Function to return first k elements of array 
*/
int * first_k (int arr[],int SIZE,int k)
{
	int i;
	int *initial;
	initial = (int *)malloc(k*sizeof(int));
	for(i=0;i<k;i++) 
	{ 
	  initial[i] = arr[i]; 
	} 
    return initial;
}

/*
    Function to subtract Set 2 from Set 1
*/
int * set_subtract(int set1[],int set2[],int n1,int n2)
{
	int k,i,j,flag,p;
	int *set3;
	set3 = (int *)malloc((n1-n2)*sizeof(int));

	k=0;
	for(i=0;i<n1;i++)
	{
		flag=1;
		for(j=0;j<n2;j++)
		{
			if(set1[i]==set2[j])
			{
				flag=0;
				break;
			}
		}
		if(flag==1)
		{
			set3[k]=set1[i];
			k++;
		}
	 }
	return set3;
}


/*
    Function to append one element to array end
*/
int * add_elem(int *pp,int v,int k)
{
	int *sset;
	int i;
	sset = (int *)malloc((k)*sizeof(int));
	for (i = 0; i < k-1; i++) 
	{
		sset[i]=pp[i];
	}
	
	sset[k-1]=v;

	return sset;

}


/*
    Function to remove particular element from array 
*/
int * remove_elem(int *pp,int u,int k)
{
	int *sset;
	int i;
	int j=0;
	sset = (int *)malloc((k-1)*sizeof(int));
	for (i = 0; i < k; i++) 
	{
		if (pp[i] != u) 
			sset[j++]=pp[i];
	}

	return sset;

}

/*
    Function to Removes all rows and columns from L matrix except those in set sset
*/
double ** find_submatrix(double **L,int *sset,int k)
{
	int i,j;
	double **L_Y;
	L_Y = malloc(k * sizeof(double *));
	for (i = 0; i < k; i++) 
	{
		(L_Y)[i] = (double*)malloc( sizeof(double) * k);
		for (j = 0; j < k; j++) 
		{
			L_Y[i][j]=L[sset[i]][sset[j]];
			//printf("(%d,%d) = %.4f %.4f\n", sset[i], sset[j], L_Y[i][j],L[sset[i]][sset[j]]);
		}
	}

	return L_Y;

}

/*
    Function to remove all elements from column v of L except those in sset
*/
double * find_column_subset(double **L,int *sset,int v,int k)
{
	int i;
	double *L_Y;
	L_Y = (double *)malloc(k*sizeof(double));
	for (i = 0; i < k; i++) 
	{
		L_Y[i]=L[sset[i]][v];
		//printf("(%d,%d) = %.4f %.4f\n", sset[i], sset[j], L_Y[i][j],L[sset[i]][sset[j]]);
	}

	return L_Y;
}

/*
    Function to read similarity kernel from text file
*/
double **read_simil(int n)
{

	double **array;
	array = malloc(n * sizeof(double *));

	int i,j;
	double x;
	FILE* pFile = fopen("similarity.txt", "r");
	if (pFile == NULL) 
	{
		printf("Error.");
		exit(1);
	}

	for (i = 0; i < n; i++) 
	{
		(array)[i] = (double*)malloc( sizeof(double) * n);
		for (j = 0; j < n; j++) 
		{	  
			fscanf(pFile, "%lf,", &x);		  
			array[i][j] = x;
			//printf("(%d,%d) = %.4f %.4f\n", i, j, x,array[i][j]);
		}
	}

	// Close files
	fclose(pFile);
	
	return array;


}
/*
	Function to calculate x'Lx
*/

double matr_product(double **L,double *x,int k)
{
	double result=0.0;
	int q,r;
		for(q=0;q<k;q++)
		{
			for(r=0;r<k;r++)
				result +=x[r]*L[r][q]*x[q];
		}
	return result;
}


int main()
{

//SIZE is total number of items and k is number of elements to sample
	int SIZE=100;
	int max_nb_iterations=1000;
	int k=10;

	int i;
	int rem=SIZE-k;
	int *p;
	int *arr;
	int *remset;
	int *sset;
	int u,v;
	double **L;
	double **L_Y;
	double **L_Y_inv;
	double detm,c_u,c_v;
	double *b_v;
	double *b_u;
	double rndn;
	double prob_v,prob_u,prob;

	b_v = (double *)malloc((k-1)*sizeof(double));
	b_u = (double *)malloc((k-1)*sizeof(double));
	p = (int *)malloc(k*sizeof(int));
	remset=(int *)malloc(rem*sizeof(int));
	arr= (int *)malloc(SIZE*sizeof(int));

	
	for (i = 0; i < SIZE; i++) 
		arr[i] = i;

	/*Select k random elements in p and remset is remaining set*/
	randomize (arr, SIZE);
	p=first_k(arr,SIZE,k);
	remset=set_subtract(arr,p,SIZE,k);

	/*L is similarity kernel*/
	L=read_simil(SIZE);
	//printMatrix (L,SIZE);

	for(i=0;i<max_nb_iterations;i++)
	{		
		u=p[rand()%k];
		v=remset[rand()%rem];

		sset=remove_elem(p,u,k);
		//printArray(sset, k-1);
		L_Y=find_submatrix(L,sset,k-1);
		//printMatrix (L_Y,k-1);

		/*Calculate determinant*/
		detm=Determinant(L_Y,k-1);

		if(detm!=0)
		{
			L_Y_inv=CoFactor(L_Y,k-1);
			inv_Transpose(L_Y_inv,k-1,detm);
			c_u=L[u][u];
			c_v=L[v][v];
			b_v=find_column_subset(L,sset,v,k-1);
			b_u=find_column_subset(L,sset,u,k-1);

			prob_v=c_v-matr_product(L_Y_inv,b_v,k-1);
			prob_u=c_u-matr_product(L_Y_inv,b_u,k-1);
			prob=prob_v/prob_u;
			
			if(prob>1)
				prob=1;
			rndn=(double)rand() / (double)RAND_MAX;
			//printf("\n%lf\n",prob);
			if(rndn<=prob)
			{
				p=add_elem(sset,v,k);
				printArray(p, k);
				remset=set_subtract(arr,p,SIZE,k);
			}
		
		}
		/*else
		{
			randomize (arr, SIZE);
			p=first_k(arr,SIZE,k);
		}*/

	}



	return(0);
}
