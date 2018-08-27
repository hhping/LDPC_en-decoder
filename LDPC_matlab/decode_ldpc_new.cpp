
/*LDPC Decoder*/

#include "mex.h"
#include "matrix.h"			
#include "math.h"
#include <stdlib.h>			
#include "decodeutil_new.h"

#define INF 1000  

// for fixed point
#define  wordLen_din_frag       2
#define  wordLen_tanh_frag      6
#define  wordLen_tanh_mul       6
#define  wordLen_atanh_frag     4
#define  wordLen_atanh_real     3


double atanh(double x);


void decode(double max_iter, double *vhat, int mrows, int ncols, double *iter_out, double *gamma_n,
			double *check_node_ones,double max_check_degree, double BIGVALUE_COLS,
			double *variable_node_ones,double max_variable_degree,double BIGVALUE_ROWS)

{//function braces
	
	int i=0, j=0;
	
	double **sg_array,**sa_array, **bitmessage_temp ;
	double **H, **check_node_ones_matrix,**variable_node_ones_matrix;
	
	sg_array		= matrix(0,mrows-1,0,ncols-1);	
	sa_array		= matrix(0,mrows-1,0,ncols-1);	
	H				= matrix(0,mrows-1,0,ncols-1);	

	check_node_ones_matrix=matrix(0,mrows-1,0,max_check_degree-1);	
	variable_node_ones_matrix=matrix(0,max_variable_degree-1,0,ncols-1);



		for (i=0;i<max_check_degree;i++)
		{
			for(j=0;j<mrows;j++)
			{
				check_node_ones_matrix[j][i]=*(check_node_ones++); 
			
			}
		}


		for(i=0;i<ncols;i++)
		{
			for(j=0;j<max_variable_degree;j++)
			{
				variable_node_ones_matrix[j][i]=*(variable_node_ones++);
			}
			
		}
	
		
		double element=0,col_index=0,row_index=0;
		int temp_row_index=0;

		//initializing the matrices
		for ( i=0;i<ncols;i++)
		{
            gamma_n[i] =gamma_n[i]/pow(2,wordLen_din_frag);   
            
			for(int j=0;j<max_variable_degree;j++)
			{
				row_index=variable_node_ones_matrix[j][i];
				if(row_index==BIGVALUE_ROWS)
					break;

				temp_row_index=(int)row_index;
                
                sg_array[temp_row_index][i]=gamma_n[i];
				sa_array[temp_row_index][i]=0;
				H[temp_row_index][i]=1;
			}

		}


	bitmessage_temp = matrix(0,mrows-1,0,ncols-1);
	double *sum_of_b;
	sum_of_b=vector(0,ncols-1);
	double *vhat_temp;
	vhat_temp =vector(0,ncols-1);
    double tt1, tt2;
    
	for (int iteration=0;iteration<max_iter;iteration++)  
	{	//main iteration loop
	
		//bit-to-check messages
		for ( i=0;i<mrows;i++)
		{
			for (j=0;j<max_check_degree;j++)
			{

				col_index=check_node_ones_matrix[i][j];

				if (col_index==BIGVALUE_COLS)					
				{
					break;
				}

                tt1 = tanh (   (     (sg_array[i][(int)col_index]) - sa_array[i][(int)col_index]     ) / 2  );
                
                tt1 = floor(tt1 * pow(2,wordLen_tanh_frag));                  
                tt1 = tt1/ pow(2,wordLen_tanh_frag);           
                            
                bitmessage_temp[i][(int)col_index] =tt1;
			}

		}



		for (int u=0;u<mrows;u++)
		{//across mrows
			double temp=1;
			
			for (int v =0;v<max_check_degree;v++)	
			{//across columns
				
				element=check_node_ones_matrix[u][v];  

				if (element==BIGVALUE_COLS)					
				{
					break;
				}
				
				for (int w=0;w<max_check_degree;w++)
					{//accross columns again
						
						col_index=check_node_ones_matrix[u][w];

						if(check_node_ones_matrix[u][w]!=element && check_node_ones_matrix[u][w]!=BIGVALUE_COLS)
						{//second if
							temp=temp* bitmessage_temp[u][(int)col_index]; 
						}//second if
                        


					}//accross columns again
                    

                    temp = floor(temp * pow(2,wordLen_tanh_mul));                   
                    temp = temp/ pow(2,wordLen_tanh_mul);     
                                  
                    tt2=2*atanh(temp);
                   

                    tt2 = floor(tt2 * pow(2,wordLen_atanh_frag));                    
                    tt2 = tt2/ pow(2,wordLen_atanh_frag);
                    
                    if (tt2>pow(2,wordLen_atanh_real))
                    {
                        tt2 = pow(2,wordLen_atanh_real);
                    }
                    
                    sa_array[u][(int)element]=tt2;
					temp=1;
				
			
			}//across columns
		
		}//accross mrows    

		for ( i=0;i<ncols;i++)
		{
			double temp=0;
			for(int j=0;j<max_variable_degree;j++)
			{	
				row_index=variable_node_ones_matrix[j][i];
				if(row_index==BIGVALUE_ROWS)
					break;


				temp=temp+sa_array[(int)row_index][i];
			}
		
			sum_of_b[i]=temp;

		}

		for ( i=0;i<ncols;i++)
		{
			for(int j=0;j<max_variable_degree;j++)
			{
				
				row_index=variable_node_ones_matrix[j][i];
				if(row_index==BIGVALUE_ROWS)
					break;

				sg_array[(int)row_index][i]=sum_of_b[i]+gamma_n[i];
			}

			vhat_temp[i]=sum_of_b[i]+gamma_n[i];
		

			if (vhat_temp[i]>=0)
                *(vhat+i)=0;
			else
                *(vhat+i)=1;
		}
		
	    int parity=0,cumsum=0;

		for ( j=0;j<mrows;j++)
		{
			
			for (i=0;i<ncols;i++)
			{
				*(vhat_temp+i)=H[j][i]* *(vhat+i);
			}
			
			parity=0;
			for(i=0;i<ncols;i++)
			{
				parity=(parity + (int)(*(vhat_temp+i)) )%2;
				 
			}
		
			cumsum=parity+cumsum;
			
			if (cumsum==1) 
			{	break;	}
		}

		if (cumsum==0) 
		{	
			*iter_out=iteration+1;
			return; 
		}

	}	//main iteration loop
		
	*iter_out=iteration;

}//function braces



double atanh(double x) {
  double epsilon;	
  epsilon = pow(10,-16);
  
  if(x>(1-epsilon)) return INF;
  if(x<(-1+epsilon)) return -INF;
  return 0.5*log((1+x)/(1-x));
}


void mexFunction( int nlhs, mxArray *plhs[], 
				  int nrhs, const mxArray*prhs[] )
{
	double *vhat, *iter_out,*gamma_n, *check_node_ones, *variable_node_ones; /*pointer variables for input Matrices*/
	double max_iter,max_check_degree,max_variable_degree, BIGVALUE_COLS,BIGVALUE_ROWS; 
	int mrows,ncols;


	gamma_n  = mxGetPr(prhs[1]); 
	check_node_ones=mxGetPr(prhs[2]);  
	max_check_degree=mxGetScalar(prhs[3]); 
	BIGVALUE_COLS=mxGetScalar(prhs[4]);

	variable_node_ones=mxGetPr(prhs[5]); 
	max_variable_degree=mxGetScalar(prhs[6]); 
	BIGVALUE_ROWS=mxGetScalar(prhs[7]);


	max_iter = mxGetScalar(prhs[0]); 
	
	  mrows = mxGetScalar(prhs[8]);
	  ncols = mxGetScalar(prhs[9]); 

	 plhs[0] = mxCreateDoubleMatrix(1,ncols, mxREAL); 
	 vhat = mxGetPr(plhs[0]);	
	 plhs[1] = mxCreateDoubleScalar(0);
	 iter_out = mxGetPr(plhs[1]);

	 decode(max_iter,vhat,mrows,ncols,iter_out,gamma_n,check_node_ones,max_check_degree,BIGVALUE_COLS,variable_node_ones,max_variable_degree,BIGVALUE_ROWS);


}