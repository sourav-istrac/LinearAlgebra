#include <stddef.h>
#include<stdio.h>
#include<math.h>
#include<string.h>

int rref(int m,int n,double A[m][n],double R[m][n],double O[n][m]);
int matrixMult(int m,int n,double A[m][n],int p,int q,double B[p][q],double C[m][m]);
void readMatrix(int opt,int rows,int cols,double M[rows][cols],FILE *fp);
void cpyMatrix(int rows,int cols,double A[rows][cols],double B[rows][cols]);

int main(int arg,char *argc[]){
    int rows,cols,opt;
    FILE *fp;
    if(arg>1)
        opt=1;
    else
        opt=0;
    if(opt==1){
        char inp[100];
        fp=fopen(argc[1],"r");
        if(fp!=NULL){
            char inp[100];
            fgets(inp,100,fp);
            sscanf(inp,"%d %d",&rows,&cols);
        }else {
            opt = 0;   
        }
    }
    if(opt!=1){
        printf("Please Enter rows and cols respectively : ");
        scanf("%d %d",&rows,&cols);
        fp=NULL;
    }
    double M[rows][cols];
    readMatrix(opt,rows,cols,M,fp);
    if(fp!=NULL)
        fclose(fp);
    double R[rows][cols],O[rows][rows];
    rref(rows,cols,M,R,O);
    printf("Before rref : \n");
    int i,j;
    for(i=0;i<rows;i++){
        for(j=0;j<cols;j++){
            printf("%15.3lf ",M[i][j]);
        }
        printf("\n");
    }
    printf("After rref : \n");
    for(i=0;i<rows;i++){
        for(j=0;j<cols;j++){
            printf("%15.3lf ",R[i][j]);
        }
        printf("\n");
    }
    printf("Total Operations : \n");
    for(i=0;i<rows;i++){
        for(j=0;j<rows;j++){
            printf("%15.3lf ",O[i][j]);
        }
        printf("\n");
    }
}
void readMatrix(int opt,int rows,int cols,double M[rows][cols],FILE *fp){
    int i,j;
    if(opt==1 && fp!=NULL){
        for(i=0;i<rows;i++){
            for(j=0;j<cols;j++){
                fscanf(fp,"%lf",&M[i][j]);
            }
        }
    }else{
        printf("Please enter the %d elements of %dX%d Matrix:\n",rows*cols,rows,cols);
        for(i=0;i<rows;i++){
            for(j=0;j<cols;j++){
                scanf("%lf",&M[i][j]);
            }
        }

    }
}
int matrixMult(int m,int n,double A[m][n],int p,int q,double B[p][q],double C[m][q]){
    if(n!=p){
        printf("Matrix Multiplication not possible!!\n");
        return 1;
    }
    int i,j,k;
    for(i=0;i<m;i++){
        for(j=0;j<q;j++){
            C[i][j]=0;
            for(k=0;k<n;k++){
                C[i][j]+=A[i][k]*B[k][j];
            }
        }
    }
    return 0;
}
// copy Matrix B to A
// rows,cols,dest. M,Sou. M
void cpyMatrix(int rows,int cols,double A[rows][cols],double B[rows][cols]){
    int i,j;
    for(i=0;i<rows;i++){
        for(j=0;j<cols;j++){
            A[i][j]=B[i][j];
        }
    }
}

int rref(int m,int n,double A[m][n],double R[m][n],double O[m][m]){
    int p=0; // piviot
    int maxP;
    if(m>n)
        maxP = n;
    else
        maxP = m;
    int i,j,k;
    double opr[m][m],oprTemp[m][m],temp[m][n];
    for(j=0;j<m;j++){
        for(k=0;k<m;k++){
            O[j][k]=0;
            if(j==k)
                O[j][k]=1;
        }
    }
    for(j=0;j<m;j++){
        for(k=0;k<n;k++){
            R[j][k]=A[j][k];
        }
    }
    while(p<maxP){
        if(R[p][p]==0){
            for(i=p+1;i<m;i++){
                if(R[i][p]!=0){
                    for(j=0;j<m;j++){
                        for(k=0;k<m;k++){
                            opr[j][k]=0;
                            if(j==k)
                                opr[j][k]=1;
                        }
                    }
                    opr[i][i]=0;
                    opr[i][p]=1;
                    opr[p][p]=0;
                    opr[p][i]=1;
                    matrixMult(m,m,opr,m,n,R,temp);
                    cpyMatrix(m,n,R,temp);
                    matrixMult(m,m,opr,m,m,O,oprTemp);
                    cpyMatrix(m,m,O,oprTemp);
                    printf("Operations changing rows : \n");
                    for(i=0;i<m;i++){
                        for(j=0;j<m;j++){
                            printf("%15.3lf ",opr[i][j]);
                        }
                        printf("\n");
                    }
                    printf("--------------------\n");
                    for(i=0;i<m;i++){
                        for(j=0;j<n;j++){
                            printf("%15.3lf ",R[i][j]);
                        }
                        printf("\n");
                    }
                }
                
            }

        }
        if(R[p][p]!=0){
            for(j=0;j<m;j++){
                for(k=0;k<m;k++){
                    opr[j][k]=0;
                    if(j==k)
                        opr[j][k]=1;
                }
            }
            opr[p][p]=1/R[p][p];
            for(i=0;i<m;i++){
                if(i!=p && R[i][p]!=0){
                    opr[i][p]=(-1*R[i][p]/R[p][p]);
                }
            }
            matrixMult(m,m,opr,m,n,R,temp);
            cpyMatrix(m,n,R,temp);
            matrixMult(m,m,opr,m,m,O,oprTemp);
            cpyMatrix(m,m,O,oprTemp);
            printf("Operations reducing : \n");
            for(i=0;i<m;i++){
                for(j=0;j<n;j++){
                    printf("%15.3lf ",R[i][j]);
                }
                printf("\n");
            }

        }
        p++;
    }
    return 0;
}