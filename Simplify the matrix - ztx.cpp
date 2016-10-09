/****************************************\
* Author : ztx
* Title  : Simplify the matrix
* ALG    :
* CMT    : 
* Time   :
\****************************************/

#include <cstdio>
#define Rep(i,l,r) for(i=(l);i<=(r);i++)
#define rep(i,l,r) for(i=(l);i< (r);i++)
#define Rev(i,r,l) for(i=(r);i>=(l);i--)
#define rev(i,r,l) for(i=(r);i> (l);i--)
typedef long long ll ;
typedef double lf ;
int CH , NEG ;
template <typename TP>inline void read(TP& ret) {
    ret = NEG = 0 ; while (CH=getchar() , CH<'!') ;
    if (CH == '-') NEG = true , CH = getchar() ;
    while (ret = ret*10+CH-'0' , CH=getchar() , CH>'!') ;
    if (NEG) ret = -ret ;
}
template <typename TP>inline void readc(TP& ret) {
    while (ret=getchar() , ret<'!') ;
    while (CH=getchar() , CH>'!') ;
}
template <typename TP>inline void reads(TP *ret) {
    ret[0]=0;while (CH=getchar() , CH<'!') ;
    while (ret[++ret[0]]=CH,CH=getchar(),CH>'!') ;
    ret[ret[0]+1]=0;
}

#include <algorithm>
#include <cstdlib>

using namespace std;

#define  maxn  100

int n, m;
int a[maxn][maxn];
int b[maxn][maxn];

inline void out() {
    int i, j;
    Rep (i,1,n) {
        printf("[");
        Rep (j,1,m) {
            int gcd = __gcd(a[i][j],b[i][j]);
            if (a[i][j] == 0) printf("%7d ",0);
            else if (b[i][j]/gcd == 1) printf("%7d ",a[i][j]/gcd);
            else {
                if (b[i][j]/gcd < 0)
                    printf("%3d/%-3d ",-a[i][j]/gcd,-b[i][j]/gcd);
                else printf("%3d/%-3d ",a[i][j]/gcd,b[i][j]/gcd);
            }
        }
        printf("   ]\n\n");
    }
    puts("");
}

int main() {
//    #define READ
    #ifdef  READ
        freopen(".in" ,"r",stdin ) ;
        freopen(".out","w",stdout) ;
    #endif
    int i, j, k, l;
    system("echo off");
    system("cls");
    system("color 0B");
    system("title SimplifyTheMatrix-ztx");
    puts("********************************************************************************\n\n");
    puts("                       Welcome to my program ! QwQ\n\n\n");
    puts("* Feature:                 matrix  => row echelon form\n");
    puts("                  row echelon form => reduced row echelon form\n");
    puts("*   Usage: input the N*M matrix, N is the row of it, M is ths column of it\n");
    puts("*          then input all the entries in the order of each row and column\n");
    puts("*  Notice: This program just to help you check for errors(maybe my program is wrong QAQ\n");
    puts("*          Do not run away from calculations :)\n");
    puts("* Version: 1.0\n");
    puts("* by FDU-JK1 ztx , Some rights reserved");
    puts("********************************************************************************");
    while (true) {
        //system("cls");
        printf("\n\nSolve N*M matrix, input the N and M:\n\n");
        scanf("%d%d", &n, &m);
        if (n<1 || m<1) {
            puts(">>> Error: N or M is illegal!\n\n");
            continue;
        }
        printf("Input the elements:\n\n");
        Rep (i,1,n) Rep (j,1,m) scanf("%d", &a[i][j]), b[i][j] = 1;
        k = -1;
        Rep (i,1,n) if (a[i][1]) { k = i ; break; }
        if (k < 0) {
            puts(">>> Error: The first column can not all be zero\n\n");
            continue ;
        }
system("cls");
        puts("Now solving:");
        puts("");
        out();
        puts("");
        int step = 0;
        int rownow = 1;
        Rep (i,1,n) {
            if (a[rownow][i] == 0) {
                k = -1;
                Rep (j,rownow+1,n) if (a[j][i] != 0) {
                    k = j; break;
                }
                if (k == -1) {
                    printf("The %d col don\'t need option\n\n", i);
                    continue;
                } else {
                    printf("Step %d: swap row %d & %d\n\n",++step,i,k);
                    Rep (j,1,m) swap(a[i][j],a[k][j]), swap(b[i][j],b[k][j]);
                    out();
                }
            }
            printf("Step %d: solve col %d\n\n",++step,i);
            int fenzi,fenmu;
            Rep (k,rownow+1,n) {
                if (k!=rownow && a[k][i]) {
                    fenzi = a[k][i]*b[rownow][i];
                    fenmu = a[rownow][i]*b[k][i];
                    int gcd=__gcd(fenzi,fenmu);
                    fenzi /= gcd;
                    fenzi = -fenzi;
                    fenmu /= gcd;
                    Rep (j,i,m) if (a[rownow][j]) {
                        int tmpzi,tmpmu,ttzi,ttmu;
                        tmpzi = a[rownow][j]*fenzi;
                        tmpmu = b[rownow][j]*fenmu;
                        ttzi = tmpmu*a[k][j]+tmpzi*b[k][j];
                        ttmu = tmpmu*b[k][j];
                        gcd = __gcd(ttzi,ttmu);
                        if (gcd) {
                            ttzi /= gcd;
                            ttmu /= gcd;
                        } else {
                            ttzi = 0;
                            ttmu = 1;
                        }
                        a[k][j] = ttzi;
                        b[k][j] = ttmu;
                    }
                }
            }
            out();
            rownow++;
            if (rownow > n) break;
        }
        printf("Step %d: Row echelon form >_<\n\n",++step);
        //out();
        Rep (i,1,n) {
            Rep (k,1,m)
                if (a[i][k] != 0) {
                    if ((a[i][k]*b[i][k]<0)||(a[i][k]*b[i][k]>0&&a[i][k]<0)) {
                        Rep (j,k,m) if (a[i][j]) a[i][j] = -a[i][j];
                    }
                    break;
                }
        }
        out();
        printf("Step %d: Reduced row echelon form >_<\n\n",++step);
        Rep (i,1,n) {
            Rep (k,1,m)
                if (a[i][k] != 0) {
                    Rep (j,k+1,m) if (a[i][j]) {
                        a[i][j] *= b[i][k];
                        b[i][j] *= a[i][k];
                    }
                    a[i][k] = 1, b[i][k] = 1;
                    break;
                }
        }
        //out();
        Rep (i,1,n) {
            Rep (l,1,m)
                if (a[i][l] != 0) {
                    Rep (k,1,n) if (a[k][l] != 0 && k != i) {
                        int fenzi = -a[k][l], fenmu = b[k][l];
                        int ttzi, ttmu, gcd;
                        a[k][l] = 0; b[k][l] = 1;
                        Rep (j,l+1,m) {
                            ttzi = a[i][j]*fenzi*b[k][j]+a[k][j]*b[i][j]*fenmu;
                            ttmu = b[i][j]*fenmu*b[k][j];
                            gcd = __gcd(ttzi,ttmu);
                            if (ttzi == 0) { a[k][j]=0;b[k][j]=1; }
                            else { a[k][j] = ttzi/gcd;b[k][j]=ttmu/gcd; }
                        }
                    }
                    break;
                }
        }
        out();
        //system("pause");
    }
    #ifdef  READ
        fclose(stdin) ; fclose(stdout) ;
    #else
        getchar() ; getchar() ;
    #endif
    return 0 ;
}