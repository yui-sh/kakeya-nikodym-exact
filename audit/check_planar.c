/* AUDIT (fresh code): enumerate all 2^25 subsets S of AG(2,5).
   Outputs (to stdout):
     SEC k r u c1..cr | t1..t5 | e1..e5      distinct tuples over all choices of r full lines in distinct directions
                                             (e = sorted sections of the unmarked direction, only when r==5)
     FT  a1..a6 min                           min |S| over S whose sorted full-line-count profile dominates (a1..a6) is
                                             computed in Python from the raw profile minima printed as
     PROF p1..p6 min                          min |S| with sorted full-line counts exactly (p1..p6)
   Independent of lp/src: own point/line construction. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define Q 5
#define NP 25
#define NL 30
static uint32_t lmask[NL]; static int ldir[NL];
typedef struct { unsigned char k,r,u,c[7],t[6],e[6]; } tup; /* 26 bytes */
#define HS (1u<<20)
static tup tab[HS]; static unsigned char used[HS]; static long ntup=0;
static uint64_t hsh(const tup*x){uint64_t h=1469598103934665603ULL;const unsigned char*p=(const unsigned char*)x;for(size_t i=0;i<sizeof(tup);i++){h^=p[i];h*=1099511628211ULL;}return h;}
static void ins(const tup*x){uint64_t h=hsh(x)&(HS-1);while(used[h]){if(!memcmp(&tab[h],x,sizeof(tup)))return;h=(h+1)&(HS-1);}used[h]=1;tab[h]=*x;ntup++;}
static int profmin[6][6][6][6][6][6];
int main(void){
  /* lines: direction d<5: y = d*x + c ; d==5: x = c.  point index 5x+y */
  for(int d=0;d<=Q;d++) for(int c=0;c<Q;c++){ int li=d*Q+c; ldir[li]=d; uint32_t m=0;
    for(int x=0;x<Q;x++){ int y; if(d<Q){ y=(d*x+c)%Q; m|=1u<<(Q*x+y);} else { m|=1u<<(Q*c+x); } }
    lmask[li]=m; }
  /* sanity: each line 5 points, two lines of different direction meet in exactly one point */
  for(int a=0;a<NL;a++){ if(__builtin_popcount(lmask[a])!=5){fprintf(stderr,"bad line\n");return 1;}
    for(int b=0;b<NL;b++) if(a!=b){ int c=__builtin_popcount(lmask[a]&lmask[b]); if(ldir[a]==ldir[b]?c!=0:c!=1){fprintf(stderr,"bad incidence %d %d\n",a,b);return 1;} } }
  memset(profmin,0x3f,sizeof profmin);
  for(uint32_t S=0;S<(1u<<NP);S++){
    int sec[NL]; int full[6][5]; int nfull[6]={0};
    tup x; memset(&x,0,sizeof x); x.k=__builtin_popcount(S);
    for(int l=0;l<NL;l++){ sec[l]=__builtin_popcount(S&lmask[l]); if(sec[l]) x.t[sec[l]]++; if(sec[l]==Q) full[ldir[l]][nfull[ldir[l]]++]=l; }
    { int p[6]; for(int d=0;d<6;d++) p[d]=nfull[d]; /* sort ascending */
      for(int i=0;i<6;i++) for(int j=i+1;j<6;j++) if(p[j]<p[i]){int tt=p[i];p[i]=p[j];p[j]=tt;}
      if(x.k<profmin[p[0]][p[1]][p[2]][p[3]][p[4]][p[5]]) profmin[p[0]][p[1]][p[2]][p[3]][p[4]][p[5]]=x.k; }
    /* r = 0 */
    x.r=0; x.u=0; ins(&x);
    int fd[6], nd=0; for(int d=0;d<6;d++) if(nfull[d]) fd[nd++]=d;
    for(int sub=1; sub<(1<<nd); sub++){
      int ds[6], r=0; for(int i=0;i<nd;i++) if(sub>>i&1) ds[r++]=fd[i];
      int idx[6]={0};
      while(1){
        int cnt[NP]={0}; uint32_t un=0;
        for(int i=0;i<r;i++){ uint32_t lm=lmask[full[ds[i]][idx[i]]]; un|=lm; for(int p=0;p<NP;p++) if(lm>>p&1) cnt[p]++; }
        tup y=x; y.r=r; y.u=__builtin_popcount(un); memset(y.c,0,sizeof y.c); memset(y.e,0,sizeof y.e);
        for(int p=0;p<NP;p++) if(cnt[p]) y.c[cnt[p]]++;
        if(r==5){ int ud=-1; for(int d=0;d<6;d++){ int m=0; for(int i=0;i<r;i++) if(ds[i]==d) m=1; if(!m) ud=d; }
          int e[5]; for(int c=0;c<5;c++) e[c]=sec[ud*Q+c];
          for(int i=0;i<5;i++) for(int j=i+1;j<5;j++) if(e[j]<e[i]){int tt=e[i];e[i]=e[j];e[j]=tt;}
          for(int i=0;i<5;i++) y.e[i]=e[i]; }
        ins(&y);
        int i=0; while(i<r){ idx[i]++; if(idx[i]<nfull[ds[i]]) break; idx[i]=0; i++; } if(i==r) break;
      }
    }
  }
  for(uint32_t h=0;h<HS;h++) if(used[h]){ tup*x=&tab[h]; printf("SEC %d %d %d",x->k,x->r,x->u); for(int i=1;i<=x->r;i++) printf(" %d",x->c[i]); printf(" |"); for(int j=1;j<=Q;j++) printf(" %d",x->t[j]); if(x->r==5){ printf(" |"); for(int i=0;i<5;i++) printf(" %d",x->e[i]); } printf("\n"); }
  for(int a=0;a<6;a++)for(int b=a;b<6;b++)for(int c=b;c<6;c++)for(int d=c;d<6;d++)for(int e=d;e<6;e++)for(int f=e;f<6;f++) if(profmin[a][b][c][d][e][f]<0x3f3f3f3f) printf("PROF %d %d %d %d %d %d %d\n",a,b,c,d,e,f,profmin[a][b][c][d][e][f]);
  fprintf(stderr,"distinct SEC tuples: %ld\n",ntup);
  return 0; }
