/* Session VI: enumerate all subsets S of AG(2,q) (q=5: 2^25) and record the distinct tuples
   (k, r, u, c_1..c_r, t_1..t_q) where k=|S|, t_j = # lines of AG(2,q) meeting S in exactly j points,
   and (u; c) is the crossing profile of some choice of r full lines (lines contained in S) with pairwise
   distinct directions (u = |union|, c_i = # points on exactly i of the chosen lines). r ranges over 0..#directions
   having a full line.  Output: one tuple per line to stdout.  Usage: sections_enum q */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define MAXQ 5
int RMIN=5;
int q, np, nd, nl;
uint32_t lmask[64]; int ldir[64];
/* hash set of tuples */
typedef struct { int k, r, u; signed char c[8]; signed char t[8]; signed char e[8][8]; int ne; } tup;
#define HS (1<<22)
tup *tab; unsigned char *used; long ntup=0;
static uint64_t hsh(const tup*x){ uint64_t h=1469598103934665603ULL; const unsigned char*p=(const unsigned char*)x; for(size_t i=0;i<sizeof(tup);i++){h^=p[i];h*=1099511628211ULL;} return h; }
void ins(tup*x){ uint64_t h=hsh(x)&(HS-1); while(used[h]){ if(!memcmp(&tab[h],x,sizeof(tup))) return; h=(h+1)&(HS-1);} used[h]=1; tab[h]=*x; ntup++; }
int main(int argc,char**argv){
  q=atoi(argv[1]); if(argc>2) RMIN=atoi(argv[2]); np=q*q; nd=q+1; nl=q*(q+1);
  tab=calloc(HS,sizeof(tup)); used=calloc(HS,1);
  { int rp, rl; if(scanf("%d %d",&rp,&rl)!=2||rp!=np||rl!=nl){fprintf(stderr,"bad header\n");return 1;}
    for(int l=0;l<nl;l++){ unsigned m; if(scanf("%d %u",&ldir[l],&m)!=2) return 1; lmask[l]=m; } }
  if(0){
  /* lines: direction index d: d<q -> (1,d), d==q -> (0,1) */
  int li=0;
  for(int d=0;d<nd;d++){ int a=(d<q)?1:0, b=(d<q)?d:1;
    for(int s=0;s<q;s++){ /* lines of direction d: start point (0,s) if a==1 else (s,0) */
      int x0=(a==1)?0:s, y0=(a==1)?s:0; uint32_t m=0;
      for(int t=0;t<q;t++){ int x=(x0+t*a)%q, y=(y0+t*b)%q; m|=1u<<(q*x+y); }
      lmask[li]=m; ldir[li]=d; li++; } } }
  uint32_t total=1u<<np;
  for(uint32_t S=0;S<total;S++){
    tup x; memset(&x,0,sizeof x); x.k=__builtin_popcount(S);
    int full[64], nf=0; int fulldir[8]; int ndirs=0; int dirfull[8][8], dircnt[8]; memset(dircnt,0,sizeof dircnt);
    for(int l=0;l<nl;l++){ int c=__builtin_popcount(S&lmask[l]); if(c) x.t[c]++; if(c==q){ dirfull[ldir[l]][dircnt[ldir[l]]++]=l; } }
    for(int d=0;d<nd;d++) if(dircnt[d]) fulldir[ndirs++]=d;
    /* r=0 */
    x.r=0; x.u=0; ins(&x);
    /* r>=1: subsets of fulldir, one line per direction */
    for(int sub=1; sub<(1<<ndirs); sub++){
      int r=__builtin_popcount(sub); int ds[8], m=0; for(int i=0;i<ndirs;i++) if(sub>>i&1) ds[m++]=fulldir[i];
      /* iterate over product of choices */
      int idx[8]; memset(idx,0,sizeof idx);
      while(1){
        int cnt[32]; memset(cnt,0,sizeof cnt);
        for(int i=0;i<r;i++){ uint32_t lm=lmask[dirfull[ds[i]][idx[i]]]; for(int p=0;p<np;p++) if(lm>>p&1) cnt[p]++; }
        tup y=x; y.r=r; y.u=0; memset(y.c,0,sizeof y.c);
        for(int p=0;p<np;p++) if(cnt[p]){ y.u++; y.c[cnt[p]-1]++; }
        if(r>=RMIN){ /* profiles of all unmarked directions, each sorted, then the list sorted lexicographically */
          int used_d[8]={0}; for(int i=0;i<r;i++) used_d[ds[i]]=1; y.ne=0;
          for(int d=0;d<nd;d++) if(!used_d[d]){ int ne=0; signed char *pr=y.e[y.ne];
            for(int l=0;l<nl;l++) if(ldir[l]==d) pr[ne++]=__builtin_popcount(S&lmask[l]);
            for(int a=0;a<ne;a++) for(int b=a+1;b<ne;b++) if(pr[b]<pr[a]){ signed char tmp=pr[a]; pr[a]=pr[b]; pr[b]=tmp; }
            y.ne++; }
          for(int a=0;a<y.ne;a++) for(int b=a+1;b<y.ne;b++) if(memcmp(y.e[b],y.e[a],8)<0){ signed char tmp[8]; memcpy(tmp,y.e[a],8); memcpy(y.e[a],y.e[b],8); memcpy(y.e[b],tmp,8); }
        }
        ins(&y);
        int i=0; while(i<r){ idx[i]++; if(idx[i]<dircnt[ds[i]]) break; idx[i]=0; i++; } if(i==r) break;
      }
    }
  }
  for(long h=0;h<HS;h++) if(used[h]){ tup*x=&tab[h]; printf("%d %d %d",x->k,x->r,x->u); for(int i=0;i<x->r;i++) printf(" %d",x->c[i]); printf(" |"); for(int j=1;j<=q;j++) printf(" %d",x->t[j]); for(int a=0;a<x->ne;a++){ printf(" |"); for(int j=0;j<q;j++) printf(" %d",x->e[a][j]); } printf("\n"); }
  fprintf(stderr,"distinct tuples: %ld\n",ntup);
  return 0; }
