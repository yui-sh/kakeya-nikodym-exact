/* AUDIT (fresh): catalogue of all planes of AG(2,5) with r=5 or 6 inside lines (one per chosen direction) and any set
   of extra points: for each (lines, S) output key (k,r,u,t1..t5,e1..e5) and the sorted multiset of per-point line-incidence
   masks (which lines through each point of S). Distinct (key, mask-multiset) pairs only. No orbit reduction. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define Q 5
#define NP 25
#define NL 30
static uint32_t lmask[NL]; static int ldir[NL];
typedef struct { unsigned char k,r,u,t[6],e[6]; unsigned char mm[26]; } rec; /* mm = sorted masks (as bytes, 6-bit) */
#define HS (1u<<22)
static rec *tab; static unsigned char *used; static long nrec=0;
static uint64_t hsh(const rec*x){uint64_t h=1469598103934665603ULL;const unsigned char*p=(const unsigned char*)x;for(size_t i=0;i<sizeof(rec);i++){h^=p[i];h*=1099511628211ULL;}return h;}
static void ins(const rec*x){uint64_t h=hsh(x)&(HS-1);while(used[h]){if(!memcmp(&tab[h],x,sizeof(rec)))return;h=(h+1)&(HS-1);}used[h]=1;tab[h]=*x;nrec++;}
static int cmpb(const void*a,const void*b){return (int)*(const unsigned char*)a-(int)*(const unsigned char*)b;}
int main(void){
  tab=calloc(HS,sizeof(rec)); used=calloc(HS,1);
  for(int d=0;d<=Q;d++) for(int c=0;c<Q;c++){ int li=d*Q+c; ldir[li]=d; uint32_t m=0;
    for(int x=0;x<Q;x++){ int y; if(d<Q){ y=(d*x+c)%Q; m|=1u<<(Q*x+y);} else { m|=1u<<(Q*c+x); } } lmask[li]=m; }
  for(int r=5;r<=6;r++){
    for(int dsub=0;dsub<64;dsub++){ if(__builtin_popcount(dsub)!=r) continue;
      int ds[6],nd=0; for(int d=0;d<6;d++) if(dsub>>d&1) ds[nd++]=d;
      int ud=-1; for(int d=0;d<6;d++) if(!(dsub>>d&1)) ud=d;
      int idx[6]={0};
      while(1){
        int ls[6]; uint32_t U=0; for(int i=0;i<r;i++){ ls[i]=ds[i]*Q+idx[i]; U|=lmask[ls[i]]; }
        int out[NP],no=0; for(int p=0;p<NP;p++) if(!(U>>p&1)) out[no++]=p;
        for(uint32_t bits=0; bits<(1u<<no); bits++){
          uint32_t S=U; for(int i=0;i<no;i++) if(bits>>i&1) S|=1u<<out[i];
          rec x; memset(&x,0,sizeof x); x.k=__builtin_popcount(S); x.r=r; x.u=__builtin_popcount(U);
          for(int l=0;l<NL;l++){ int c=__builtin_popcount(S&lmask[l]); if(c) x.t[c]++; }
          if(r==5){ int e[5]; for(int c=0;c<5;c++) e[c]=__builtin_popcount(S&lmask[ud*Q+c]);
            for(int i=0;i<5;i++) for(int j=i+1;j<5;j++) if(e[j]<e[i]){int tt=e[i];e[i]=e[j];e[j]=tt;} for(int i=0;i<5;i++) x.e[i]=e[i]; }
          int n=0; for(int p=0;p<NP;p++) if(S>>p&1){ unsigned char m=0; for(int i=0;i<r;i++) if(lmask[ls[i]]>>p&1) m|=1<<i; x.mm[n++]=m; }
          qsort(x.mm,n,1,cmpb);
          ins(&x);
        }
        int i=0; while(i<r){ idx[i]++; if(idx[i]<Q) break; idx[i]=0; i++; } if(i==r) break;
      }
    }
  }
  for(uint32_t h=0;h<HS;h++) if(used[h]){ rec*x=&tab[h]; printf("%d %d %d |",x->k,x->r,x->u); for(int j=1;j<=Q;j++) printf(" %d",x->t[j]); printf(" |"); if(x->r==5) for(int i=0;i<5;i++) printf(" %d",x->e[i]); printf(" |"); for(int i=0;i<x->k;i++) printf(" %d",x->mm[i]); printf("\n"); }
  fprintf(stderr,"distinct (key, mask-multiset): %ld\n",nrec); return 0; }
