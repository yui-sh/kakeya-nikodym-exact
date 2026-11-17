// Exact finite DFS. Pruning uses integer lower bounds over every compatible dm.
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <vector>
#include <climits>
using namespace std;
constexpr long long INF=LLONG_MAX/8;
int D,k,r,tid,geo; bool first,done; long long base,best,nodes,leaves;
vector<int> h,mask,ds,bd,bg,pow6; array<int,6> state;
vector<long long> lower;vector<int> fullgid;vector<vector<int>> suffix;
void subhist(int j,int code,const vector<int>&cnt,long long c){
 if(j==D){lower[code]=min(lower[code],c);return;}
 for(int a=0;a<=cnt[j];a++)subhist(j+1,code+a*pow6[j],cnt,c);
}
void dfs(int p){
 ++nodes;if(done)return;
 long long bound=base;
 for(int j=0;j<r;j++){if(lower[state[j]]==INF)return;bound+=lower[state[j]];}
 if(bound>=best)return;
 if(p==k){
  ++leaves;best=bound;bd=ds;bg.clear();for(int j=0;j<r;j++){if(fullgid[state[j]]<0)abort();bg.push_back(fullgid[state[j]]);}
  if(first)done=true;
  return;
 }
 int lo=max(1,__builtin_popcount((unsigned)mask[p]));
 if(p && mask[p]==mask[p-1])lo=max(lo,ds[p-1]);
 for(int d=lo;d<=D;d++)if(h[d-1]){
  ds[p]=d;--h[d-1];bool ok=true;int cap=0;
  for(int j=D;j>=1;j--){cap+=h[j-1];if(suffix[p+1][j]>cap){ok=false;break;}}
  if(ok){
   for(int j=0;j<r;j++)if(mask[p]>>j&1)state[j]+=pow6[d-1];
   dfs(p+1);
   for(int j=0;j<r;j++)if(mask[p]>>j&1)state[j]-=pow6[d-1];
  }
  ++h[d-1];if(done)return;
 }
}
int main(int argc,char**argv){
 if(argc!=3)return 2;
 ifstream in(argv[1]);ofstream out(argv[2]);int jobs,ff;in>>jobs>>ff;first=ff;
 long long total=0,allLeaves=0;
 for(int job=0;job<jobs;job++){
  in>>tid>>geo>>D>>k>>r>>base;h.resize(D);mask.resize(k);ds.resize(k);for(int&v:h)in>>v;for(int&v:mask)in>>v;
  pow6.assign(D+1,1);for(int j=1;j<=D;j++)pow6[j]=pow6[j-1]*6;
  lower.assign(pow6[D],INF);fullgid.assign(pow6[D],-1);
  int n;in>>n;for(int z=0;z<n;z++){int g,d,code=0;long long c;in>>g>>c;vector<int> cnt(D,0);for(int j=0;j<5;j++){in>>d;cnt[d-1]++;code+=pow6[d-1];}fullgid[code]=g;subhist(0,0,cnt,c);}
  suffix.assign(k+1,vector<int>(D+1,0));for(int p=k-1;p>=0;p--){suffix[p]=suffix[p+1];for(int j=1;j<=D;j++)if(__builtin_popcount((unsigned)mask[p])>=j)suffix[p][j]++;}
  best=INF;done=false;nodes=leaves=0;state.fill(0);dfs(0);total+=nodes;allLeaves+=leaves;
  if(best!=INF){out<<tid<<' '<<geo<<' '<<best;for(int v:bg)out<<' '<<v;for(int v:bd)out<<' '<<v;out<<'\n';}
 }
 cerr<<"jobs "<<jobs<<" first "<<first<<" nodes "<<total<<" improving_leaves "<<allLeaves<<'\n';
}
