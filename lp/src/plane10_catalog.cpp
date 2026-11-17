// Materialize every distinct joint dm multiset, with one geometric degree witness.
// Input is the same finite jobs as plane10_enum. Costs are ignored.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>
using namespace std;
struct Hash {size_t operator()(const array<uint16_t,6>&a)const{size_t h=0;for(auto x:a)h=h*1000003+x;return h;}};
int D,k,r;vector<int> h,mask,ds,pow6,fullgid;array<int,6> state;
unordered_map<array<uint16_t,6>,uint64_t,Hash> columns;long long nodes,leaves;
void dfs(int p){
 ++nodes;
 if(p==k){
  ++leaves;array<uint16_t,6> key;key.fill(65535);
  for(int j=0;j<r;j++){if(fullgid[state[j]]<0)return;key[j]=fullgid[state[j]];}
  sort(key.begin(),key.end());uint64_t packed=0;for(int j=0;j<k;j++)packed|=uint64_t(ds[j])<<(3*j);
  columns.emplace(key,packed);return;
 }
 int lo=max(1,__builtin_popcount((unsigned)mask[p]));if(p && mask[p]==mask[p-1])lo=max(lo,ds[p-1]);
 for(int d=lo;d<=D;d++)if(h[d-1]){
  ds[p]=d;--h[d-1];bool ok=true;
  for(int j=0;j<r;j++)if(mask[p]>>j&1){state[j]+=pow6[d-1];int count=0,v=state[j];while(v){count+=v%6;v/=6;}if(count==5 && fullgid[state[j]]<0)ok=false;}
  if(ok)dfs(p+1);
  for(int j=0;j<r;j++)if(mask[p]>>j&1)state[j]-=pow6[d-1];
  ++h[d-1];
 }
}
int main(int argc,char**argv){
 if(argc!=3)return 2;
 ifstream in(argv[1]);ofstream out(argv[2],ios::binary);int jobs,mode;in>>jobs>>mode;uint64_t total=0;
 for(int job=0;job<jobs;job++){
  uint32_t tid,geo;long long base;in>>tid>>geo>>D>>k>>r>>base;if(k>21 || D>7)return 3;
  h.resize(D);mask.resize(k);ds.resize(k);for(auto&v:h)in>>v;for(auto&v:mask)in>>v;
  pow6.assign(D+1,1);for(int j=1;j<=D;j++)pow6[j]=6*pow6[j-1];fullgid.assign(pow6[D],-1);
  int n;in>>n;for(int z=0;z<n;z++){int g,code=0,d;long long cost;in>>g>>cost;for(int j=0;j<5;j++){in>>d;code+=pow6[d-1];}fullgid[code]=g;}
  state.fill(0);columns.clear();nodes=leaves=0;dfs(0);
  for(auto& [key,packed]:columns){out.write((char*)&tid,4);out.write((char*)&geo,4);out.write((char*)key.data(),12);out.write((char*)&packed,8);}
  total+=columns.size();cerr<<"type "<<tid<<" geometry "<<geo<<" nodes "<<nodes<<" leaves "<<leaves<<" distinct "<<columns.size()<<'\n';
 }
 cerr<<"CATALOG_COMPLETE jobs "<<jobs<<" distinct "<<total<<'\n';
}
