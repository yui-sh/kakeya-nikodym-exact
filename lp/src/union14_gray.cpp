#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;
int main(int argc,char**argv){
 assert(argc==3);ifstream in(argv[1]);ofstream out(argv[2]);vector<int> through[25];
 for(int l=0;l<30;++l){uint32_t mask;in>>mask;for(int p=0;p<25;++p)if(mask>>p&1)through[p].push_back(l);}
 for(auto &v:through)assert(v.size()==6);
 array<int,30> missing;missing.fill(5);array<int,6> counts{};array<int,46656> best;best.fill(100);array<uint32_t,46656>witness{};
 uint32_t S=0;int size=0,fixed_best=100;uint64_t admitted=0;int target[6]={1,1,2,2,3,3};
 for(uint32_t i=1;i<(1u<<25);++i){int p=__builtin_ctz(i);bool add=!(S>>p&1);S^=1u<<p;size+=add?1:-1;
  for(int l:through[p]){if(add){if(--missing[l]==0)++counts[l/5];}else{if(missing[l]++==0)--counts[l/5];}}
  bool ok=true;for(int d=0;d<6;++d)ok&=counts[d]>0;if(!ok)continue;++admitted;
  bool fixed=true;for(int d=0;d<6;++d)fixed&=counts[d]>=target[d];if(fixed)fixed_best=min(fixed_best,size);
  auto b=counts;sort(b.begin(),b.end());int idx=0;for(int n:b)idx=6*idx+n;if(size<best[idx]){best[idx]=size;witness[idx]=S;}
 }
 for(int i=0;i<46656;++i)if(best[i]<100)out<<i<<' '<<best[i]<<' '<<witness[i]<<'\n';
 cout<<"GRAY_PLANAR14_ENUMERATED "<<(1u<<25)<<" admitted "<<admitted<<" fixed_order_min "<<fixed_best<<endl;
}
