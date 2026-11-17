#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
using namespace std;
int main(int argc,char**argv){
 assert(argc==3);ifstream in(argv[1]);ofstream out(argv[2]);uint32_t lines[6][5];for(auto &dir:lines)for(auto &l:dir)in>>l;
 array<int,46656> best;best.fill(100);array<uint32_t,46656> witness{};uint64_t admitted=0;
 for(uint32_t S=0;S<(1u<<25);++S){
  array<int,6> counts{};bool ok=true;
  for(int d=0;d<6;++d){for(int j=0;j<5;++j)counts[d]+=(S&lines[d][j])==lines[d][j];if(!counts[d]){ok=false;break;}}
  if(!ok)continue;++admitted;sort(counts.begin(),counts.end());int idx=0;for(int n:counts)idx=6*idx+n;
  int k=__builtin_popcount(S);if(k<best[idx]){best[idx]=k;witness[idx]=S;}
 }
 for(int i=0;i<46656;++i)if(best[i]<100)out<<i<<' '<<best[i]<<' '<<witness[i]<<'\n';
 cout<<"ALL_PLANAR_SUBSETS14_ENUMERATED "<<(1u<<25)<<" admitted "<<admitted<<endl;
}
