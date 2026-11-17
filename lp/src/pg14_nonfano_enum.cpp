#include <bit>
#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>
#include <functional>
int main(){
  int n,D,L;if(!(std::cin>>n>>D>>L)||n>63||D>n)return 1;
  std::vector<uint64_t> lines(L);for(auto &x:lines)std::cin>>x;
  std::vector<uint64_t> count(D+1),lo_mask(D+1),hi_mask(D+1);
  std::vector<int> lo(D+1,std::numeric_limits<int>::max()),hi(D+1,-1);
  std::function<void(int,int,uint64_t)> visit=[&](int start,int j,uint64_t mask){
    int triples=0;
    if(j>=3)for(auto line:lines){int c=std::popcount(line&mask);triples+=c*(c-1)*(c-2)/6-c*(c-1)*(c-2)*(c-3)/24;}
    ++count[j];
    if(triples<lo[j]){lo[j]=triples;lo_mask[j]=mask;}
    if(triples>hi[j]){hi[j]=triples;hi_mask[j]=mask;}
    if(j<D)for(int i=start;i<n;++i)visit(i+1,j+1,mask|(uint64_t(1)<<i));
  };
  visit(0,0,0);
  std::cout<<"[";
  for(int j=0;j<=D;++j){if(j)std::cout<<",";std::cout<<"{\"degree\":"<<j<<",\"subsets\":"<<count[j]<<",\"minimum\":"<<lo[j]<<",\"maximum\":"<<hi[j]<<",\"min_mask\":"<<lo_mask[j]<<",\"max_mask\":"<<hi_mask[j]<<"}";}
  std::cout<<"]\n";
}
