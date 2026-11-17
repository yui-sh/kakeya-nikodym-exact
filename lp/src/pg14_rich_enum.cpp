#include <bit>
#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>
#include <functional>
int main(){
 int n,D,L;if(!(std::cin>>n>>D>>L)||n>63||D>n)return 1;
 std::vector<uint64_t> lines(L);for(auto &x:lines)std::cin>>x;
 std::vector<uint64_t> count(D+1),capcount(D+1),lm(D+1),hm(D+1),clm(D+1),chm(D+1);
 std::vector<int> lo(D+1,999),hi(D+1,-1),clo(D+1,999),chi(D+1,-1);
 std::function<void(int,int,uint64_t)> visit=[&](int start,int j,uint64_t mask){
  int value=0,maxb=0;for(auto line:lines){int b=std::popcount(line&mask);value+=b>=3;maxb=std::max(maxb,b);}
  ++count[j];if(value<lo[j]){lo[j]=value;lm[j]=mask;}if(value>hi[j]){hi[j]=value;hm[j]=mask;}
  if(maxb<=4){++capcount[j];if(value<clo[j]){clo[j]=value;clm[j]=mask;}if(value>chi[j]){chi[j]=value;chm[j]=mask;}}
  if(j<D)for(int i=start;i<n;++i)visit(i+1,j+1,mask|(uint64_t(1)<<i));
 };visit(0,0,0);std::cout<<"[";
 for(int j=0;j<=D;++j){if(j)std::cout<<",";std::cout<<"{\"degree\":"<<j<<",\"subsets\":"<<count[j]<<",\"minimum\":"<<lo[j]<<",\"maximum\":"<<hi[j]<<",\"min_mask\":"<<lm[j]<<",\"max_mask\":"<<hm[j]<<",\"cap_subsets\":"<<capcount[j]<<",\"cap_minimum\":"<<clo[j]<<",\"cap_maximum\":"<<chi[j]<<",\"cap_min_mask\":"<<clm[j]<<",\"cap_max_mask\":"<<chm[j]<<"}";}
 std::cout<<"]\n";
}
