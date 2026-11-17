// Exact integer dynamic programming for histogram pencil completion/pricing.
// No floating arithmetic or solver. Each repeated item occurrence gets a choice;
// ordered occurrences only duplicate columns, never omit a multiset.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <unordered_map>
#include <vector>
using namespace std;
using I = int64_t;
const I INF = (I(1)<<60);
struct Choice { int id; array<int,7> h{}; I cost; };
struct Val { I cost; int choice; };
struct Node { int item=-1, next=-1; array<int,7> lo{}, hi{}; unordered_map<uint64_t,Val> memo; };
int D;
vector<vector<Choice>> groups;
vector<Node> nodes;
uint64_t pack(const array<int,7>& h) { uint64_t z=0; for(int j=0;j<D;j++){assert(h[j]>=0 && h[j]<64);z|=uint64_t(h[j])<<(6*j);} return z; }
Val solve(int s, const array<int,7>& rem) {
    if(s==0) { for(int j=0;j<D;j++) if(rem[j]) return {INF,-1}; return {0,-1}; }
    auto &nd=nodes[s];
    for(int j=0;j<D;j++) if(rem[j]<nd.lo[j] || rem[j]>nd.hi[j]) return {INF,-1};
    uint64_t key=pack(rem); auto it=nd.memo.find(key); if(it!=nd.memo.end()) return it->second;
    Val best{INF,-1};
    for(int c=0;c<int(groups[nd.item].size());c++) {
        const auto &ch=groups[nd.item][c]; array<int,7> rest=rem; bool ok=true;
        for(int j=0;j<D;j++){rest[j]-=ch.h[j]; if(rest[j]<0) {ok=false;break;}}
        if(!ok) continue;
        auto v=solve(nd.next,rest); if(v.cost!=INF && v.cost+ch.cost<best.cost) best={v.cost+ch.cost,c};
    }
    nd.memo.emplace(key,best); return best;
}
int main(int argc,char**argv) {
    assert(argc==3); ifstream in(argv[1]); ofstream out(argv[2]); int batches; in>>D>>batches; assert(D>0&&D<=7);
    for(int b=0;b<batches;b++) {
        int ng,np; in>>ng>>np; groups.assign(ng,{}); nodes.clear(); nodes.push_back(Node{});
        for(auto &g:groups) {int n;in>>n;for(int a=0;a<n;a++){Choice c;in>>c.id>>c.cost;assert(abs(c.cost)<INF/100);for(int j=0;j<D;j++)in>>c.h[j];g.push_back(c);}}
        map<pair<int,int>,int> intern;
        struct Pattern {int id,s;I cost;array<int,7> target{};}; vector<Pattern> patterns;
        for(int p=0;p<np;p++) {
            Pattern pat;int n;in>>pat.id>>pat.cost>>n;vector<int> seq(n);for(int &v:seq)in>>v;for(int j=0;j<D;j++)in>>pat.target[j];
            int s=0;for(int a=n-1;a>=0;a--){pair<int,int> key{seq[a],s};auto it=intern.find(key);if(it!=intern.end()){s=it->second;continue;}
                Node nd;nd.item=seq[a];nd.next=s;
                for(int j=0;j<D;j++){int lo=1000,hi=-1000;for(auto &c:groups[nd.item]){lo=min(lo,c.h[j]);hi=max(hi,c.h[j]);}nd.lo[j]=lo+nodes[s].lo[j];nd.hi[j]=hi+nodes[s].hi[j];}
                int ns=nodes.size();nodes.push_back(move(nd));intern[key]=ns;s=ns;
            }pat.s=s;patterns.push_back(pat);
        }
        size_t feasible=0;
        for(auto &p:patterns){auto v=solve(p.s,p.target);out<<p.id<<' ';if(v.cost==INF){out<<"NONE\n";continue;}feasible++;out<<v.cost+p.cost;int s=p.s;auto rem=p.target;
            while(s){auto a=solve(s,rem);auto &nd=nodes[s];auto &c=groups[nd.item][a.choice];out<<' '<<c.id;for(int j=0;j<D;j++)rem[j]-=c.h[j];s=nd.next;}out<<'\n';
        }
        size_t states=0;for(auto &nd:nodes)states+=nd.memo.size();cerr<<"batch "<<b<<" patterns "<<np<<" feasible "<<feasible<<" states "<<states<<'\n';
    }
    assert(in.good());
}
