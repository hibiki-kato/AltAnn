#include "altann.hpp"

void run_case(const string &left, const string &right, bool invalid, bool overlap=false) {
    string prefix(31, 'C');
    string first = overlap ? "AT" : "ATG" + string(120, 'C') + left;
    string intron = "GT" + string(36, 'C') + "AG";
    string second = right + string(120, 'C') + "TAA";
    string dna = prefix + first + intron + second + string(35, 'C');
    vector<char> seq(dna.begin(), dna.end());
    vector<array<float, NUM_STATES>> emissions(seq.size());
    vector<double> site(seq.size(), 0);
    auto trans = init_transitions();
    ModelInputs inputs{{emissions.data(),0,(int)seq.size()}, {site.data(),0,(int)seq.size()},
        {site.data(),0,(int)seq.size()}, {site.data(),0,(int)seq.size()},
        {site.data(),0,(int)seq.size()},seq,trans};
    int start=prefix.size(), donor=start+first.size(), acceptor=donor+intron.size()-2;
    int end=acceptor+2+second.size()-1;
    int state=0; PathMetadata metadata; metadata.inter_len=1;
    bool rejected=false;
    for(int p=1;p<(int)seq.size();++p) {
        int to=state;
        if(p==start+2) to=1+start%3;
        if(p==donor+1) to=4+start%3;
        if(p==acceptor+1) to=1+(start+(int)intron.size())%3;
        if(p==end) to=0;
        auto step=evaluate_transition(p,state,to,metadata,inputs);
        if(!step.allowed) {rejected=true; break;}
        metadata=step.next_metadata; state=to;
    }
    if(rejected!=invalid) throw runtime_error("unexpected spliced validity: "+left+"/"+right);
}
int main(){
    for(string stop: {"TAA","TAG","TGA"}) for(int split: {1,2})
        run_case(stop.substr(0,split),stop.substr(split),true);
    for(string coding: {"CAA","TAC","TGG","AAA"}) for(int split: {1,2})
        run_case(coding.substr(0,split),coding.substr(split),false);
    run_case("","CCC",false);
    run_case("","C",true,true);
    cout << "16 spliced codon/start-boundary transition cases passed\n";
}
