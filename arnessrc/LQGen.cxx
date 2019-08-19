#include "UHH2/LQToTopMu/include/LQGen.h"

using namespace std;
using namespace uhh2;

LQGen::LQGen(const vector<GenParticle> & genparticles, bool throw_on_failure)/*: m_type(e_notfound)*/ {    
    int n_LQ = 0, n_antiLQ = 0;
    for(unsigned int i=0; i<genparticles.size(); ++i) {
        const GenParticle & genp = genparticles[i];
        if (abs(genp.pdgId()) == 42){ // 42 = LQ's
            auto top = genp.daughter(&genparticles, 1);
            auto mu = genp.daughter(&genparticles, 2);
            if(!top || !mu){
                if(throw_on_failure) throw runtime_error("LQGen: LQ has not ==2 daughters");
                return;
            }
            if(abs(top->pdgId()) != 6){
                std::swap(top, mu);
            }
            if(abs(top->pdgId()) != 6){
                if(throw_on_failure) throw runtime_error("LQGen: LQ has no top daughter");
                return;
            }
            
            if(abs(mu->pdgId()) != 13){
                if(throw_on_failure) throw runtime_error("LQGen: LQ has no muon daughter");
                return;
            }
            // now get W daughters:
            auto topd1 = top->daughter(&genparticles, 1);
            auto topd2 = top->daughter(&genparticles, 2);
            if(!topd1 || !topd2){
                if(throw_on_failure) throw runtime_error("LQGen: top has not ==2 daughters");
                return;
            }
            
            // now that we collected everything, fill the member variables. 
            // Use different member variables according to LQ charge.
            if(genp.pdgId() == 42){
                m_LQ = genp;
                m_TopLQ = *top;
                m_muLQ = *mu;
                m_Topdecay1 = *topd1;
                m_Topdecay2 = *topd2;
                ++n_LQ;
            }
            else{
                m_AntiLQ = genp;
                m_TopAntiLQ = *top;
                m_muAntiLQ = *mu;
                m_Antitopdecay1 = *topd1;
                m_Antitopdecay2 = *topd2;
                ++n_antiLQ;
            }
        }
    }
    if(n_LQ != 1 || n_antiLQ != 1){
        if(throw_on_failure)  throw runtime_error("LQGen: did not find exactly one LQ and one antiLQ in the event");
        return;
    }
}   


GenParticle LQGen::LQ() const{
    return m_LQ;
}

GenParticle LQGen::AntiLQ() const{
    return m_AntiLQ;
} 

GenParticle LQGen::TopLQ() const{
    return m_TopLQ;
}

GenParticle LQGen::TopAntiLQ() const{
    return m_TopAntiLQ;
}

GenParticle LQGen::muLQ() const{
    return m_muLQ;
}

GenParticle LQGen::muAntiLQ() const{
    return m_muAntiLQ;
} 

GenParticle LQGen::Topdecay1() const{
    return m_Topdecay1;
} 

GenParticle LQGen::Topdecay2() const{
    return m_Topdecay2;
} 

GenParticle LQGen::Antitopdecay1() const{
    return m_Antitopdecay1;
} 

GenParticle LQGen::Antitopdecay2() const{
    return m_Antitopdecay2;
} 


LQGenProducer::LQGenProducer(uhh2::Context & ctx, const std::string & name, bool throw_on_failure_): throw_on_failure(throw_on_failure_){
    h_LQLQbargen = ctx.get_handle<LQGen>(name);
}

bool LQGenProducer::process(Event & event){
    event.set(h_LQLQbargen, LQGen(*event.genparticles, throw_on_failure));
    return true;
}
