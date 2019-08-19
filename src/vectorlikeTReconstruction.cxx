#include "UHH2/vectorlikeT/arnesinclude/SingleVLQReconstruction.h"
#include "UHH2/vectorlikeT/include/vectorlikeTReconstruction.h"
#include "UHH2/core/include/LorentzVector.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/core/include/Utils.h"
#include "UHH2/common/include/JetIds.h"

#include <cassert>

using namespace uhh2;
using namespace std;




HighMassvectorlikeTReconstruction::HighMassvectorlikeTReconstruction(Context & ctx, const NeutrinoReconstructionMethod & neutrinofunction, CSVBTag::wp wp_btag_): m_neutrinofunction(neutrinofunction), wp_btag(wp_btag_) {

	h_recohyps = ctx.get_handle<vector<vectorlikeTReconstructionHypothesis>>("TprimeHypotheses");
	h_is_tprime_reco = ctx.get_handle<bool>("is_tprime_reco");

}

HighMassvectorlikeTReconstruction::~HighMassvectorlikeTReconstruction() {}

bool HighMassvectorlikeTReconstruction::process(uhh2::Event & event) {
	assert(event.jets);
	assert(event.met);

	//find primary charged lepton
	const Particle & lepton = event.muons->at(0);
	LorentzVector lepton_v4 = lepton.v4();
	std::vector<vectorlikeTReconstructionHypothesis> recoHyps;

	//reconstruct neutrino
	std::vector<LorentzVector> neutrinos = m_neutrinofunction( lepton.v4(), event.met->v4());
	std::vector<Jet> btags;
	CSVBTag Btag = CSVBTag(wp_btag);
	for (unsigned int i =0; i<event.jets->size(); ++i) {
		if(Btag(event.jets->at(i),event)) btags.push_back(event.jets->at(i));
	}
	unsigned int n_jets = btags.size();
	// if(n_jets>3) n_jets=3;
	const unsigned int max_j = pow(3, n_jets);

	//loop over neutrino solutions and jet assignments to fill hyotheses
	for(const auto & neutrino_p4 : neutrinos) {
		const LorentzVector wlep_v4 = lepton.v4() + neutrino_p4;
		for (unsigned int j=0; j < max_j; j++) {
			LorentzVector higgs_v4;
			LorentzVector toplep_v4 = wlep_v4;
			int hadjets=0;
			int lepjets=0;
			int num = j;
			vectorlikeTReconstructionHypothesis hyp;

			vector<Particle> had_jets, lep_jets;
			vector<Jet> hadrjets, leptjets;
			for (unsigned int k=0; k<n_jets; k++) {
				if(num%3==0) {
					higgs_v4 = higgs_v4 + btags.at(k).v4();
					had_jets.push_back(btags.at(k));
					hadrjets.push_back(btags.at(k));
					hadjets++;
				}

				if(num%3==1) {
					toplep_v4 = toplep_v4 + btags.at(k).v4();
					lep_jets.push_back(btags.at(k));
					leptjets.push_back(btags.at(k));
					lepjets++;
				}
				//in case num%3==2 do not take this jet at all
				//shift the trigits of num to the right:
				num /= 3;
			}



			//fill only hypotheses with at least one jet assigned to each top quark
			//if(hadjets>0 && lepjets>0 && hadbtag && lepbtag) {
			if(hadjets==2 && lepjets==1) {

				hyp.set_higgs_v4(higgs_v4);
				hyp.set_toplep_v4(toplep_v4);
				hyp.set_lepton(lepton);
				hyp.set_lepton_v4(lepton_v4);
				hyp.set_neutrino_v4(neutrino_p4);
				hyp.set_higgs_jets(had_jets);
				hyp.set_toplep_jets(lep_jets);

				recoHyps.emplace_back(move(hyp));
			} // if at least 1 jet is assigned to each top quark
		} // 3^n_jets jet combinations * n_muon muon combinations
	} // neutrinos

	event.set(h_recohyps, move(recoHyps));
	event.set(h_is_tprime_reco, true);
	return true;
}




std::vector<LorentzVector> vectorlikeTNeutrinoReconstruction(const LorentzVector & lepton, const LorentzVector & met) {
	TVector3 lepton_pT = toVector(lepton);
	lepton_pT.SetZ(0);
	TVector3 neutrino_pT = toVector(met);
	neutrino_pT.SetZ(0);
	constexpr float mass_w = 80.399f;
	float mu = mass_w * mass_w / 2 + lepton_pT * neutrino_pT;
	float A = - (lepton_pT * lepton_pT);
	float B = mu * lepton.pz();
	float C = mu * mu - lepton.e() * lepton.e() * (neutrino_pT * neutrino_pT);
	float discriminant = B * B - A * C;
	std::vector<LorentzVector> solutions;
	if (0 >= discriminant) {
		// Take only real part of the solution for pz:
		LorentzVectorXYZE solution (met.Px(),met.Py(),-B / A,0);
		solution.SetE(solution.P());
		solutions.emplace_back(toPtEtaPhi(solution));
	}
	else {
		discriminant = sqrt(discriminant);
		LorentzVectorXYZE solution (met.Px(),met.Py(),(-B - discriminant) / A,0);
		solution.SetE(solution.P());
		solutions.emplace_back(toPtEtaPhi(solution));

		LorentzVectorXYZE solution2 (met.Px(),met.Py(),(-B + discriminant) / A,0);
		solution2.SetE(solution2.P());
		solutions.emplace_back(toPtEtaPhi(solution2));
	}
	return solutions;
}
