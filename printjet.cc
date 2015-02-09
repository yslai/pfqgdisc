#include <cstdio>
#include <cmath>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>

namespace {

	void q_g_discr_raw(float *jet_multiplicity,
					   float *jet_multiplicity_charged,
					   float *jet_multiplicity_neutral,
					   float *jet_sigma_1, float *jet_sigma_2,
					   float *jet_ptd, float *jet_pull,
					   float *jet_fraction_max,
					   const float *jet_pseudorapidity,
					   const float *jet_azimuth, const size_t njet,
					   const int *pf_id, const float *pf_perp,
					   const float *pf_pseudorapidity,
					   const float *pf_azimuth, const size_t npf,
					   const float jet_r)
	{
		for (size_t i = 0; i < njet; i++) {
			jet_multiplicity[i] = 0;
			jet_multiplicity_charged[i] = 0;
			jet_multiplicity_neutral[i] = 0;

			float sum_pf_perp = 0;
			float sum_pf_perp_square = 0;

			float m11 = 0;
			float m22 = 0;
			float m12 = 0;

			float sum_r_pseudorapidity = 0;
			float sum_r_azimuth = 0;

			float max_pf_perp = 0;

			for (size_t j = 0; j < npf; j++) {
				const float dpseudorapidity =
					pf_pseudorapidity[j] - jet_pseudorapidity[i];
				const float dazimuth = pf_azimuth[j] - jet_azimuth[i];
				const float dpseudorapidity_square =
					dpseudorapidity * dpseudorapidity;
				const float dazimuth_square = dazimuth * dazimuth;
				const float dr_square =
					dpseudorapidity_square + dazimuth_square;

				if (dr_square < jet_r * jet_r) {
					const float pf_perp_square =
						pf_perp[j] * pf_perp[j];

					jet_multiplicity[i]++;
					if (pf_id[j] >= 1 && pf_id[j] <= 3) {
						jet_multiplicity_charged[i]++;
					}
					else {
						jet_multiplicity_neutral[i]++;
					}

					sum_pf_perp += pf_perp[j];
					sum_pf_perp_square += pf_perp_square;

					m11 += pf_perp_square * dpseudorapidity_square;
					m22 += pf_perp_square * dazimuth_square;
					m12 += -pf_perp_square * dpseudorapidity * dazimuth;

					const float dr = sqrt(dr_square);

					sum_r_pseudorapidity +=
						pf_perp_square * dr * dpseudorapidity;
					sum_r_azimuth += pf_perp_square * dr * dazimuth;

					max_pf_perp = std::max(max_pf_perp, pf_perp[j]);
				}
			}

			const float b = m11 + m22;
			const float c = m11 * m22 - m12 * m12;
			const float q = 0.5F * (b + copysign(sqrt(std::max(
				0.0F, b * b - 4.0F * c)), b));

			// Note that by construction, b >= 0, in Vieta's form,
			// lambda_2 = c / q

			const float lambda_1 = q;
			const float lambda_2 = q == 0 ? 0.0F : c / q;

			jet_sigma_1[i] = sum_pf_perp_square > 0 ?
				sqrt(std::max(0.0F, lambda_1) / sum_pf_perp_square) :
				0;
			jet_sigma_2[i] = sum_pf_perp_square > 0 ?
				sqrt(std::max(0.0F, lambda_2) / sum_pf_perp_square) :
				0;
			jet_ptd[i] = sum_pf_perp > 0 ?
				sqrt(sum_pf_perp_square) / sum_pf_perp :
				0;

			jet_pull[i] = sum_pf_perp_square > 0 ?
				sqrt(sum_r_pseudorapidity * sum_r_pseudorapidity +
					 sum_r_azimuth * sum_r_azimuth) /
				sum_pf_perp_square :
				0;
			jet_fraction_max[i] = sum_pf_perp > 0 ?
				max_pf_perp / sum_pf_perp : 0;
		}
	}
}

void printjet(TString filename = "")
{
	fprintf(stderr, "<F %s\n", (const char *)filename);

	TFile *root_file = TFile::Open(filename);

	TTree *skim_tree = dynamic_cast<TTree *>(
		root_file->Get("skimanalysis/HltTree"));
	TTree *event_tree = dynamic_cast<TTree *>(
		root_file->Get("hiEvtAnalyzer/HiTree"));
	TTree *pf_tree = dynamic_cast<TTree *>(
		root_file->Get("pfcandAnalyzer/pfTree"));

	int npf;
	int pf_id[32768];
	float pf_perp[32768];
	float pf_pseudorapidity[32768];
	float pf_azimuth[32768];

	if (pf_tree != NULL) {
		pf_tree->SetBranchAddress("nPFpart", &npf);
		pf_tree->SetBranchAddress("pfId", pf_id);
		pf_tree->SetBranchAddress("pfPt", pf_perp);
		pf_tree->SetBranchAddress("pfEta", pf_pseudorapidity);
		pf_tree->SetBranchAddress("pfPhi", pf_azimuth);
	}

	int collision_event_selection;
	int hb_he_noise_filter;
	int centrality_bin_cms;

	skim_tree->SetBranchAddress("pcollisionEventSelection",
								&collision_event_selection);
	skim_tree->SetBranchAddress("pHBHENoiseFilter", &hb_he_noise_filter);
	event_tree->SetBranchAddress("hiBin", &centrality_bin_cms);

	fprintf(stderr, "%s:%d:\n", __FILE__, __LINE__);

	const size_t nanalyzer = 42;
	const char *jet_analyzer[nanalyzer] = {
		"akVs1CaloJetAnalyzer",
		"akPu1CaloJetAnalyzer",
		"akVs1PFJetAnalyzer",
		"akPu1PFJetAnalyzer",
		"ak1CaloJetAnalyzer",
		"ak1PFJetAnalyzer",
		"akVs2CaloJetAnalyzer",
		"akPu2CaloJetAnalyzer",
		"akVs2PFJetAnalyzer",
		"akPu2PFJetAnalyzer",
		"ak2CaloJetAnalyzer",
		"ak2PFJetAnalyzer",
		"akVs3CaloJetAnalyzer",
		"akPu3CaloJetAnalyzer",
		"akVs3PFJetAnalyzer",
		"akPu3PFJetAnalyzer",
		"ak3CaloJetAnalyzer",
		"ak3PFJetAnalyzer",
		"akVs4CaloJetAnalyzer",
		"akPu4CaloJetAnalyzer",
		"akVs4PFJetAnalyzer",
		"akPu4PFJetAnalyzer",
		"ak4CaloJetAnalyzer",
		"ak4PFJetAnalyzer",
		"akVs5CaloJetAnalyzer",
		"akPu5CaloJetAnalyzer",
		"akVs5PFJetAnalyzer",
		"akPu5PFJetAnalyzer",
		"ak5CaloJetAnalyzer",
		"ak5PFJetAnalyzer",
		"akVs6CaloJetAnalyzer",
		"akPu6CaloJetAnalyzer",
		"akVs6PFJetAnalyzer",
		"akPu6PFJetAnalyzer",
		"ak6CaloJetAnalyzer",
		"ak6PFJetAnalyzer",
		"akVs7CaloJetAnalyzer",
		"akPu7CaloJetAnalyzer",
		"akVs7PFJetAnalyzer",
		"akPu7PFJetAnalyzer",
		"ak7CaloJetAnalyzer",
		"ak7PFJetAnalyzer",
	};

	TTree *jet_tree[nanalyzer];
	float jet_r[nanalyzer];

	for (size_t i = 0; i < nanalyzer; i++) {
		char buf[BUFSIZ];

		snprintf(buf, BUFSIZ, "%s/t", jet_analyzer[i]);
		jet_tree[i] = dynamic_cast<TTree *>(root_file->Get(buf));

		for (const char *p = jet_analyzer[i]; p != '\0'; p++) {
			if (*p >= '0' && *p <= '9') {
				jet_r[i] = 0.1F * static_cast<float>(*p - '0');
				break;
			}
		}

		fprintf(stderr, "%s:%d: %s %p %g\n", __FILE__, __LINE__,
				buf, jet_tree[i], jet_r[i]);
	}

	float jet_perp_calib[nanalyzer][4096];
	float jet_perp[nanalyzer][4096];
	float jet_reference_perp[nanalyzer][4096];
	float jet_pseudorapidity[nanalyzer][4096];
	float jet_azimuth[nanalyzer][4096];
	int jet_reference_flavor[nanalyzer][4096];
	int jet_subevent_id[nanalyzer][4096];
	int njet[nanalyzer];

	for (size_t i = 0; i < nanalyzer; i++) {
		if (jet_tree[i] != NULL) {
			jet_tree[i]->SetBranchAddress("jtpt", jet_perp_calib[i]);
			jet_tree[i]->SetBranchAddress("rawpt", jet_perp[i]);
			jet_tree[i]->SetBranchAddress("refpt",
										  jet_reference_perp[i]);
			jet_tree[i]->SetBranchAddress("jteta",
										  jet_pseudorapidity[i]);
			jet_tree[i]->SetBranchAddress("jtphi", jet_azimuth[i]);
			jet_tree[i]->SetBranchAddress("refparton_flavor",
										  jet_reference_flavor[i]);
			jet_tree[i]->SetBranchAddress("subid",
										  jet_subevent_id[i]);
			jet_tree[i]->SetBranchAddress("nref", &njet[i]);
		}
	}

	fprintf(stderr, "%s:%d:\n", __FILE__, __LINE__);

	long int nevent = skim_tree->GetEntries();

	for (long int i = 0; i < nevent; i++) {
		skim_tree->GetEntry(i);
		event_tree->GetEntry(i);

		if (//collision_event_selection &&
			hb_he_noise_filter) {
			fprintf(stderr, "<E %g %d %d\n",
					centrality_bin_cms * 0.5,
					collision_event_selection, hb_he_noise_filter);
			if (pf_tree != NULL) {
				pf_tree->GetEntry(i);
			}
			for (size_t j = 0; j < nanalyzer; j++) {
				if (jet_tree[j] != NULL) {
					jet_tree[j]->GetEntry(i);

					float jet_multiplicity[4096];
					float jet_multiplicity_charged[4096];
					float jet_multiplicity_neutral[4096];
					float jet_sigma_1[4096];
					float jet_sigma_2[4096];
					float jet_ptd[4096];
					float jet_pull[4096];
					float jet_fraction_max[4096];

					q_g_discr_raw(jet_multiplicity,
								  jet_multiplicity_charged,
								  jet_multiplicity_neutral,
								  jet_sigma_1,
								  jet_sigma_2,
								  jet_ptd,
								  jet_pull,
								  jet_fraction_max,
								  jet_pseudorapidity[j],
								  jet_azimuth[j], njet[j],
								  pf_id,
								  pf_perp, pf_pseudorapidity,
								  pf_azimuth,
								  pf_tree == NULL ? 0 : npf,
								  jet_r[j]);
 
					size_t matched_count = 0;

					for (int k = 0; k < njet[j]; k++) {
						if (jet_reference_perp[j][k] >= 0 &&
							jet_subevent_id[j][k] == 0) {
							matched_count++;
						}
					}

					if (true || matched_count > 0) {
						fprintf(stderr, "<A %lu %d\n", j, njet[j]);
						for (int k = 0; k < njet[j]; k++) {
							if ((jet_perp_calib[j][k] >= 30 ||
								 jet_reference_perp[j][k] >= 0) &&
								jet_subevent_id[j][k] == 0) {
								fprintf(stderr,
										"<J %.8e %.8e %.8e %.8e %d "
										"%d\n",
										jet_perp[j][k],
										jet_reference_perp[j][k],
										jet_pseudorapidity[j][k],
										jet_azimuth[j][k],
										jet_reference_flavor[j][k],
										jet_subevent_id[j][k]);
								fprintf(stderr,
										"<C %.8e\n",
										jet_perp_calib[j][k]);
								fprintf(stderr,
										"<S %g %g %g %.8e %.8e %.8e "
										"%.8e %.8e\n",
										jet_multiplicity[k],
										jet_multiplicity_charged[k],
										jet_multiplicity_neutral[k],
										jet_sigma_1[k],
										jet_sigma_2[k],
										jet_ptd[k],
										jet_pull[k],
										jet_fraction_max[k]);
							}
						}
					}
				}
			}
		}
	}

	gSystem->Exit(0);
}
