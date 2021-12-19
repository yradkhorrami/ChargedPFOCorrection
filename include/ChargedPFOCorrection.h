#ifndef ChargedPFOCorrection_h
#define ChargedPFOCorrection_h 1
#include <marlin/Processor.h>
#include <marlin/Global.h>
#include "UTIL/LCRelationNavigator.h"
#include "UTIL/LCIterator.h"
#include "UTIL/Operators.h"
#include <UTIL/BitField64.h>
#include "UTIL/LCTrackerConf.h"
#include <UTIL/ILDConf.h>
#include <UTIL/BitSet32.h>
#include "EVENT/LCStrVec.h"
#include <EVENT/MCParticle.h>
#include <EVENT/Track.h>
#include "IMPL/LCCollectionVec.h"
#include <IMPL/ReconstructedParticleImpl.h>
#include "lcio.h"
#include <string>
#include <vector>
#include <math.h>
#include <set>
#include <vector>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMatrixD.h"
#include "TStyle.h"
#include "TCanvas.h"

class TFile;
class TDirectory;
class TH1F;
class TH1I;
class TH2I;
class TH2F;
class TTree;
class TF1;

using namespace lcio ;
using namespace marlin ;

class ChargedPFOCorrection : public Processor
{
	public:

		virtual Processor*  newProcessor()
		{
			return new ChargedPFOCorrection;
		}
		ChargedPFOCorrection();
		virtual ~ChargedPFOCorrection() = default;
		ChargedPFOCorrection(const ChargedPFOCorrection&) = delete;
		ChargedPFOCorrection& operator=(const ChargedPFOCorrection&) = delete;
		virtual void init();
		virtual void Clear();
		virtual void processRunHeader();
		virtual void processEvent( EVENT::LCEvent *pLCEvent );
		std::vector<double> getTruthTrkID( EVENT::LCEvent *pLCEvent , EVENT::Track* inputTrk );
		int getTrackIndex( EVENT::LCCollection *TrackCollection , EVENT::Track* inputTrk );
		TLorentzVector getMCPFourMomentum( EVENT::LCEvent *pLCEvent, EVENT::Track* inputTrk );
		TLorentzVector getTrackFourMomentum( EVENT::Track* inputTrk , double trackMass );
		std::vector<float> getPFOResidual( TLorentzVector pfoFourMomentum , TLorentzVector mcpFourMomentum );
		std::vector<float> UpdateChargedPFOCovMat( EVENT::Track* inputTrk , float trackMass );
		virtual void updatePFO( EVENT::ReconstructedParticle* inputPFO , ReconstructedParticleImpl* outputPFO , std::vector<Track*> outputPFOtrkvec , TLorentzVector pfoFourMomentum , std::vector<float> pfoCovMat , int pfoType );
		std::vector<float> getAngularUncertainties( TLorentzVector pfoFourMomentum , std::vector<float> pfoCovMat );
//		virtual void InitializeHistogram( TH1F *histogram , std::int scale , std::int color , std::int lineWidth , std::int markerSize , std::int markerStyle );
		virtual void InitializeHistogram( TH1F *histogram , int scale , int color , int lineWidth , int markerSize , int markerStyle );
		virtual void doProperGaussianFit( TH1F *histogram , float fitMin , float fitMax , float fitRange );
		virtual void MakeRatioPlots( TCanvas *Canvas , TH1F *histogram1 , TH1F *histogram2 , int scale1 , int scale2 , int color1 , int color2 , float Ratio_min , float Ratio_max , float XTitleSize , float XTitleOffset );
		virtual void check( EVENT::LCEvent *pLCEvent );
		virtual void end();

	private:

		typedef std::vector<int>		IntVector;
		typedef std::vector<double>		DoubleVector;
		typedef std::vector<float>		FloatVector;

		std::string				m_mcParticleCollection{};
		std::string				m_inputPfoCollection{};
		std::string				m_MarlinTrkTracks{};
		std::string				m_MarlinTrkTracksKAON{};
		std::string				m_MarlinTrkTracksPROTON{};
		std::string				m_TrackMCTruthLinkCollection{};
		std::string				m_MCTruthTrackLinkCollection{};
		std::string				m_outputPfoCollection{};
		std::string				m_rootFile{};

		bool					m_updatePFOwithOneTrack = true;
		bool					m_updatePFOwithTwoTrack = true;
		bool					m_updatePFOwithMoreTrack = true;
		bool					m_updatePFOwithPionTrack = true;
		bool					m_fillRootTree = false;

		int					m_nRun;
		int					m_nEvt;
		int					m_nRunSum;
		int					m_nEvtSum;
		int					n_OldPFOS1trk_NormalizedResidualPx;
		int					n_OldPFOS1trk_NormalizedResidualPxPy;
		int					n_OldPFOS1trk_NormalizedResidualPy;
		int					n_OldPFOS1trk_NormalizedResidualPxPz;
		int					n_OldPFOS1trk_NormalizedResidualPyPz;
		int					n_OldPFOS1trk_NormalizedResidualPz;
		int					n_OldPFOS1trk_NormalizedResidualPxE;
		int					n_OldPFOS1trk_NormalizedResidualPyE;
		int					n_OldPFOS1trk_NormalizedResidualPzE;
		int					n_OldPFOS1trk_NormalizedResidualE;
		int					n_OldPFOS1trk_NormalizedResidualTheta;
		int					n_OldPFOS1trk_NormalizedResidualPhi;
		int					n_NewPFOS1trk_NormalizedResidualPx;
		int					n_NewPFOS1trk_NormalizedResidualPxPy;
		int					n_NewPFOS1trk_NormalizedResidualPy;
		int					n_NewPFOS1trk_NormalizedResidualPxPz;
		int					n_NewPFOS1trk_NormalizedResidualPyPz;
		int					n_NewPFOS1trk_NormalizedResidualPz;
		int					n_NewPFOS1trk_NormalizedResidualPxE;
		int					n_NewPFOS1trk_NormalizedResidualPyE;
		int					n_NewPFOS1trk_NormalizedResidualPzE;
		int					n_NewPFOS1trk_NormalizedResidualE;
		int					n_NewPFOS1trk_NormalizedResidualTheta;
		int					n_NewPFOS1trk_NormalizedResidualPhi;
		int					n_NewPFOS2trk_NormalizedResidualPx;
		int					n_NewPFOS2trk_NormalizedResidualPxPy;
		int					n_NewPFOS2trk_NormalizedResidualPy;
		int					n_NewPFOS2trk_NormalizedResidualPxPz;
		int					n_NewPFOS2trk_NormalizedResidualPyPz;
		int					n_NewPFOS2trk_NormalizedResidualPz;
		int					n_NewPFOS2trk_NormalizedResidualPxE;
		int					n_NewPFOS2trk_NormalizedResidualPyE;
		int					n_NewPFOS2trk_NormalizedResidualPzE;
		int					n_NewPFOS2trk_NormalizedResidualE;
		int					n_NewPFOS2trk_NormalizedResidualTheta;
		int					n_NewPFOS2trk_NormalizedResidualPhi;
		int					n_TrueProtonLinkWeight;
		int					n_StdTrkProton_NormalizedResidualPx;
		int					n_StdTrkProton_NormalizedResidualPxPy;
		int					n_StdTrkProton_NormalizedResidualPy;
		int					n_StdTrkProton_NormalizedResidualPxPz;
		int					n_StdTrkProton_NormalizedResidualPyPz;
		int					n_StdTrkProton_NormalizedResidualPz;
		int					n_StdTrkProton_NormalizedResidualPxE;
		int					n_StdTrkProton_NormalizedResidualPyE;
		int					n_StdTrkProton_NormalizedResidualPzE;
		int					n_StdTrkProton_NormalizedResidualE;
		int					n_StdTrkProton_NormalizedResidualTheta;
		int					n_StdTrkProton_NormalizedResidualPhi;
		int					n_RFTrkProton_NormalizedResidualPx;
		int					n_RFTrkProton_NormalizedResidualPxPy;
		int					n_RFTrkProton_NormalizedResidualPy;
		int					n_RFTrkProton_NormalizedResidualPxPz;
		int					n_RFTrkProton_NormalizedResidualPyPz;
		int					n_RFTrkProton_NormalizedResidualPz;
		int					n_RFTrkProton_NormalizedResidualPxE;
		int					n_RFTrkProton_NormalizedResidualPyE;
		int					n_RFTrkProton_NormalizedResidualPzE;
		int					n_RFTrkProton_NormalizedResidualE;
		int					n_RFTrkProton_NormalizedResidualTheta;
		int					n_RFTrkProton_NormalizedResidualPhi;
		int					n_TrueKaonLinkWeight;
		int					n_StdTrkKaon_NormalizedResidualPx;
		int					n_StdTrkKaon_NormalizedResidualPxPy;
		int					n_StdTrkKaon_NormalizedResidualPy;
		int					n_StdTrkKaon_NormalizedResidualPxPz;
		int					n_StdTrkKaon_NormalizedResidualPyPz;
		int					n_StdTrkKaon_NormalizedResidualPz;
		int					n_StdTrkKaon_NormalizedResidualPxE;
		int					n_StdTrkKaon_NormalizedResidualPyE;
		int					n_StdTrkKaon_NormalizedResidualPzE;
		int					n_StdTrkKaon_NormalizedResidualE;
		int					n_StdTrkKaon_NormalizedResidualTheta;
		int					n_StdTrkKaon_NormalizedResidualPhi;
		int					n_RFTrkKaon_NormalizedResidualPx;
		int					n_RFTrkKaon_NormalizedResidualPxPy;
		int					n_RFTrkKaon_NormalizedResidualPy;
		int					n_RFTrkKaon_NormalizedResidualPxPz;
		int					n_RFTrkKaon_NormalizedResidualPyPz;
		int					n_RFTrkKaon_NormalizedResidualPz;
		int					n_RFTrkKaon_NormalizedResidualPxE;
		int					n_RFTrkKaon_NormalizedResidualPyE;
		int					n_RFTrkKaon_NormalizedResidualPzE;
		int					n_RFTrkKaon_NormalizedResidualE;
		int					n_RFTrkKaon_NormalizedResidualTheta;
		int					n_RFTrkKaon_NormalizedResidualPhi;
		int					n_StdTrk_NormalizedResidualPx;
		int					n_StdTrk_NormalizedResidualPxPy;
		int					n_StdTrk_NormalizedResidualPy;
		int					n_StdTrk_NormalizedResidualPxPz;
		int					n_StdTrk_NormalizedResidualPyPz;
		int					n_StdTrk_NormalizedResidualPz;
		int					n_StdTrk_NormalizedResidualPxE;
		int					n_StdTrk_NormalizedResidualPyE;
		int					n_StdTrk_NormalizedResidualPzE;
		int					n_StdTrk_NormalizedResidualE;
		int					n_StdTrk_NormalizedResidualTheta;
		int					n_StdTrk_NormalizedResidualPhi;
		int					n_RFTrk_NormalizedResidualPx;
		int					n_RFTrk_NormalizedResidualPxPy;
		int					n_RFTrk_NormalizedResidualPy;
		int					n_RFTrk_NormalizedResidualPxPz;
		int					n_RFTrk_NormalizedResidualPyPz;
		int					n_RFTrk_NormalizedResidualPz;
		int					n_RFTrk_NormalizedResidualPxE;
		int					n_RFTrk_NormalizedResidualPyE;
		int					n_RFTrk_NormalizedResidualPzE;
		int					n_RFTrk_NormalizedResidualE;
		int					n_RFTrk_NormalizedResidualTheta;
		int					n_RFTrk_NormalizedResidualPhi;
		IntVector				m_pfoCharge{};
		IntVector				m_nTracksOfPFO{};
		FloatVector				m_TrkToMCPLinkWeight{};
		FloatVector				m_oldPFO_Px{};
		FloatVector				m_oldPFO_Py{};
		FloatVector				m_oldPFO_Pz{};
		FloatVector				m_oldPFO_E{};
		FloatVector				m_oldPFO_Theta{};
		FloatVector				m_oldPFO_Phi{};
		FloatVector				m_newPFO_Px{};
		FloatVector				m_newPFO_Py{};
		FloatVector				m_newPFO_Pz{};
		FloatVector				m_newPFO_E{};
		FloatVector				m_newPFO_Theta{};
		FloatVector				m_newPFO_Phi{};
		FloatVector				m_oldPFO_ResidualPx{};
		FloatVector				m_oldPFO_ResidualPy{};
		FloatVector				m_oldPFO_ResidualPz{};
		FloatVector				m_oldPFO_ResidualE{};
		FloatVector				m_oldPFO_ResidualTheta{};
		FloatVector				m_oldPFO_ResidualPhi{};
		FloatVector				m_newPFO_ResidualPx{};
		FloatVector				m_newPFO_ResidualPy{};
		FloatVector				m_newPFO_ResidualPz{};
		FloatVector				m_newPFO_ResidualE{};
		FloatVector				m_newPFO_ResidualTheta{};
		FloatVector				m_newPFO_ResidualPhi{};
		FloatVector				m_oldPFO_NormalizedResidualPx{};
		FloatVector				m_oldPFO_NormalizedResidualPy{};
		FloatVector				m_oldPFO_NormalizedResidualPz{};
		FloatVector				m_oldPFO_NormalizedResidualE{};
		FloatVector				m_oldPFO_NormalizedResidualTheta{};
		FloatVector				m_oldPFO_NormalizedResidualPhi{};
		FloatVector				m_newPFO_NormalizedResidualPx{};
		FloatVector				m_newPFO_NormalizedResidualPy{};
		FloatVector				m_newPFO_NormalizedResidualPz{};
		FloatVector				m_newPFO_NormalizedResidualE{};
		FloatVector				m_newPFO_NormalizedResidualTheta{};
		FloatVector				m_newPFO_NormalizedResidualPhi{};
//		IntVector				m_nClustersOfPFO{};
//		IntVector				m_nSubDetectorHits{};
		float					m_Bfield;
		double					c;
		double					mm2m;
		double					eV2GeV;
		double					eB;
		float					m_MinWeightTrackMCTruthLink;
		float					m_MinWeightMCTruthTrackLink;
		float					m_pion_mass;
		float					m_proton_mass;
		float					m_kaon_mass;
		TFile					*m_pTFile{};
		TTree					*m_pTTree{};
	        TTree					*m_pTTree_1trk{};
		TTree					*m_pTTree_2trk{};
		TTree					*m_pTTree_ntrk{};
		TDirectory				*m_Histograms{};
		TDirectory				*m_Plots{};
		TDirectory				*m_OldPFOs_1Trk{};
		TDirectory				*m_NewPFOs_1Trk{};
		TDirectory				*m_NewPFOs_2Trk{};
		TDirectory				*m_NewPFOs_nTrk{};
		TDirectory				*m_TrueProtons_1Trk{};
		TDirectory				*m_TrueKaons_1Trk{};
		TDirectory				*m_IndividualTracks_2Trk{};
		TCanvas					*c_1trk_NormalizedResidualPx{};
		TCanvas					*c_1trk_NormalizedResidualPxPy{};
		TCanvas					*c_1trk_NormalizedResidualPy{};
		TCanvas					*c_1trk_NormalizedResidualPxPz{};
		TCanvas					*c_1trk_NormalizedResidualPyPz{};
		TCanvas					*c_1trk_NormalizedResidualPz{};
		TCanvas					*c_1trk_NormalizedResidualPxE{};
		TCanvas					*c_1trk_NormalizedResidualPyE{};
		TCanvas					*c_1trk_NormalizedResidualPzE{};
		TCanvas					*c_1trk_NormalizedResidualE{};
		TCanvas					*c_1trk_NormalizedResidualTheta{};
		TCanvas					*c_1trk_NormalizedResidualPhi{};
		TCanvas					*c_Protons_NormalizedResidualPx{};
		TCanvas					*c_Protons_NormalizedResidualPxPy{};
		TCanvas					*c_Protons_NormalizedResidualPy{};
		TCanvas					*c_Protons_NormalizedResidualPxPz{};
		TCanvas					*c_Protons_NormalizedResidualPyPz{};
		TCanvas					*c_Protons_NormalizedResidualPz{};
		TCanvas					*c_Protons_NormalizedResidualPxE{};
		TCanvas					*c_Protons_NormalizedResidualPyE{};
		TCanvas					*c_Protons_NormalizedResidualPzE{};
		TCanvas					*c_Protons_NormalizedResidualE{};
		TCanvas					*c_Protons_NormalizedResidualTheta{};
		TCanvas					*c_Protons_NormalizedResidualPhi{};
		TCanvas					*c_Kaons_NormalizedResidualPx{};
		TCanvas					*c_Kaons_NormalizedResidualPxPy{};
		TCanvas					*c_Kaons_NormalizedResidualPy{};
		TCanvas					*c_Kaons_NormalizedResidualPxPz{};
		TCanvas					*c_Kaons_NormalizedResidualPyPz{};
		TCanvas					*c_Kaons_NormalizedResidualPz{};
		TCanvas					*c_Kaons_NormalizedResidualPxE{};
		TCanvas					*c_Kaons_NormalizedResidualPyE{};
		TCanvas					*c_Kaons_NormalizedResidualPzE{};
		TCanvas					*c_Kaons_NormalizedResidualE{};
		TCanvas					*c_Kaons_NormalizedResidualTheta{};
		TCanvas					*c_Kaons_NormalizedResidualPhi{};
		TH2I					*h_nClusters_nTracks{};
		TH2I					*h_pfoCharge_nTracks{};
		TH2F					*h_InnermostRadiusHit_Neutral{};
		TH2F					*h_InnermostRadiusHit_Charged{};
		TH2I					*h_FirstSubDet_Charged{};
		TH1F					*h_OldPFOS1trk_NormalizedResidualPx{};
		TH1F					*h_OldPFOS1trk_NormalizedResidualPxPy{};
		TH1F					*h_OldPFOS1trk_NormalizedResidualPy{};
		TH1F					*h_OldPFOS1trk_NormalizedResidualPxPz{};
		TH1F					*h_OldPFOS1trk_NormalizedResidualPyPz{};
		TH1F					*h_OldPFOS1trk_NormalizedResidualPz{};
		TH1F					*h_OldPFOS1trk_NormalizedResidualPxE{};
		TH1F					*h_OldPFOS1trk_NormalizedResidualPyE{};
		TH1F					*h_OldPFOS1trk_NormalizedResidualPzE{};
		TH1F					*h_OldPFOS1trk_NormalizedResidualE{};
		TH1F					*h_OldPFOS1trk_NormalizedResidualTheta{};
		TH1F					*h_OldPFOS1trk_NormalizedResidualPhi{};
		TH2F					*h_OldPFOS1trk_ResidualPxPy{};
		TH2F					*h_OldPFOS1trk_ResidualPxPz{};
		TH2F					*h_OldPFOS1trk_ResidualPyPz{};
		TH2F					*h_OldPFOS1trk_ResidualPxE{};
		TH2F					*h_OldPFOS1trk_ResidualPyE{};
		TH2F					*h_OldPFOS1trk_ResidualPzE{};
		TH1F					*h_NewPFOS1trk_NormalizedResidualPx{};
		TH1F					*h_NewPFOS1trk_NormalizedResidualPxPy{};
		TH1F					*h_NewPFOS1trk_NormalizedResidualPy{};
		TH1F					*h_NewPFOS1trk_NormalizedResidualPxPz{};
		TH1F					*h_NewPFOS1trk_NormalizedResidualPyPz{};
		TH1F					*h_NewPFOS1trk_NormalizedResidualPz{};
		TH1F					*h_NewPFOS1trk_NormalizedResidualPxE{};
		TH1F					*h_NewPFOS1trk_NormalizedResidualPyE{};
		TH1F					*h_NewPFOS1trk_NormalizedResidualPzE{};
		TH1F					*h_NewPFOS1trk_NormalizedResidualE{};
		TH1F					*h_NewPFOS1trk_NormalizedResidualTheta{};
		TH1F					*h_NewPFOS1trk_NormalizedResidualPhi{};
		TH2F					*h_NewPFOS1trk_ResidualPxPy{};
		TH2F					*h_NewPFOS1trk_ResidualPxPz{};
		TH2F					*h_NewPFOS1trk_ResidualPyPz{};
		TH2F					*h_NewPFOS1trk_ResidualPxE{};
		TH2F					*h_NewPFOS1trk_ResidualPyE{};
		TH2F					*h_NewPFOS1trk_ResidualPzE{};
		TH1F					*h_NewPFOS2trk_NormalizedResidualPx{};
		TH1F					*h_NewPFOS2trk_NormalizedResidualPxPy{};
		TH1F					*h_NewPFOS2trk_NormalizedResidualPy{};
		TH1F					*h_NewPFOS2trk_NormalizedResidualPxPz{};
		TH1F					*h_NewPFOS2trk_NormalizedResidualPyPz{};
		TH1F					*h_NewPFOS2trk_NormalizedResidualPz{};
		TH1F					*h_NewPFOS2trk_NormalizedResidualPxE{};
		TH1F					*h_NewPFOS2trk_NormalizedResidualPyE{};
		TH1F					*h_NewPFOS2trk_NormalizedResidualPzE{};
		TH1F					*h_NewPFOS2trk_NormalizedResidualE{};
		TH1F					*h_NewPFOS2trk_NormalizedResidualTheta{};
		TH1F					*h_NewPFOS2trk_NormalizedResidualPhi{};
		TH2F					*h_NewPFOS2trk_ResidualPxPy{};
		TH2F					*h_NewPFOS2trk_ResidualPxPz{};
		TH2F					*h_NewPFOS2trk_ResidualPyPz{};
		TH2F					*h_NewPFOS2trk_ResidualPxE{};
		TH2F					*h_NewPFOS2trk_ResidualPyE{};
		TH2F					*h_NewPFOS2trk_ResidualPzE{};
		TH1F					*h_TrueProtonLinkWeight{};
		TH1F					*h_StdTrkProton_NormalizedResidualPx{};
		TH1F					*h_StdTrkProton_NormalizedResidualPxPy{};
		TH1F					*h_StdTrkProton_NormalizedResidualPy{};
		TH1F					*h_StdTrkProton_NormalizedResidualPxPz{};
		TH1F					*h_StdTrkProton_NormalizedResidualPyPz{};
		TH1F					*h_StdTrkProton_NormalizedResidualPz{};
		TH1F					*h_StdTrkProton_NormalizedResidualPxE{};
		TH1F					*h_StdTrkProton_NormalizedResidualPyE{};
		TH1F					*h_StdTrkProton_NormalizedResidualPzE{};
		TH1F					*h_StdTrkProton_NormalizedResidualE{};
		TH1F					*h_StdTrkProton_NormalizedResidualTheta{};
		TH1F					*h_StdTrkProton_NormalizedResidualPhi{};
		TH2F					*h_StdTrkProton_ResidualPxPy{};
		TH2F					*h_StdTrkProton_ResidualPxPz{};
		TH2F					*h_StdTrkProton_ResidualPyPz{};
		TH2F					*h_StdTrkProton_ResidualPxE{};
		TH2F					*h_StdTrkProton_ResidualPyE{};
		TH2F					*h_StdTrkProton_ResidualPzE{};
		TH1F					*h_RFTrkProton_NormalizedResidualPx{};
		TH1F					*h_RFTrkProton_NormalizedResidualPxPy{};
		TH1F					*h_RFTrkProton_NormalizedResidualPy{};
		TH1F					*h_RFTrkProton_NormalizedResidualPxPz{};
		TH1F					*h_RFTrkProton_NormalizedResidualPyPz{};
		TH1F					*h_RFTrkProton_NormalizedResidualPz{};
		TH1F					*h_RFTrkProton_NormalizedResidualPxE{};
		TH1F					*h_RFTrkProton_NormalizedResidualPyE{};
		TH1F					*h_RFTrkProton_NormalizedResidualPzE{};
		TH1F					*h_RFTrkProton_NormalizedResidualE{};
		TH1F					*h_RFTrkProton_NormalizedResidualTheta{};
		TH1F					*h_RFTrkProton_NormalizedResidualPhi{};
		TH2F					*h_RFTrkProton_ResidualPxPy{};
		TH2F					*h_RFTrkProton_ResidualPxPz{};
		TH2F					*h_RFTrkProton_ResidualPyPz{};
		TH2F					*h_RFTrkProton_ResidualPxE{};
		TH2F					*h_RFTrkProton_ResidualPyE{};
		TH2F					*h_RFTrkProton_ResidualPzE{};
		TH1F					*h_TrueKaonLinkWeight{};
		TH1F					*h_StdTrkKaon_NormalizedResidualPx{};
		TH1F					*h_StdTrkKaon_NormalizedResidualPxPy{};
		TH1F					*h_StdTrkKaon_NormalizedResidualPy{};
		TH1F					*h_StdTrkKaon_NormalizedResidualPxPz{};
		TH1F					*h_StdTrkKaon_NormalizedResidualPyPz{};
		TH1F					*h_StdTrkKaon_NormalizedResidualPz{};
		TH1F					*h_StdTrkKaon_NormalizedResidualPxE{};
		TH1F					*h_StdTrkKaon_NormalizedResidualPyE{};
		TH1F					*h_StdTrkKaon_NormalizedResidualPzE{};
		TH1F					*h_StdTrkKaon_NormalizedResidualE{};
		TH1F					*h_StdTrkKaon_NormalizedResidualTheta{};
		TH1F					*h_StdTrkKaon_NormalizedResidualPhi{};
		TH2F					*h_StdTrkKaon_ResidualPxPy{};
		TH2F					*h_StdTrkKaon_ResidualPxPz{};
		TH2F					*h_StdTrkKaon_ResidualPyPz{};
		TH2F					*h_StdTrkKaon_ResidualPxE{};
		TH2F					*h_StdTrkKaon_ResidualPyE{};
		TH2F					*h_StdTrkKaon_ResidualPzE{};
		TH1F					*h_RFTrkKaon_NormalizedResidualPx{};
		TH1F					*h_RFTrkKaon_NormalizedResidualPxPy{};
		TH1F					*h_RFTrkKaon_NormalizedResidualPy{};
		TH1F					*h_RFTrkKaon_NormalizedResidualPxPz{};
		TH1F					*h_RFTrkKaon_NormalizedResidualPyPz{};
		TH1F					*h_RFTrkKaon_NormalizedResidualPz{};
		TH1F					*h_RFTrkKaon_NormalizedResidualPxE{};
		TH1F					*h_RFTrkKaon_NormalizedResidualPyE{};
		TH1F					*h_RFTrkKaon_NormalizedResidualPzE{};
		TH1F					*h_RFTrkKaon_NormalizedResidualE{};
		TH1F					*h_RFTrkKaon_NormalizedResidualTheta{};
		TH1F					*h_RFTrkKaon_NormalizedResidualPhi{};
		TH2F					*h_RFTrkKaon_ResidualPxPy{};
		TH2F					*h_RFTrkKaon_ResidualPxPz{};
		TH2F					*h_RFTrkKaon_ResidualPyPz{};
		TH2F					*h_RFTrkKaon_ResidualPxE{};
		TH2F					*h_RFTrkKaon_ResidualPyE{};
		TH2F					*h_RFTrkKaon_ResidualPzE{};
		TH1F					*h_StdTrk_NormalizedResidualPx{};
		TH1F					*h_StdTrk_NormalizedResidualPxPy{};
		TH1F					*h_StdTrk_NormalizedResidualPy{};
		TH1F					*h_StdTrk_NormalizedResidualPxPz{};
		TH1F					*h_StdTrk_NormalizedResidualPyPz{};
		TH1F					*h_StdTrk_NormalizedResidualPz{};
		TH1F					*h_StdTrk_NormalizedResidualPxE{};
		TH1F					*h_StdTrk_NormalizedResidualPyE{};
		TH1F					*h_StdTrk_NormalizedResidualPzE{};
		TH1F					*h_StdTrk_NormalizedResidualE{};
		TH1F					*h_StdTrk_NormalizedResidualTheta{};
		TH1F					*h_StdTrk_NormalizedResidualPhi{};
		TH1F					*h_RFTrk_NormalizedResidualPx{};
		TH1F					*h_RFTrk_NormalizedResidualPxPy{};
		TH1F					*h_RFTrk_NormalizedResidualPy{};
		TH1F					*h_RFTrk_NormalizedResidualPxPz{};
		TH1F					*h_RFTrk_NormalizedResidualPyPz{};
		TH1F					*h_RFTrk_NormalizedResidualPz{};
		TH1F					*h_RFTrk_NormalizedResidualPxE{};
		TH1F					*h_RFTrk_NormalizedResidualPyE{};
		TH1F					*h_RFTrk_NormalizedResidualPzE{};
		TH1F					*h_RFTrk_NormalizedResidualE{};
		TH1F					*h_RFTrk_NormalizedResidualTheta{};
		TH1F					*h_RFTrk_NormalizedResidualPhi{};

};

#endif
