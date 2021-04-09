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
		int getTruthTrkID( EVENT::LCEvent *pLCEvent , EVENT::Track* inputTrk );
		int getTrackIndex( EVENT::LCEvent *pLCEvent , EVENT::Track* inputTrk );
		TLorentzVector getMCPFourMomentum( EVENT::LCEvent *pLCEvent, EVENT::Track* inputTrk );
		TLorentzVector getTrackFourMomentum( EVENT::Track* inputTrk , double trackMass );
		std::vector<float> getPFOResidual( TLorentzVector pfoFourMomentum , TLorentzVector mcpFourMomentum );
		std::vector<float> UpdateChargedPFOCovMat( EVENT::Track* inputTrk , float trackMass );
		virtual void updatePFO( EVENT::ReconstructedParticle* inputPFO , ReconstructedParticleImpl* outputPFO , std::vector<Track*> outputPFOtrkvec , TLorentzVector pfoFourMomentum , std::vector<float> pfoCovMat );
		std::vector<float> getAngularUncertainties( TLorentzVector pfoFourMomentum , std::vector<float> pfoCovMat );
//		virtual void InitializeHistogram( TH1F *histogram , std::int scale , std::int color , std::int lineWidth , std::int markerSize , std::int markerStyle );
		virtual void InitializeHistogram( TH1F *histogram , int scale , int color , int lineWidth , int markerSize , int markerStyle );
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
		bool					m_fillRootTree = true;

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
		float					m_pion_mass;
		float					m_proton_mass;
		float					m_kaon_mass;
		TFile					*m_pTFile{};
		TTree					*m_pTTree{};
	        TTree					*m_pTTree_1trk{};
		TTree					*m_pTTree_2trk{};
		TTree					*m_pTTree_ntrk{};
		TDirectory				*m_Histograms{};
		TDirectory				*m_OldPFOs_1Trk{};
		TDirectory				*m_NewPFOs_1Trk{};
		TDirectory				*m_NewPFOs_2Trk{};
		TDirectory				*m_NewPFOs_nTrk{};
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

};

#endif
