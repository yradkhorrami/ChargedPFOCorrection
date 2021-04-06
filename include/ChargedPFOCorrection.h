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
//		MCParticle* getTruthTrkID(EVENT::LCEvent *pLCEvent, EVENT::Track* inputTrk);
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
		IntVector				m_pfoCharge{};
		IntVector				m_nTracksOfPFO{};
//		IntVector				m_nClustersOfPFO{};
//		IntVector				m_nSubDetectorHits{};
		float					m_Bfield;
		double					c;
		double					mm2m;
		double					eV2GeV;
		double					eB;
		double					proton_mass;
		double					kaon_mass;
		double					pion_mass;
		float					m_MinWeightTrackMCTruthLink;
		TFile					*m_pTFile{};
		TTree					*m_pTTree{};
	        TTree					*m_pTTree_1trk{};
		TTree					*m_pTTree_2trk{};
		TTree					*m_pTTree_ntrk{};
		TDirectory				*m_Histograms;
		TH2I					*h_nClusters_nTracks{};
		TH2I					*h_pfoCharge_nTracks{};
		TH2F					*h_InnermostRadiusHit_Neutral{};
		TH2F					*h_InnermostRadiusHit_Charged{};
		TH2I					*h_FirstSubDet_Charged{};

};

#endif
