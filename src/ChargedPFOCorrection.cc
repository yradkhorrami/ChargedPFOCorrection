#include "ChargedPFOCorrection.h"
#include <iostream>
#include <EVENT/LCCollection.h>
#include "EVENT/LCCollection.h"
#include "IMPL/LCCollectionVec.h"
#include "EVENT/MCParticle.h"
#include "EVENT/ReconstructedParticle.h"
#include <IMPL/ReconstructedParticleImpl.h>
#include "IMPL/ParticleIDImpl.h"
#include "UTIL/PIDHandler.h"
#include "marlin/VerbosityLevels.h"
#include <GeometryUtil.h>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1I.h"
#include "TH2I.h"
#include "TTree.h"

using namespace lcio ;
using namespace marlin ;

ChargedPFOCorrection aChargedPFOCorrection;

ChargedPFOCorrection::ChargedPFOCorrection() :

	Processor("ChargedPFOCorrection"),
	m_nRun(0),
	m_nEvt(0),
	m_nRunSum(0),
	m_nEvtSum(0),
	m_Bfield(0.f),
	c(0.),
	mm2m(0.),
	eV2GeV(0.),
	eB(0.),
	proton_mass(0.),
	kaon_mass(0.),
	pion_mass(0.)
{
	_description = "ChargedPFOCorrection creates new RECONSTRUCTEDPARTICLE collection that PFOs are updated using tracks refitted with true mass for protons and kaons";

	registerInputCollection(	LCIO::MCPARTICLE,
					"MCParticleCollection" ,
					"Name of the MCParticle collection"  ,
					m_mcParticleCollection,
					std::string("MCParticle")
				);

	registerInputCollection(	LCIO::RECONSTRUCTEDPARTICLE,
					"PfoCollection",
					"Name of input pfo collection",
					m_inputPfoCollection,
					std::string("PandoraPFOs")
				);

	registerInputCollection(	LCIO::TRACK,
					"MarlinTrkTracksCollection" ,
					"Name of the MarlinTrkTracks collection"  ,
					m_MarlinTrkTracks,
					std::string("MarlinTrkTracks")
				);

	registerInputCollection(	LCIO::TRACK,
					"MarlinTrkTracksCollectionKaon" ,
					"Name of the MarlinTrkTracks collection"  ,
					m_MarlinTrkTracksKAON,
					std::string("MarlinTrkTracksKaon")
				);

	registerInputCollection(	LCIO::TRACK,
					"MarlinTrkTracksCollectionProton" ,
					"Name of the MarlinTrkTracks collection"  ,
					m_MarlinTrkTracksPROTON,
					std::string("MarlinTrkTracksProton")
				);

	registerInputCollection(	LCIO::LCRELATION,
					"TrackMCTruthLinkCollection",
					"Name of input TrackMCTruthLink Collection",
					m_TrackMCTruthLinkCollection,
					std::string("MarlinTrkTracksMCTruthLink")
				);

	registerInputCollection(	LCIO::LCRELATION,
					"MCTruthTrackLinkCollection",
					"Name of input MCTruthTrackLink Collection",
					m_MCTruthTrackLinkCollection,
					std::string("MCTruthMarlinTrkTracksLink")
				);

	registerOutputCollection(	LCIO::RECONSTRUCTEDPARTICLE,
					"updatedPfoCollection",
					"Name of output pfo collection",
					m_outputPfoCollection,
					std::string("updatedPfoCollection")
				);

	registerProcessorParameter(	"updatePFOwithOneTrack",
					"Update PFOs with one track.",
					m_updatePFOwithOneTrack,
					bool(true)
				);

	registerProcessorParameter(	"updatePFOwithTwoTrack",
					"Update PFOs with two track.",
					m_updatePFOwithTwoTrack,
					bool(true)
				);

	registerProcessorParameter(	"updatePFOwithMoreTrack",
					"Update PFOs with More than Two track.",
					m_updatePFOwithMoreTrack,
					bool(true)
				);

	registerProcessorParameter(	"updatePFOwithPionTrack",
					"Update PFOs with Pion track.",
					m_updatePFOwithPionTrack,
					bool(true)
				);

	registerProcessorParameter(	"MinWeightTrackMCTruthLink" ,
					"Minimum acceptable weight for Track <--> MCParticle Link"  ,
					m_MinWeightTrackMCTruthLink ,
					float(0.9f)
				);

	registerProcessorParameter(	"RootFile",
	                                "Name of the output root file",
					m_rootFile,
					std::string("Output.root")
				);

}

void ChargedPFOCorrection::init()
{
	streamlog_out(DEBUG) << "	init called  " << std::endl;
//	m_Bfield = MarlinUtil::getBzAtOrigin();
	m_Bfield = 3.5;
	c = 2.99792458e8;
	mm2m = 1e-3;
	eV2GeV = 1e-9;
	eB = m_Bfield * c * mm2m * eV2GeV;
	printParameters();
	streamlog_out(DEBUG) << "	BField =  "<< m_Bfield << " Tesla" << std::endl ;
	m_nRun = 0 ;
	m_nEvt = 0 ;
	m_nRunSum = 0;
	m_nEvtSum = 0;
	proton_mass = 0.938272088;
	kaon_mass = 0.493677;
	pion_mass = 0.13957018;

	m_pTFile = new TFile(m_rootFile.c_str(), "recreate");

	m_pTTree = new TTree("PFOswithRFT", "PFOswithRFT");
	m_pTTree->SetDirectory(m_pTFile);
	m_pTTree->Branch("run", &m_nRun, "run/I");
	m_pTTree->Branch("event", &m_nEvt, "event/I");
	m_Histograms = m_pTFile->mkdir("Histograms");
	h_nClusters_nTracks = new TH2I("All PFOs", "; n_{Clusters}; n_{Tracks}", 11, -0.5, 10.5, 11, -0.5, 10.5);
	h_pfoCharge_nTracks = new TH2I("Neutral PFOs", "; PFO Charge; n_{Tracks}", 5, -2.5, 2.5, 11, -0.5, 10.5);
	h_InnermostRadiusHit_Neutral = new TH2F("Neutral PFOs with 2 tracks", "; r_{innermost hit}^{Track 1} [mm]; r_{innermost hit}^{Track 2} [mm]", 180, 0.0, 1800.0, 180, 0.0, 1800.0);
	h_InnermostRadiusHit_Charged = new TH2F("Charged PFOs with 2 tracks", "; r_{innermost hit}^{Track 1} [mm]; r_{innermost hit}^{Track 2} [mm]", 180, 0.0, 1800.0, 180, 0.0, 1800.0);
	h_FirstSubDet_Charged = new TH2I("First SubDet of Charged PFOs with 2 tracks", "; First SubDet^{Track 1}; First SubDet^{Track 2}", 11, -0.5, 10.5, 11, -0.5, 10.5);

}

void ChargedPFOCorrection::Clear()
{

}

void ChargedPFOCorrection::processRunHeader()
{
	m_nRun = 0;
	m_nEvt = 0;
	++m_nRunSum;
}

void ChargedPFOCorrection::processEvent( EVENT::LCEvent *pLCEvent )
{

	m_nRun = pLCEvent->getRunNumber();
	m_nEvt = pLCEvent->getEventNumber();
	++m_nEvtSum;

	LCCollection *inputPfoCollection{};
	LCCollection *MarlinTrkTracks{};
	LCCollection *MarlinTrkTracksKAON{};
	LCCollection *MarlinTrkTracksPROTON{};
	LCCollection *MCParticleCollection{};
	LCCollectionVec *outputPfoCollection{};
	int n_PFO = -1;
	int n_TRK = -1;
	int n_TRKp = -1;
	int n_TRKk = -1;
	int n_MCP = -1;
	this->Clear();
	streamlog_out(MESSAGE) << "" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////	Processing event 	" << m_nEvt << "	////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;

        try
        {
		inputPfoCollection = pLCEvent->getCollection(m_inputPfoCollection);
		MarlinTrkTracks = pLCEvent->getCollection(m_MarlinTrkTracks);
		MarlinTrkTracksKAON = pLCEvent->getCollection(m_MarlinTrkTracksKAON);
		MarlinTrkTracksPROTON = pLCEvent->getCollection(m_MarlinTrkTracksPROTON);
		MCParticleCollection = pLCEvent->getCollection(m_mcParticleCollection);
		outputPfoCollection = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

		n_PFO = inputPfoCollection->getNumberOfElements();
		n_TRK = MarlinTrkTracks->getNumberOfElements();
		n_TRKk = MarlinTrkTracksKAON->getNumberOfElements();
		n_TRKp = MarlinTrkTracksPROTON->getNumberOfElements();
		n_MCP = MCParticleCollection->getNumberOfElements();
		if ( n_PFO == -1 ) streamlog_out(DEBUG4) << "	Input PFO collection (" << m_inputPfoCollection << ") has no element (PFO) " << std::endl;
		if ( n_TRK == -1 ) streamlog_out(DEBUG4) << "	Input TRACK collection (" << m_MarlinTrkTracks << ") has no element (Track) " << std::endl;
		if ( n_TRK != n_TRKk ) streamlog_out(DEBUG4) << "	Input TRACK collection (" << m_MarlinTrkTracksKAON << ") has in-equal number of elements with Main Track Collection (" << m_MarlinTrkTracks << ")!" << std::endl;
		if ( n_TRK != n_TRKp ) streamlog_out(DEBUG4) << "	Input TRACK collection (" << MarlinTrkTracksPROTON << ") has in-equal number of elements with Main Track Collection (" << m_MarlinTrkTracks << ")!" << std::endl;

		streamlog_out(DEBUG4) << "	Total Number of PFOs: " << n_PFO << std::endl;
		streamlog_out(DEBUG4) << "	Total Number of Tracks: " << n_TRK << std::endl;
		streamlog_out(DEBUG4) << "	Total Number of KaonTracks: " << n_TRKk << std::endl;
		streamlog_out(DEBUG4) << "	Total Number of ProtonTracks: " << n_TRKp << std::endl;
		streamlog_out(DEBUG4) << "	Total Number of MCParticles: " << n_MCP << std::endl;

		double trackMass = 0.0;
		int linkWeight = 100 * m_MinWeightTrackMCTruthLink;

		for (int i_pfo = 0; i_pfo < n_PFO ; ++i_pfo)
		{
			ReconstructedParticle* inputPFO = dynamic_cast<ReconstructedParticle*>( inputPfoCollection->getElementAt( i_pfo ) );
			ReconstructedParticleImpl* outputPFO = new ReconstructedParticleImpl;
			const EVENT::TrackVec& inputPFOtrkvec = inputPFO->getTracks();
			TLorentzVector oldpfoFourMomentum( inputPFO->getMomentum()[ 0 ] , inputPFO->getMomentum()[ 1 ] , inputPFO->getMomentum()[ 2 ] , inputPFO->getEnergy() );
			std::vector<Track*> outputPFOtrkvec;
			int nTRKsofPFO = inputPFOtrkvec.size();
			int nClustersofPFO = (inputPFO->getClusters()).size();
			h_nClusters_nTracks->Fill( nClustersofPFO , nTRKsofPFO );
			Track *refittedTrack = NULL;
			double maxweightTRKtoMCP = 0.;
			bool m_updatePFO = true;
			bool foundLinkedMCP = true;
			float pfoMass = 0.0;
			TLorentzVector pfoFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
			TLorentzVector mcpFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
			std::vector<float> oldResiduals( 6 , 0.0 );
			std::vector<float> newResiduals( 6 , 0.0 );
			std::vector<float> newPFOCovMat( 10, 0.0 );
			std::vector<float> oldPFOCovMat = inputPFO->getCovMatrix();
			std::vector<float> trkCovMat( 10, 0.0 );
			std::vector<float> RadiusOfInnermostHit( nTRKsofPFO , 0.0 );

			streamlog_out(DEBUG3) << " " << std::endl;
			streamlog_out(DEBUG3) << "	**********************************************************************************" << std::endl;
			if ( inputPFO->getCharge() == 0.0 )
			{
				streamlog_out(DEBUG3) << "	Processing a neutral PFO at index: 	" << i_pfo << " 	with 	" << nTRKsofPFO << " 	track(s)" << std::endl;
			}
			else
			{
				streamlog_out(DEBUG3) << "	Processing a charged PFO at index: 	" << i_pfo << " 	with 	" << nTRKsofPFO << " 	track(s)" << std::endl;
			}
			h_pfoCharge_nTracks->Fill( inputPFO->getCharge() , nTRKsofPFO );

			if ( nTRKsofPFO == 0 )
			{
				m_updatePFO = false;
				foundLinkedMCP = false;
			}
			else if ( nTRKsofPFO == 1 )
			{
				Track *inputTrk = (Track*)inputPFOtrkvec.at(0);
				int TrackID = this->getTruthTrkID( pLCEvent, inputTrk );
				int TrackIndex = this->getTrackIndex( pLCEvent, inputTrk );
				if ( TrackID != 0 )
				{
					streamlog_out(DEBUG1) << "	A MCParticle linked to Track with weight higher than " << linkWeight << "% in both direction" << std::endl;
					mcpFourMomentum = this->getMCPFourMomentum( pLCEvent, inputTrk );
					foundLinkedMCP = true;
				}
				else
				{
					streamlog_out(DEBUG1) << "	Couldn't find a MCParticle linked to Track with weight higher than " << linkWeight << "% in both direction!" << std::endl;
					foundLinkedMCP = false;
					m_updatePFO = false;
				}
				streamlog_out(DEBUG2) << "	TrackID: 	" << TrackID << std::endl;
				streamlog_out(DEBUG2) << "	TrackIndex: 	" << TrackIndex << std::endl;

				if ( abs( TrackID ) == 2212 )
				{
					refittedTrack = dynamic_cast<EVENT::Track*>( MarlinTrkTracksPROTON->getElementAt( TrackIndex ) );
					outputPFOtrkvec.push_back( refittedTrack );
					trackMass = proton_mass;
					m_updatePFO = true;
					streamlog_out(DEBUG2) << "	Default track is replace with track refitted with proton mass" << std::endl;
				}
				else if ( abs( TrackID ) == 321 )
				{
					refittedTrack = dynamic_cast<EVENT::Track*>( MarlinTrkTracksKAON->getElementAt( TrackIndex ) );
					outputPFOtrkvec.push_back( refittedTrack );
					trackMass = kaon_mass;
					m_updatePFO = true;
					streamlog_out(DEBUG2) << "	Default track is replace with track refitted with kaon mass" << std::endl;
				}
				else
				{
					refittedTrack = dynamic_cast<EVENT::Track*>( MarlinTrkTracks->getElementAt( TrackIndex ) );
					outputPFOtrkvec.push_back( refittedTrack );
					trackMass = pion_mass;
					m_updatePFO = m_updatePFOwithPionTrack;
					streamlog_out(DEBUG2) << "	Default track is used for updating PFO" << std::endl;
				}
				pfoFourMomentum = this->getTrackFourMomentum( refittedTrack , trackMass );
				newPFOCovMat = this->UpdateChargedPFOCovMat( refittedTrack , trackMass );
			}
			else if ( nTRKsofPFO == 2 )
			{
/*
				const EVENT::IntVec& SubDetectorHitNumbers_trk1 = ( (Track*)inputPFOtrkvec.at( 0 ) )->getSubdetectorHitNumbers();
				const EVENT::IntVec& SubDetectorHitNumbers_trk2 = ( (Track*)inputPFOtrkvec.at( 1 ) )->getSubdetectorHitNumbers();
				for ( int i_trk = 0 ; i_trk < nTRKsofPFO ; ++i_trk )
				{
					Track *inputTrk = (Track*)inputPFOtrkvec.at( i_trk );
					const EVENT::IntVec& SubDetectorHitNumbers = inputTrk->getSubdetectorHitNumbers();
					RadiusOfInnermostHit[ i_trk ] = inputTrk->getRadiusOfInnermostHit();
					streamlog_out(DEBUG3) << "	Track at index 	" << i_trk << " 	has innermost hit with the radius of: 	" << RadiusOfInnermostHit[ i_trk ] << std::endl;
					streamlog_out(DEBUG3) << "	SubDetectorHitNumbers: [ " << SubDetectorHitNumbers[ 0 ];
					for ( int i_subDet = 1 ; i_subDet < SubDetectorHitNumbers.size() ; ++i_subDet )
					{
						streamlog_out(DEBUG3) << "," << SubDetectorHitNumbers[ i_subDet ];
					}
					streamlog_out(DEBUG3) << " ]" << std::endl;
				}
*/
				if ( inputPFO->getCharge() == 0.0 )
				{
					h_InnermostRadiusHit_Neutral->Fill( RadiusOfInnermostHit[ 0 ] , RadiusOfInnermostHit[ 0 ] );
					for ( int i_trk = 0 ; i_trk < nTRKsofPFO ; ++i_trk )
					{
						streamlog_out(DEBUG2) << "	---------------------" << std::endl;
						streamlog_out(DEBUG2) << "	Adding Track " << i_trk << " to PFO " << std::endl;
						streamlog_out(DEBUG2) << "	---------------------" << std::endl;
						Track *inputTrk = (Track*)inputPFOtrkvec.at( i_trk );
						int TrackID = this->getTruthTrkID( pLCEvent, inputTrk );
						int TrackIndex = this->getTrackIndex( pLCEvent, inputTrk );
						if ( TrackID != 0 )
						{
							streamlog_out(DEBUG1) << "	A MCParticle linked to Track[" << i_trk << "] with weight higher than " << linkWeight << "% in both direction" << std::endl;
							mcpFourMomentum += this->getMCPFourMomentum( pLCEvent, inputTrk );
							if ( foundLinkedMCP ) foundLinkedMCP = true;
						}
						else
						{
							streamlog_out(DEBUG1) << "	Couldn't find a MCParticle linked to Track[" << i_trk << "] with weight higher than " << linkWeight << "% in both direction!" << std::endl;
							if ( foundLinkedMCP ) foundLinkedMCP = false;
						}
						streamlog_out(DEBUG2) << "	TrackID: 	" << TrackID << std::endl;
						streamlog_out(DEBUG2) << "	TrackIndex: 	" << TrackIndex << std::endl;

						if ( abs( TrackID ) == 2212 )
						{
							refittedTrack = dynamic_cast<EVENT::Track*>( MarlinTrkTracksPROTON->getElementAt( TrackIndex ) );
							outputPFOtrkvec.push_back( refittedTrack );
							trackMass = proton_mass;
							m_updatePFO = true;
							streamlog_out(DEBUG2) << "	Default track is replace with track refitted with proton mass" << std::endl;
						}
						else if ( abs( TrackID ) == 321 )
						{
							refittedTrack = dynamic_cast<EVENT::Track*>( MarlinTrkTracksKAON->getElementAt( TrackIndex ) );
							outputPFOtrkvec.push_back( refittedTrack );
							trackMass = kaon_mass;
							m_updatePFO = true;
							streamlog_out(DEBUG2) << "	Default track is replace with track refitted with kaon mass" << std::endl;
						}
						else
						{
							refittedTrack = dynamic_cast<EVENT::Track*>( MarlinTrkTracks->getElementAt( TrackIndex ) );
							outputPFOtrkvec.push_back( refittedTrack );
							trackMass = pion_mass;
							m_updatePFO = m_updatePFOwithPionTrack;
							streamlog_out(DEBUG2) << "	Default track is used for updating PFO" << std::endl;
						}
						pfoFourMomentum += this->getTrackFourMomentum( refittedTrack , trackMass );
						trkCovMat = this->UpdateChargedPFOCovMat( refittedTrack , trackMass );
						for ( int i_element = 0 ; i_element < 10 ; ++i_element )
						{
							newPFOCovMat[ i_element ] += trkCovMat[ i_element ];
						}

					}
				}
				else
				{
/*
					h_InnermostRadiusHit_Charged->Fill( RadiusOfInnermostHit[ 0 ] , RadiusOfInnermostHit[ 1 ] );
					int firstSubDet_trk1 = 41;
					int firstSubDet_trk2 = 41;
					for ( int i_subDet = 0 ; i_subDet < SubDetectorHitNumbers_trk1.size() ; ++i_subDet )
					{
						if ( SubDetectorHitNumbers_trk1[ i_subDet ] != 0 )
						{
							firstSubDet_trk1 = i_subDet;
							continue;
						}
					}
					for ( int i_subDet = 0 ; i_subDet < SubDetectorHitNumbers_trk2.size() ; ++i_subDet )
					{
						if ( SubDetectorHitNumbers_trk2[ i_subDet ] != 0 )
						{
							firstSubDet_trk2 = i_subDet;
							continue;
						}
					}
					if ( RadiusOfInnermostHit[ 0 ] > 200 )// RadiusOfInnermostHit[ 1 ] )
					{
						h_FirstSubDet_Charged->Fill( firstSubDet_trk1 , firstSubDet_trk2 );
					}
*/
				}
			}
			if( m_updatePFO )
			{
				this->updatePFO( inputPFO , outputPFO , outputPFOtrkvec , pfoFourMomentum , newPFOCovMat );
			}
			else
			{
				this->updatePFO( inputPFO , outputPFO , inputPFOtrkvec , oldpfoFourMomentum , oldPFOCovMat );
			}

			streamlog_out(DEBUG3) << "	------------------------------" << std::endl;
			streamlog_out(DEBUG3) << "	Getting Residuals for inputPFO" << std::endl;
			if ( foundLinkedMCP ) oldResiduals = this->getPFOResidual( oldpfoFourMomentum , mcpFourMomentum );
			streamlog_out(DEBUG3) << "	-------------------------------" << std::endl;
			streamlog_out(DEBUG3) << "	Getting Residuals for outputPFO" << std::endl;
			TLorentzVector newpfoFourMomentum( outputPFO->getMomentum()[ 0 ] , outputPFO->getMomentum()[ 1 ] , outputPFO->getMomentum()[ 2 ] , outputPFO->getEnergy() );
			if ( foundLinkedMCP ) newResiduals = this->getPFOResidual( newpfoFourMomentum , mcpFourMomentum );
			outputPfoCollection->addElement( outputPFO );
		}


		pLCEvent->addCollection( outputPfoCollection , m_outputPfoCollection );
		m_pTTree->Fill();
	}
	catch(DataNotAvailableException &e)
        {
        	streamlog_out(MESSAGE) << "	Input collection not found in event " << m_nEvt << std::endl;
        }

}

int ChargedPFOCorrection::getTruthTrkID( EVENT::LCEvent *pLCEvent , EVENT::Track* inputTrk )
{
	std::vector<float> TRKMCPlink;
	int TrackID = 0;
	MCParticle *linkedMCP{};
	try
	{
		LCRelationNavigator TrackMCParticleNav( pLCEvent->getCollection( m_TrackMCTruthLinkCollection ) );
		LCRelationNavigator MCParticleTrackNav( pLCEvent->getCollection( m_MCTruthTrackLinkCollection ) );
	}
	catch (DataNotAvailableException &e)
	{
		streamlog_out(WARNING) << "Could not find the LCRelationNavigator for Track <-> MCParticle" << std::endl;
		return TrackID;
	}
	LCRelationNavigator TrackMCParticleNav( pLCEvent->getCollection( m_TrackMCTruthLinkCollection ) );
	const EVENT::LCObjectVec& mcpvec = TrackMCParticleNav.getRelatedToObjects(inputTrk);
	const EVENT::FloatVec&  mcpweightvec = TrackMCParticleNav.getRelatedToWeights(inputTrk);
	double maxweightTRKtoMCP = 0.;
	int iTRKtoMCPmax = -1;
	for ( unsigned int i_mcp = 0; i_mcp < mcpvec.size(); i_mcp++ )
	{
		double mcp_weight = mcpweightvec.at(i_mcp);
		MCParticle *testMCP = (MCParticle *) mcpvec.at(i_mcp);
		if ( mcp_weight > maxweightTRKtoMCP && mcp_weight >= m_MinWeightTrackMCTruthLink )
		{
			maxweightTRKtoMCP = mcp_weight;
			iTRKtoMCPmax = i_mcp;
			streamlog_out(DEBUG0) << "	MCParticle at index: " << i_mcp << " has PDG: " << testMCP->getPDG() << " and MCParticle to Track link weight is " << mcp_weight << std::endl;
		}
	}
	if ( iTRKtoMCPmax != -1 )
	{
		MCParticle *testMCP = (MCParticle *) mcpvec.at(iTRKtoMCPmax);
		streamlog_out(DEBUG0) << "	MCParticle at index: " << iTRKtoMCPmax << " with PDG code: " << testMCP->getPDG() << " has highest link weight to input track" << std::endl;
		LCRelationNavigator MCParticleTrackNav( pLCEvent->getCollection( m_MCTruthTrackLinkCollection ) );
		const EVENT::LCObjectVec& trkvec = MCParticleTrackNav.getRelatedToObjects(testMCP);
		const EVENT::FloatVec&  trkweightvec = MCParticleTrackNav.getRelatedToWeights(testMCP);
		double maxweightMCPtoTRK = 0.;
		int iMCPtoTRKmax = -1;
		for ( unsigned int i_trk = 0; i_trk < trkvec.size(); i_trk++ )
		{
			double track_weight = trkweightvec.at( i_trk );
			if ( track_weight > maxweightMCPtoTRK && track_weight >= m_MinWeightTrackMCTruthLink )
			{
				maxweightMCPtoTRK = track_weight;
				iMCPtoTRKmax = i_trk;
				streamlog_out(DEBUG0) << "	Track at index: " << i_trk << " has highest link weight to MCParticle = " << track_weight << std::endl;
			}
		}
		if ( iMCPtoTRKmax != -1 )
		{
			Track *linkedTrack = (Track *) trkvec.at( iMCPtoTRKmax );
			if ( linkedTrack == inputTrk )
			{
				streamlog_out(DEBUG0) << "	Track linked to MCParticle with highest weight" << std::endl;
				linkedMCP = (MCParticle *) mcpvec.at(iTRKtoMCPmax);
				TrackID = linkedMCP->getPDG();
			}
		}
	}
	return TrackID;
}

int ChargedPFOCorrection::getTrackIndex( EVENT::LCEvent *pLCEvent , EVENT::Track* inputTrk )
{
	LCCollection *MarlinTrkTracks{};
	try
	{
		MarlinTrkTracks = pLCEvent->getCollection( m_MarlinTrkTracks );
	}
	catch (DataNotAvailableException &e)
	{
		streamlog_out(WARNING) << "	Could not find the " << m_MarlinTrkTracks << " Collection" << std::endl;
		return -1;
	}
	unsigned int nTRKs = MarlinTrkTracks->getNumberOfElements();
	int track_index = -1;
	for (unsigned int i_trk = 0; i_trk < nTRKs;  ++i_trk )
	{
		Track* pionTrack = dynamic_cast<EVENT::Track*>( MarlinTrkTracks->getElementAt( i_trk ) );
		if ( pionTrack == inputTrk )
		{
			track_index = i_trk;
		}
	}
	if ( track_index == -1)
	{
		streamlog_out(DEBUG1) << "	Coudln't find track_index!!!  " << std::endl;
	}
	else
	{
		streamlog_out(DEBUG1) << "	Track index in " << m_MarlinTrkTracks << " collection is " << track_index << std::endl;
	}
	return track_index;
}

TLorentzVector ChargedPFOCorrection::getMCPFourMomentum( EVENT::LCEvent *pLCEvent , EVENT::Track* inputTrk )
{
	std::vector<float> TRKMCPlink;
	MCParticle *linkedMCP{};
	TLorentzVector mcpFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	try
	{
		LCRelationNavigator TrackMCParticleNav( pLCEvent->getCollection( m_TrackMCTruthLinkCollection ) );
		LCRelationNavigator MCParticleTrackNav( pLCEvent->getCollection( m_MCTruthTrackLinkCollection ) );
	}
	catch (DataNotAvailableException &e)
	{
		streamlog_out(WARNING) << "Could not find the LCRelationNavigator for Track <-> MCParticle" << std::endl;
		return mcpFourMomentum;
	}
	LCRelationNavigator TrackMCParticleNav( pLCEvent->getCollection( m_TrackMCTruthLinkCollection ) );
	const EVENT::LCObjectVec& mcpvec = TrackMCParticleNav.getRelatedToObjects(inputTrk);
	const EVENT::FloatVec&  mcpweightvec = TrackMCParticleNav.getRelatedToWeights(inputTrk);
	double maxweightTRKtoMCP = 0.;
	int iTRKtoMCPmax = -1;
	for ( unsigned int i_mcp = 0; i_mcp < mcpvec.size(); i_mcp++ )
	{
		double mcp_weight = mcpweightvec.at(i_mcp);
		MCParticle *testMCP = (MCParticle *) mcpvec.at(i_mcp);
		if ( mcp_weight > maxweightTRKtoMCP && mcp_weight >= m_MinWeightTrackMCTruthLink )
		{
			maxweightTRKtoMCP = mcp_weight;
			iTRKtoMCPmax = i_mcp;
			streamlog_out(DEBUG0) << "	MCParticle at index: " << i_mcp << " has PDG: " << testMCP->getPDG() << " and MCParticle to Track link weight is " << mcp_weight << std::endl;
		}
	}
	if ( iTRKtoMCPmax != -1 )
	{
		MCParticle *testMCP = (MCParticle *) mcpvec.at(iTRKtoMCPmax);
		streamlog_out(DEBUG0) << "	MCParticle at index: " << iTRKtoMCPmax << " with PDG code: " << testMCP->getPDG() << " has highest link weight to input track" << std::endl;
		LCRelationNavigator MCParticleTrackNav( pLCEvent->getCollection( m_MCTruthTrackLinkCollection ) );
		const EVENT::LCObjectVec& trkvec = MCParticleTrackNav.getRelatedToObjects(testMCP);
		const EVENT::FloatVec&  trkweightvec = MCParticleTrackNav.getRelatedToWeights(testMCP);
		double maxweightMCPtoTRK = 0.;
		int iMCPtoTRKmax = -1;
		for ( unsigned int i_trk = 0; i_trk < trkvec.size(); i_trk++ )
		{
			double track_weight = trkweightvec.at( i_trk );
			if ( track_weight > maxweightMCPtoTRK && track_weight >= m_MinWeightTrackMCTruthLink )
			{
				maxweightMCPtoTRK = track_weight;
				iMCPtoTRKmax = i_trk;
				streamlog_out(DEBUG0) << "	Track at index: " << i_trk << " has highest link weight to MCParticle = " << track_weight << std::endl;
			}
		}
		if ( iMCPtoTRKmax != -1 )
		{
			Track *linkedTrack = (Track *) trkvec.at( iMCPtoTRKmax );
			if ( linkedTrack == inputTrk )
			{
				streamlog_out(DEBUG0) << "	Track linked to MCParticle with highest weight" << std::endl;
				linkedMCP = (MCParticle *) mcpvec.at(iTRKtoMCPmax);
				mcpFourMomentum = TLorentzVector( linkedMCP->getMomentum()[0] , linkedMCP->getMomentum()[1] , linkedMCP->getMomentum()[2] , linkedMCP->getEnergy() );
			}
		}
	}
	return mcpFourMomentum;
}

TLorentzVector ChargedPFOCorrection::getTrackFourMomentum( EVENT::Track* inputTrk , double trackMass )
{
	streamlog_out(DEBUG1) << "	------------------------------------------------" << std::endl;
	streamlog_out(DEBUG1) << "	Calculating PFO 4-momentum from track parameters" << std::endl;
	streamlog_out(DEBUG1) << "	------------------------------------------------" << std::endl;
	double Phi = inputTrk->getPhi();
	double Omega = inputTrk->getOmega();
	double tanLambda = inputTrk->getTanLambda();
	streamlog_out(DEBUG0) << "	Track parameters obtained" << std::endl;
	double pT = eB / fabs( Omega );
	double px = pT * TMath::Cos( Phi );
	double py = pT * TMath::Sin( Phi );
	double pz = pT * tanLambda;
	double E = sqrt( pow( trackMass , 2 ) + px * px + py * py + pz * pz);
	streamlog_out(DEBUG0) << "	Track parameters is converted to (p,E)" << std::endl;
	TLorentzVector trackFourMomentum( px , py , pz , E );
	return trackFourMomentum;
}

std::vector<float> ChargedPFOCorrection::getPFOResidual( TLorentzVector pfoFourMomentum , TLorentzVector mcpFourMomentum )
{
	std::vector<float> pfoResidual;

	float pfoPx		= pfoFourMomentum.Px();
	float pfoPy		= pfoFourMomentum.Py();
	float pfoPz		= pfoFourMomentum.Pz();
	float pfoE		= pfoFourMomentum.E();
	TVector3 pfoPvec( pfoPx , pfoPy , pfoPz ); pfoPvec.SetMag(1.0);
	TVector3 pfoPTvec( pfoPx , pfoPy , 0.0 ); pfoPTvec.SetMag(1.0);
	float pfoTheta		= pfoPvec.Theta();
	float pfoPhi		= pfoPvec.Phi();
	streamlog_out(DEBUG0) << "	PFO 4-Momentum ( px , py , pz , E , Theta , Phi ) = ( " << pfoPx << " , " << pfoPy << " , " << pfoPz << " , " << pfoE << " , " << pfoTheta << " , " << pfoPhi << " )" << std::endl;

	float mcpPx		= mcpFourMomentum.Px();
	float mcpPy		= mcpFourMomentum.Py();
	float mcpPz		= mcpFourMomentum.Pz();
	float mcpE		= mcpFourMomentum.E();
	TVector3 mcpPvec( mcpPx , mcpPy , mcpPz ); mcpPvec.SetMag(1.0);
	TVector3 mcpPTvec( mcpPx , mcpPy , 0.0 ); mcpPTvec.SetMag(1.0);
	float mcpTheta		= mcpPvec.Theta();
	float mcpPhi		= mcpPvec.Phi();
	streamlog_out(DEBUG0) << "	MCP 4-Momentum ( px , py , pz , E , Theta , Phi ) = ( " << mcpPx << " , " << mcpPy << " , " << mcpPz << " , " << mcpE << " , " << mcpTheta << " , " << mcpPhi << " )" << std::endl;

	float ResidualPx	= pfoPx - mcpPx;
	float ResidualPy	= pfoPy - mcpPy;
	float ResidualPz	= pfoPz - mcpPz;
	float ResidualEnergy	= pfoE - mcpE;
	float ResidualTheta	= pfoTheta - mcpTheta;
	float ResidualPhi	= 0.0;
	if ( pfoPhi > mcpPhi )
	{
		ResidualPhi	= acos( pfoPTvec.Dot( mcpPTvec ) );
	}
	else
	{
		ResidualPhi	= -acos( pfoPTvec.Dot( mcpPTvec ) );
	}
	streamlog_out(DEBUG1) << "	Residuals	( deltaPx , deltaPy , deltaPz , deltaE , deltaTheta , deltaPhi ) = ( " << ResidualPx << " , " << ResidualPy << " , " << ResidualPz << " , " << ResidualEnergy << " , " << ResidualTheta << " , " << ResidualPhi << " )" << std::endl;

	pfoResidual.push_back( ResidualPx );
	pfoResidual.push_back( ResidualPy );
	pfoResidual.push_back( ResidualPz );
	pfoResidual.push_back( ResidualEnergy );
	pfoResidual.push_back( ResidualTheta );
	pfoResidual.push_back( ResidualPhi );

	return pfoResidual;

}

std::vector<float> ChargedPFOCorrection::UpdateChargedPFOCovMat( EVENT::Track* inputTrk , float trackMass )
{

//	Obtain covariance matrix on (px,py,pz,E) from the
//	covariance matrix on track parameters.
//
//	define the jacobian as the 3x4 matrix:
//
//
//
//			Dpx/DTanLambda		Dpy/DTanLambda		Dpz/DTanLambda		DE/DTanLambda
//
//	 J =		Dpx/DPhi		Dpy/DPhi		Dpz/DPhi		DE/DPhi
//
//			Dpx/DOmega		Dpy/DOmega		Dpz/DOmega		DE/DOmega
//
//
//
//			0			0			Pt			Pz.Pt/E
//
//	J =		-Py			Px			0			0
//
//			-Px/Omega		-Py/Omega		-Pz/Omega		-P2/E.Omega
//
//
//
//	Order in the covariance matrix on helix parameters:
//
//			TanLambda.TanLambda	TanLambda.phi		TanLambda.Omega
//
//	Cov =		phi.TanLambda		phi.phi			phi.Omega
//
//			Omega.TanLambda		Omega.phi		Omega.Omega
//
//
//

	streamlog_out(DEBUG1) << "	-------------------------------------------------------" << std::endl;
	streamlog_out(DEBUG1) << "	Calculating PFO CovMat from track parameters and CovMat" << std::endl;
	streamlog_out(DEBUG1) << "	-------------------------------------------------------" << std::endl;
	const int rows			= 3; // n rows jacobian
	const int columns		= 4; // n columns jacobian
	const int kspace_time_dim	= 4;

	TMatrixD covMatrixMomenta(kspace_time_dim,kspace_time_dim);
	std::vector<float> covP;

	double trackTanLambda		=	inputTrk->getTanLambda();
	double trackPhi			=	inputTrk->getPhi();
	double trackOmega		=	inputTrk->getOmega();
	std::vector<float> trackCovMat 	=	inputTrk->getCovMatrix();
	double trackPt			=	eB / fabs( trackOmega );
	double trackPx			= 	trackPt * TMath::Cos( trackPhi );
	double trackPy			= 	trackPt * TMath::Sin( trackPhi );
	double trackPz			= 	trackPt * trackTanLambda;
	double trackP			= 	std::sqrt( pow( trackPt , 2 ) + pow( trackPz , 2 ) );
	double trackE			= 	std::sqrt( pow( trackP , 2 ) + pow( trackMass , 2 ) );

	float SigmaPhi2			=	trackCovMat[  2 ];
	float SigmaPhiSigmaOmega	=	trackCovMat[  4 ];
	float SigmaOmega2		=	trackCovMat[  5 ];
	float SigmaTanLambdaSigmaPhi	=	trackCovMat[ 11 ];
	float SigmaTanLambdaSigmaOmega 	=	trackCovMat[ 12 ];
	float SigmaTanLambda2		=	trackCovMat[ 14 ];

	streamlog_out(DEBUG0) << "	Track parametrs and CovMat obtained" << std::endl;

//	Define array with jacobian matrix elements by rows
	double jacobian_by_rows[rows*columns] =
	{
		0			,	0			,	trackPt			,	trackPz * trackPt / trackE			,
		-trackPy		,	trackPx		,	0				,	0						,
		-trackPx / trackOmega	,	-trackPy / trackOmega	,	-trackPz / trackOmega		,	-trackP * trackP / ( trackE * trackOmega )
	};

	streamlog_out(DEBUG0) << "	Jacobian array formed by rows" << std::endl;

//	construct the Jacobian using previous array ("F" if filling by columns, "C", if filling by rows, $ROOTSYS/math/matrix/src/TMatrixT.cxx)
	TMatrixD jacobian(rows,columns, jacobian_by_rows, "C");
	streamlog_out(DEBUG0) << "	Jacobian array converted to Jacobian matrix" << std::endl;

//	track covariance matrix by rows
	double track_cov_matrix_by_rows[rows*rows] =
	{
		SigmaTanLambda2			,	SigmaTanLambdaSigmaPhi		,	SigmaTanLambdaSigmaOmega	,
		SigmaTanLambdaSigmaPhi		,	SigmaPhi2			,	SigmaPhiSigmaOmega	,
		SigmaTanLambdaSigmaOmega	,	SigmaPhiSigmaOmega		,	SigmaOmega2
	};
	streamlog_out(DEBUG0) << "	Track covariance matrix array formed by rows" << std::endl;

	TMatrixD covMatrix_track(rows,rows, track_cov_matrix_by_rows, "C");
	streamlog_out(DEBUG0) << "	Track CovMat array converted to track covariance matrix" << std::endl;

	covMatrixMomenta.Mult( TMatrixD( jacobian ,
					TMatrixD::kTransposeMult ,
					covMatrix_track) ,
					jacobian
					);
	streamlog_out(DEBUG0) << "	Track covariance matrix converted to FourMomentum Covariance matrix" << std::endl;

	covP.push_back( covMatrixMomenta(0,0) ); // x-x
	covP.push_back( covMatrixMomenta(1,0) ); // y-x
	covP.push_back( covMatrixMomenta(1,1) ); // y-y
	covP.push_back( covMatrixMomenta(2,0) ); // z-x
	covP.push_back( covMatrixMomenta(2,1) ); // z-y
	covP.push_back( covMatrixMomenta(2,2) ); // z-z
	covP.push_back( covMatrixMomenta(3,0) ); // e-x
	covP.push_back( covMatrixMomenta(3,1) ); // e-y
	covP.push_back( covMatrixMomenta(3,2) ); // e-z
	covP.push_back( covMatrixMomenta(3,3) ); // e-e
	streamlog_out(DEBUG0) << "	FourMomentum CovarianceMatrix Filled succesfully" << std::endl;

	return covP;

}

void ChargedPFOCorrection::updatePFO( EVENT::ReconstructedParticle* inputPFO , ReconstructedParticleImpl* outputPFO , std::vector<Track*> outputPFOtrkvec , TLorentzVector pfoFourMomentum , std::vector<float> pfoCovMat )
{
	double Momentum[3]{ pfoFourMomentum.Px() , pfoFourMomentum.Py() , pfoFourMomentum.Pz() };
	double Energy = pfoFourMomentum.E();
	double Mass = pfoFourMomentum.M();
	outputPFO->setType(inputPFO->getType());
	outputPFO->setMomentum( Momentum );
	outputPFO->setEnergy( Energy );
	outputPFO->setCovMatrix( pfoCovMat );
	outputPFO->setMass( Mass );
	outputPFO->setCharge( inputPFO->getCharge() );
	outputPFO->setReferencePoint( inputPFO->getReferencePoint() );
	for ( unsigned int j = 0 ; j < inputPFO->getParticleIDs().size() ; ++j )
	{
		ParticleIDImpl* inPID = dynamic_cast<ParticleIDImpl*>( inputPFO->getParticleIDs()[ j ] );
		ParticleIDImpl* outPID = new ParticleIDImpl;
		outPID->setType( inPID->getType() );
		outPID->setPDG( inPID->getPDG() );
		outPID->setLikelihood( inPID->getLikelihood() );
		outPID->setAlgorithmType( inPID->getAlgorithmType() ) ;
		for ( unsigned int k = 0 ; k < inPID->getParameters().size() ; ++k ) outPID->addParameter( inPID->getParameters()[ k ] );
		outputPFO->addParticleID( outPID );
	}
	outputPFO->setParticleIDUsed( inputPFO->getParticleIDUsed() );
	outputPFO->setGoodnessOfPID( inputPFO->getGoodnessOfPID() );
	for ( unsigned int j = 0 ; j < inputPFO->getParticles().size() ; ++j )
	{
		outputPFO->addParticle( inputPFO->getParticles()[ j ] );
	}
	for ( unsigned int j = 0 ; j < inputPFO->getClusters().size() ; ++j )
	{
		outputPFO->addCluster( inputPFO->getClusters()[ j ] );
	}
	for ( unsigned int j = 0 ; j < outputPFOtrkvec.size() ; ++j)
	{
		outputPFO->addTrack( outputPFOtrkvec[ j ] );
	}
	outputPFO->setStartVertex( inputPFO->getStartVertex() );
}

void ChargedPFOCorrection::check( EVENT::LCEvent *pLCEvent )
{
	LCCollection *inputPfoCollection{};
	LCCollection *outputPfoCollection{};
	try
	{
		inputPfoCollection = pLCEvent->getCollection(m_inputPfoCollection);
		outputPfoCollection = pLCEvent->getCollection(m_outputPfoCollection);
		int n_inputPFOs = inputPfoCollection->getNumberOfElements();
		int n_outputPFOs = outputPfoCollection->getNumberOfElements();
		streamlog_out(DEBUG4) << "	CHECK : processed events: " << m_nEvtSum << " (Number of inputPFOS: " << n_inputPFOs << " , Number of outputPFOs: " << n_outputPFOs <<")" << std::endl;
	}
	catch(DataNotAvailableException &e)
        {
          streamlog_out(MESSAGE) << "	Input/Output collection not found in event " << m_nEvt << std::endl;
        }
}
void ChargedPFOCorrection::end()
{

	m_pTFile->cd();
	m_pTTree->Write();
	m_Histograms->cd();
	h_nClusters_nTracks->Write();
	h_pfoCharge_nTracks->Write();
	h_InnermostRadiusHit_Neutral->Write();
	h_InnermostRadiusHit_Charged->Write();
	h_FirstSubDet_Charged->Write();
	m_pTFile->Close();
	delete m_pTFile;

	std::cout << "" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////	processed events: 	" << m_nEvtSum << "	////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	std::cout << "" << std::endl;

}
