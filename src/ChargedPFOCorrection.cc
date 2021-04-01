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
	m_Bfield(0.f)
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

	registerOutputCollection(	LCIO::RECONSTRUCTEDPARTICLE,
					"CorrectedPfoCollection",
					"Name of output pfo collection",
					m_outputPfoCollection,
					std::string("CorrectedPfoCollection")
				);

}

void ChargedPFOCorrection::init()
{
	streamlog_out(DEBUG) << "   init called  " << std::endl;
	m_Bfield = MarlinUtil::getBzAtOrigin();
	printParameters();
	streamlog_out(DEBUG) << " BField =  "<< m_Bfield << " Tesla" << std::endl ;
	m_nRun = 0 ;
	m_nEvt = 0 ;
	m_nRunSum = 0;
	m_nEvtSum = 0;

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
	int n_PFO = -1;
	int n_TRK = -1;
	int n_TRKp = -1;
	int n_TRKk = -1;
	int n_MCP = -1;
	this->Clear();

        try
        {
		inputPfoCollection = pLCEvent->getCollection(m_inputPfoCollection);
		MarlinTrkTracks = pLCEvent->getCollection(m_MarlinTrkTracks);
		MarlinTrkTracksKAON = pLCEvent->getCollection(m_MarlinTrkTracksKAON);
		MarlinTrkTracksPROTON = pLCEvent->getCollection(m_MarlinTrkTracksPROTON);
		MCParticleCollection = pLCEvent->getCollection(m_mcParticleCollection);

		n_PFO = inputPfoCollection->getNumberOfElements();
		n_TRK = MarlinTrkTracks->getNumberOfElements();
		n_TRKk = MarlinTrkTracksKAON->getNumberOfElements();
		n_TRKp = MarlinTrkTracksPROTON->getNumberOfElements();
		n_MCP = MCParticleCollection->getNumberOfElements();
		if ( n_PFO == -1 ) streamlog_out(DEBUG) << "Input PFO collection (" << m_inputPfoCollection << ") has no element (PFO) " << std::endl;
		if ( n_TRK == -1 ) streamlog_out(DEBUG) << "Input TRACK collection (" << m_MarlinTrkTracks << ") has no element (Track) " << std::endl;
		if ( n_TRK != n_TRKk ) streamlog_out(DEBUG) << "Input TRACK collection (" << m_MarlinTrkTracksKAON << ") has in-equal number of elements with Main Track Collection (" << m_MarlinTrkTracks << ")!" << std::endl;
		if ( n_TRK != n_TRKp ) streamlog_out(DEBUG) << "Input TRACK collection (" << MarlinTrkTracksPROTON << ") has in-equal number of elements with Main Track Collection (" << m_MarlinTrkTracks << ")!" << std::endl;

		streamlog_out(DEBUG) << "Total Number of PFOs: " << n_PFO << std::endl;
		streamlog_out(DEBUG) << "Total Number of Tracks: " << n_TRK << std::endl;
		streamlog_out(DEBUG) << "Total Number of KaonTracks: " << n_TRKk << std::endl;
		streamlog_out(DEBUG) << "Total Number of ProtonTracks: " << n_TRKp << std::endl;
		streamlog_out(DEBUG) << "Total Number of MCParticles: " << n_MCP << std::endl;

		LCCollectionVec *m_col_outputPfo = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

	}
	catch(DataNotAvailableException &e)
        {
          streamlog_out(MESSAGE) << "Input collection not found in event " << m_nEvt << std::endl;
        }

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
		streamlog_out(DEBUG) << " CHECK : processed events: " << m_nEvtSum << " (Number of inputPFOS: " << n_inputPFOs << " , Number of outputPFOs: " << n_outputPFOs <<")" << std::endl;
	}
	catch(DataNotAvailableException &e)
        {
          streamlog_out(MESSAGE) << "Input/Output collection not found in event " << m_nEvt << std::endl;
        }
}
void ChargedPFOCorrection::end()
{

	std::cout << " END : processed events: " << m_nEvtSum << std::endl;

}
