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
#include "TF1.h"
#include "TTree.h"
#include "TPaveStats.h"
#include "TStyle.h"
#include "TPad.h"

using namespace lcio ;
using namespace marlin ;

ChargedPFOCorrection aChargedPFOCorrection;

ChargedPFOCorrection::ChargedPFOCorrection() :

	Processor("ChargedPFOCorrection"),
	m_nRun(0),
	m_nEvt(0),
	m_nRunSum(0),
	m_nEvtSum(0),
	n_OldPFOS1trk_NormalizedResidualPx(0),
	n_OldPFOS1trk_NormalizedResidualPxPy(0),
	n_OldPFOS1trk_NormalizedResidualPy(0),
	n_OldPFOS1trk_NormalizedResidualPxPz(0),
	n_OldPFOS1trk_NormalizedResidualPyPz(0),
	n_OldPFOS1trk_NormalizedResidualPz(0),
	n_OldPFOS1trk_NormalizedResidualPxE(0),
	n_OldPFOS1trk_NormalizedResidualPyE(0),
	n_OldPFOS1trk_NormalizedResidualPzE(0),
	n_OldPFOS1trk_NormalizedResidualE(0),
	n_OldPFOS1trk_NormalizedResidualTheta(0),
	n_OldPFOS1trk_NormalizedResidualPhi(0),
	n_NewPFOS1trk_NormalizedResidualPx(0),
	n_NewPFOS1trk_NormalizedResidualPxPy(0),
	n_NewPFOS1trk_NormalizedResidualPy(0),
	n_NewPFOS1trk_NormalizedResidualPxPz(0),
	n_NewPFOS1trk_NormalizedResidualPyPz(0),
	n_NewPFOS1trk_NormalizedResidualPz(0),
	n_NewPFOS1trk_NormalizedResidualPxE(0),
	n_NewPFOS1trk_NormalizedResidualPyE(0),
	n_NewPFOS1trk_NormalizedResidualPzE(0),
	n_NewPFOS1trk_NormalizedResidualE(0),
	n_NewPFOS1trk_NormalizedResidualTheta(0),
	n_NewPFOS1trk_NormalizedResidualPhi(0),
	n_NewPFOS2trk_NormalizedResidualPx(0),
	n_NewPFOS2trk_NormalizedResidualPxPy(0),
	n_NewPFOS2trk_NormalizedResidualPy(0),
	n_NewPFOS2trk_NormalizedResidualPxPz(0),
	n_NewPFOS2trk_NormalizedResidualPyPz(0),
	n_NewPFOS2trk_NormalizedResidualPz(0),
	n_NewPFOS2trk_NormalizedResidualPxE(0),
	n_NewPFOS2trk_NormalizedResidualPyE(0),
	n_NewPFOS2trk_NormalizedResidualPzE(0),
	n_NewPFOS2trk_NormalizedResidualE(0),
	n_NewPFOS2trk_NormalizedResidualTheta(0),
	n_NewPFOS2trk_NormalizedResidualPhi(0),
	n_StdTrkProton_NormalizedResidualPx(0),
	n_StdTrkProton_NormalizedResidualPxPy(0),
	n_StdTrkProton_NormalizedResidualPy(0),
	n_StdTrkProton_NormalizedResidualPxPz(0),
	n_StdTrkProton_NormalizedResidualPyPz(0),
	n_StdTrkProton_NormalizedResidualPz(0),
	n_StdTrkProton_NormalizedResidualPxE(0),
	n_StdTrkProton_NormalizedResidualPyE(0),
	n_StdTrkProton_NormalizedResidualPzE(0),
	n_StdTrkProton_NormalizedResidualE(0),
	n_StdTrkProton_NormalizedResidualTheta(0),
	n_StdTrkProton_NormalizedResidualPhi(0),
	n_RFTrkProton_NormalizedResidualPx(0),
	n_RFTrkProton_NormalizedResidualPxPy(0),
	n_RFTrkProton_NormalizedResidualPy(0),
	n_RFTrkProton_NormalizedResidualPxPz(0),
	n_RFTrkProton_NormalizedResidualPyPz(0),
	n_RFTrkProton_NormalizedResidualPz(0),
	n_RFTrkProton_NormalizedResidualPxE(0),
	n_RFTrkProton_NormalizedResidualPyE(0),
	n_RFTrkProton_NormalizedResidualPzE(0),
	n_RFTrkProton_NormalizedResidualE(0),
	n_RFTrkProton_NormalizedResidualTheta(0),
	n_RFTrkProton_NormalizedResidualPhi(0),
	n_StdTrkKaon_NormalizedResidualPx(0),
	n_StdTrkKaon_NormalizedResidualPxPy(0),
	n_StdTrkKaon_NormalizedResidualPy(0),
	n_StdTrkKaon_NormalizedResidualPxPz(0),
	n_StdTrkKaon_NormalizedResidualPyPz(0),
	n_StdTrkKaon_NormalizedResidualPz(0),
	n_StdTrkKaon_NormalizedResidualPxE(0),
	n_StdTrkKaon_NormalizedResidualPyE(0),
	n_StdTrkKaon_NormalizedResidualPzE(0),
	n_StdTrkKaon_NormalizedResidualE(0),
	n_StdTrkKaon_NormalizedResidualTheta(0),
	n_StdTrkKaon_NormalizedResidualPhi(0),
	n_RFTrkKaon_NormalizedResidualPx(0),
	n_RFTrkKaon_NormalizedResidualPxPy(0),
	n_RFTrkKaon_NormalizedResidualPy(0),
	n_RFTrkKaon_NormalizedResidualPxPz(0),
	n_RFTrkKaon_NormalizedResidualPyPz(0),
	n_RFTrkKaon_NormalizedResidualPz(0),
	n_RFTrkKaon_NormalizedResidualPxE(0),
	n_RFTrkKaon_NormalizedResidualPyE(0),
	n_RFTrkKaon_NormalizedResidualPzE(0),
	n_RFTrkKaon_NormalizedResidualE(0),
	n_RFTrkKaon_NormalizedResidualTheta(0),
	n_RFTrkKaon_NormalizedResidualPhi(0),
	n_StdTrk_NormalizedResidualPx(0),
	n_StdTrk_NormalizedResidualPxPy(0),
	n_StdTrk_NormalizedResidualPy(0),
	n_StdTrk_NormalizedResidualPxPz(0),
	n_StdTrk_NormalizedResidualPyPz(0),
	n_StdTrk_NormalizedResidualPz(0),
	n_StdTrk_NormalizedResidualPxE(0),
	n_StdTrk_NormalizedResidualPyE(0),
	n_StdTrk_NormalizedResidualPzE(0),
	n_StdTrk_NormalizedResidualE(0),
	n_StdTrk_NormalizedResidualTheta(0),
	n_StdTrk_NormalizedResidualPhi(0),
	n_RFTrk_NormalizedResidualPx(0),
	n_RFTrk_NormalizedResidualPxPy(0),
	n_RFTrk_NormalizedResidualPy(0),
	n_RFTrk_NormalizedResidualPxPz(0),
	n_RFTrk_NormalizedResidualPyPz(0),
	n_RFTrk_NormalizedResidualPz(0),
	n_RFTrk_NormalizedResidualPxE(0),
	n_RFTrk_NormalizedResidualPyE(0),
	n_RFTrk_NormalizedResidualPzE(0),
	n_RFTrk_NormalizedResidualE(0),
	n_RFTrk_NormalizedResidualTheta(0),
	n_RFTrk_NormalizedResidualPhi(0),
	m_Bfield(0.f),
	c(0.),
	mm2m(0.),
	eV2GeV(0.),
	eB(0.),
	m_MinWeightTrackMCTruthLink(0.),
	m_MinWeightMCTruthTrackLink(0.),
	m_pion_mass(0.),
	m_proton_mass(0.),
	m_kaon_mass(0.)
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
					"Minimum acceptable weight for Track -> MCParticle Link"  ,
					m_MinWeightTrackMCTruthLink ,
					float(0.9f)
				);

	registerProcessorParameter(	"MinWeightMCTruthTrackLink" ,
					"Minimum acceptable weight for MCParticle -> Track Link"  ,
					m_MinWeightMCTruthTrackLink ,
					float(0.9f)
				);

	registerProcessorParameter(	"PionMass" ,
					"Pion mass for calculating (p,E) and CovMat from track parameters"  ,
					m_pion_mass ,
					float(0.13957018)
				);

	registerProcessorParameter(	"ProtonMass" ,
					"Proton mass for calculating (p,E) and CovMat from track parameters"  ,
					m_proton_mass ,
					float(0.938272088)
				);

	registerProcessorParameter(	"KaonMass" ,
					"Kaon mass for calculating (p,E) and CovMat from track parameters"  ,
					m_kaon_mass ,
					float(0.493677)
				);

	registerProcessorParameter(	"fillRootTree",
					"Fill root tree to check processor performance",
					m_fillRootTree,
					bool(true)
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
/*	proton_mass = 0.938272088;
	kaon_mass = 0.493677;
	pion_mass = 0.13957018;
*/
	m_pTFile = new TFile(m_rootFile.c_str(), "recreate");

	m_pTTree = new TTree("PFOswithRFT", "PFOswithRFT");
	m_pTTree->SetDirectory(m_pTFile);
	m_pTTree->Branch("run", &m_nRun, "run/I");
	m_pTTree->Branch("event", &m_nEvt, "event/I");
	m_pTTree->Branch("pfoCharge", &m_pfoCharge);
	m_pTTree->Branch("nTracksOfPFO", &m_nTracksOfPFO);
	m_pTTree->Branch("TrkToMCPLinkWeight", &m_TrkToMCPLinkWeight);
	m_pTTree->Branch("oldPFO_Px", &m_oldPFO_Px);
	m_pTTree->Branch("oldPFO_Py", &m_oldPFO_Py);
	m_pTTree->Branch("oldPFO_Pz", &m_oldPFO_Pz);
	m_pTTree->Branch("oldPFO_E", &m_oldPFO_E);
	m_pTTree->Branch("oldPFO_Theta", &m_oldPFO_Theta);
	m_pTTree->Branch("oldPFO_Phi", &m_oldPFO_Phi);
	m_pTTree->Branch("newPFO_Px", &m_newPFO_Px);
	m_pTTree->Branch("newPFO_Py", &m_newPFO_Py);
	m_pTTree->Branch("newPFO_Pz", &m_newPFO_Pz);
	m_pTTree->Branch("newPFO_E", &m_newPFO_E);
	m_pTTree->Branch("newPFO_Theta", &m_newPFO_Theta);
	m_pTTree->Branch("newPFO_Phi", &m_newPFO_Phi);
	m_pTTree->Branch("oldPFO_ResidualPx", &m_oldPFO_ResidualPx);
	m_pTTree->Branch("oldPFO_ResidualPy", &m_oldPFO_ResidualPy);
	m_pTTree->Branch("oldPFO_ResidualPz", &m_oldPFO_ResidualPz);
	m_pTTree->Branch("oldPFO_ResidualE", &m_oldPFO_ResidualE);
	m_pTTree->Branch("oldPFO_ResidualTheta", &m_oldPFO_ResidualTheta);
	m_pTTree->Branch("oldPFO_ResidualPhi", &m_oldPFO_ResidualPhi);
	m_pTTree->Branch("newPFO_ResidualPx", &m_newPFO_ResidualPx);
	m_pTTree->Branch("newPFO_ResidualPy", &m_newPFO_ResidualPy);
	m_pTTree->Branch("newPFO_ResidualPz", &m_newPFO_ResidualPz);
	m_pTTree->Branch("newPFO_ResidualE", &m_newPFO_ResidualE);
	m_pTTree->Branch("newPFO_ResidualTheta", &m_newPFO_ResidualTheta);
	m_pTTree->Branch("newPFO_ResidualPhi", &m_newPFO_ResidualPhi);
	m_pTTree->Branch("oldPFO_NormalizedResidualPx", &m_oldPFO_NormalizedResidualPx);
	m_pTTree->Branch("oldPFO_NormalizedResidualPy", &m_oldPFO_NormalizedResidualPy);
	m_pTTree->Branch("oldPFO_NormalizedResidualPz", &m_oldPFO_NormalizedResidualPz);
	m_pTTree->Branch("oldPFO_NormalizedResidualE", &m_oldPFO_NormalizedResidualE);
	m_pTTree->Branch("oldPFO_NormalizedResidualTheta", &m_oldPFO_NormalizedResidualTheta);
	m_pTTree->Branch("oldPFO_NormalizedResidualPhi", &m_oldPFO_NormalizedResidualPhi);
	m_pTTree->Branch("newPFO_NormalizedResidualPx", &m_newPFO_NormalizedResidualPx);
	m_pTTree->Branch("newPFO_NormalizedResidualPy", &m_newPFO_NormalizedResidualPy);
	m_pTTree->Branch("newPFO_NormalizedResidualPz", &m_newPFO_NormalizedResidualPz);
	m_pTTree->Branch("newPFO_NormalizedResidualE", &m_newPFO_NormalizedResidualE);
	m_pTTree->Branch("newPFO_NormalizedResidualTheta", &m_newPFO_NormalizedResidualTheta);
	m_pTTree->Branch("newPFO_NormalizedResidualPhi", &m_newPFO_NormalizedResidualPhi);
	m_Histograms = m_pTFile->mkdir("Histograms");
	m_OldPFOs_1Trk = m_pTFile->mkdir("OldPFOs_1Trk");
	m_NewPFOs_1Trk = m_pTFile->mkdir("NewPFOs_1Trk");
	m_NewPFOs_2Trk = m_pTFile->mkdir("NewPFOs_2Trk");
	m_NewPFOs_nTrk = m_pTFile->mkdir("NewPFOs_nTrk");
	m_TrueProtons_1Trk = m_pTFile->mkdir("TrueProtons_1Trk");
	m_TrueKaons_1Trk = m_pTFile->mkdir("TrueKaons_1Trk");
	m_IndividualTracks_2Trk = m_pTFile->mkdir("IndividualTracks_2Trk");
	h_nClusters_nTracks = new TH2I("All PFOs", "; n_{Clusters}; n_{Tracks}", 11, -0.5, 10.5, 11, -0.5, 10.5);
	h_pfoCharge_nTracks = new TH2I("Neutral PFOs", "; PFO Charge; n_{Tracks}", 5, -2.5, 2.5, 11, -0.5, 10.5);
	h_InnermostRadiusHit_Neutral = new TH2F("Neutral PFOs with 2 tracks", "; r_{innermost hit}^{Track 1} [mm]; r_{innermost hit}^{Track 2} [mm]", 180, 0.0, 1800.0, 180, 0.0, 1800.0);
	h_InnermostRadiusHit_Charged = new TH2F("Charged PFOs with 2 tracks", "; r_{innermost hit}^{Track 1} [mm]; r_{innermost hit}^{Track 2} [mm]", 180, 0.0, 1800.0, 180, 0.0, 1800.0);
	h_FirstSubDet_Charged = new TH2I("First SubDet of Charged PFOs with 2 tracks", "; First SubDet^{Track 1}; First SubDet^{Track 2}", 11, -0.5, 10.5, 11, -0.5, 10.5);
	h_OldPFOS1trk_NormalizedResidualPx = new TH1F( "Std. Track (1 trk)" , "; (_{}p_{x}^{REC} - p_{x}^{MC}) / #sigma_{p_{x}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_OldPFOS1trk_NormalizedResidualPx = 0;
	h_OldPFOS1trk_NormalizedResidualPxPy = new TH1F( "Std. Track (1 trk)" , "; #sqrt{(_{}p_{x}^{REC} - p_{x}^{MC})#times(_{}p_{y}^{REC} - p_{y}^{MC}) / #sigma_{p_{x}p_{y}}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_OldPFOS1trk_NormalizedResidualPxPy = 0;
	h_OldPFOS1trk_NormalizedResidualPy = new TH1F( "Std. Track (1 trk)" , "; (_{}p_{y}^{REC} - p_{y}^{MC}) / #sigma_{p_{y}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_OldPFOS1trk_NormalizedResidualPy = 0;
	h_OldPFOS1trk_NormalizedResidualPxPz = new TH1F( "Std. Track (1 trk)" , "; #sqrt{(_{}p_{x}^{REC} - p_{x}^{MC})#times(_{}p_{z}^{REC} - p_{z}^{MC}) / #sigma_{p_{x}p_{z}}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_OldPFOS1trk_NormalizedResidualPxPz = 0;
	h_OldPFOS1trk_NormalizedResidualPyPz = new TH1F( "Std. Track (1 trk)" , "; #sqrt{(_{}p_{y}^{REC} - p_{y}^{MC})#times(_{}p_{z}^{REC} - p_{z}^{MC}) / #sigma_{p_{y}p_{z}}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_OldPFOS1trk_NormalizedResidualPyPz = 0;
	h_OldPFOS1trk_NormalizedResidualPz = new TH1F( "Std. Track (1 trk)" , "; (_{}p_{z}^{REC} - p_{z}^{MC}) / #sigma_{p_{z}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_OldPFOS1trk_NormalizedResidualPz = 0;
	h_OldPFOS1trk_NormalizedResidualPxE = new TH1F( "Std. Track (1 trk)" , "; #sqrt{(_{}p_{x}^{REC} - p_{x}^{MC})#times(_{}E^{REC} - E^{MC}) / #sigma_{p_{x}E}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_OldPFOS1trk_NormalizedResidualPxE = 0;
	h_OldPFOS1trk_NormalizedResidualPyE = new TH1F( "Std. Track (1 trk)" , "; #sqrt{(_{}p_{y}^{REC} - p_{y}^{MC})#times(_{}E^{REC} - E^{MC}) / #sigma_{p_{y}E}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_OldPFOS1trk_NormalizedResidualPyE = 0;
	h_OldPFOS1trk_NormalizedResidualPzE = new TH1F( "Std. Track (1 trk)" , "; #sqrt{(_{}p_{z}^{REC} - p_{z}^{MC})#times(_{}E^{REC} - E^{MC}) / #sigma_{p_{z}E}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_OldPFOS1trk_NormalizedResidualPzE = 0;
	h_OldPFOS1trk_NormalizedResidualE = new TH1F( "Std. Track (1 trk)" , "; (_{}E^{REC} - E^{MC}) / #sigma_{E}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_OldPFOS1trk_NormalizedResidualE = 0;
	h_OldPFOS1trk_NormalizedResidualTheta = new TH1F( "Std. Track (1 trk)" , "; (_{}#theta^{REC} - #theta^{MC}) / #sigma_{#theta}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_OldPFOS1trk_NormalizedResidualTheta = 0;
	h_OldPFOS1trk_NormalizedResidualPhi = new TH1F( "Std. Track (1 trk)" , "; (_{}#phi^{REC} - #phi^{MC}) / #sigma_{#phi}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_OldPFOS1trk_NormalizedResidualPhi = 0;
	h_NewPFOS1trk_NormalizedResidualPx = new TH1F( "refitted Track (1 trk)" , "; (_{}p_{x}^{REC} - p_{x}^{MC}) / #sigma_{p_{x}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NewPFOS1trk_NormalizedResidualPx = 0;
	h_NewPFOS1trk_NormalizedResidualPxPy = new TH1F( "refitted Track (1 trk)" , "; #sqrt{(_{}p_{x}^{REC} - p_{x}^{MC})#times(_{}p_{y}^{REC} - p_{y}^{MC}) / #sigma_{p_{x}p_{y}}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NewPFOS1trk_NormalizedResidualPxPy = 0;
	h_NewPFOS1trk_NormalizedResidualPy = new TH1F( "refitted Track (1 trk)" , "; (_{}p_{y}^{REC} - p_{y}^{MC}) / #sigma_{p_{y}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NewPFOS1trk_NormalizedResidualPy = 0;
	h_NewPFOS1trk_NormalizedResidualPxPz = new TH1F( "refitted Track (1 trk)" , "; #sqrt{(_{}p_{x}^{REC} - p_{x}^{MC})#times(_{}p_{z}^{REC} - p_{z}^{MC}) / #sigma_{p_{x}p_{z}}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NewPFOS1trk_NormalizedResidualPxPz = 0;
	h_NewPFOS1trk_NormalizedResidualPyPz = new TH1F( "refitted Track (1 trk)" , "; #sqrt{(_{}p_{y}^{REC} - p_{y}^{MC})#times(_{}p_{z}^{REC} - p_{z}^{MC}) / #sigma_{p_{y}p_{z}}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NewPFOS1trk_NormalizedResidualPyPz = 0;
	h_NewPFOS1trk_NormalizedResidualPz = new TH1F( "refitted Track (1 trk)" , "; (_{}p_{z}^{REC} - p_{z}^{MC}) / #sigma_{p_{z}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NewPFOS1trk_NormalizedResidualPz = 0;
	h_NewPFOS1trk_NormalizedResidualPxE = new TH1F( "refitted Track (1 trk)" , "; #sqrt{(_{}p_{x}^{REC} - p_{x}^{MC})#times(_{}E^{REC} - E^{MC}) / #sigma_{p_{x}E}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NewPFOS1trk_NormalizedResidualPxE = 0;
	h_NewPFOS1trk_NormalizedResidualPyE = new TH1F( "refitted Track (1 trk)" , "; #sqrt{(_{}p_{y}^{REC} - p_{y}^{MC})#times(_{}E^{REC} - E^{MC}) / #sigma_{p_{y}E}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NewPFOS1trk_NormalizedResidualPyE = 0;
	h_NewPFOS1trk_NormalizedResidualPzE = new TH1F( "refitted Track (1 trk)" , "; #sqrt{(_{}p_{z}^{REC} - p_{z}^{MC})#times(_{}E^{REC} - E^{MC}) / #sigma_{p_{z}E}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NewPFOS1trk_NormalizedResidualPzE = 0;
	h_NewPFOS1trk_NormalizedResidualE = new TH1F( "refitted Track (1 trk)" , "; (_{}E^{REC} - E^{MC}) / #sigma_{E}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NewPFOS1trk_NormalizedResidualE = 0;
	h_NewPFOS1trk_NormalizedResidualTheta = new TH1F( "refitted Track (1 trk)" , "; (_{}#theta^{REC} - #theta^{MC}) / #sigma_{#theta}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NewPFOS1trk_NormalizedResidualTheta = 0;
	h_NewPFOS1trk_NormalizedResidualPhi = new TH1F( "refitted Track (1 trk)" , "; (_{}#phi^{REC} - #phi^{MC}) / #sigma_{#phi}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NewPFOS1trk_NormalizedResidualPhi = 0;
	h_NewPFOS2trk_NormalizedResidualPx = new TH1F( "refitted Track (2 trk)" , "; (_{}p_{x}^{REC} - p_{x}^{MC}) / #sigma_{p_{x}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NewPFOS2trk_NormalizedResidualPx = 0;
	h_NewPFOS2trk_NormalizedResidualPxPy = new TH1F( "refitted Track (2 trk)" , "; #sqrt{(_{}p_{x}^{REC} - p_{x}^{MC})#times(_{}p_{y}^{REC} - p_{y}^{MC}) / #sigma_{p_{x}p_{y}}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NewPFOS2trk_NormalizedResidualPxPy = 0;
	h_NewPFOS2trk_NormalizedResidualPy = new TH1F( "refitted Track (2 trk)" , "; (_{}p_{y}^{REC} - p_{y}^{MC}) / #sigma_{p_{y}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NewPFOS2trk_NormalizedResidualPy = 0;
	h_NewPFOS2trk_NormalizedResidualPxPz = new TH1F( "refitted Track (2 trk)" , "; #sqrt{(_{}p_{x}^{REC} - p_{x}^{MC})#times(_{}p_{z}^{REC} - p_{z}^{MC}) / #sigma_{p_{x}p_{z}}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NewPFOS2trk_NormalizedResidualPxPz = 0;
	h_NewPFOS2trk_NormalizedResidualPyPz = new TH1F( "refitted Track (2 trk)" , "; #sqrt{(_{}p_{y}^{REC} - p_{y}^{MC})#times(_{}p_{z}^{REC} - p_{z}^{MC}) / #sigma_{p_{y}p_{z}}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NewPFOS2trk_NormalizedResidualPyPz = 0;
	h_NewPFOS2trk_NormalizedResidualPz = new TH1F( "refitted Track (2 trk)" , "; (_{}p_{z}^{REC} - p_{z}^{MC}) / #sigma_{p_{z}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NewPFOS2trk_NormalizedResidualPz = 0;
	h_NewPFOS2trk_NormalizedResidualPxE = new TH1F( "refitted Track (2 trk)" , "; #sqrt{(_{}p_{x}^{REC} - p_{x}^{MC})#times(_{}E^{REC} - E^{MC}) / #sigma_{p_{x}E}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NewPFOS2trk_NormalizedResidualPxE = 0;
	h_NewPFOS2trk_NormalizedResidualPyE = new TH1F( "refitted Track (2 trk)" , "; #sqrt{(_{}p_{y}^{REC} - p_{y}^{MC})#times(_{}E^{REC} - E^{MC}) / #sigma_{p_{y}E}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NewPFOS2trk_NormalizedResidualPyE = 0;
	h_NewPFOS2trk_NormalizedResidualPzE = new TH1F( "refitted Track (2 trk)" , "; #sqrt{(_{}p_{z}^{REC} - p_{z}^{MC})#times(_{}E^{REC} - E^{MC}) / #sigma_{p_{z}E}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NewPFOS2trk_NormalizedResidualPzE = 0;
	h_NewPFOS2trk_NormalizedResidualE = new TH1F( "refitted Track (2 trk)" , "; (_{}E^{REC} - E^{MC}) / #sigma_{E}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NewPFOS2trk_NormalizedResidualE = 0;
	h_NewPFOS2trk_NormalizedResidualTheta = new TH1F( "refitted Track (2 trk)" , "; (_{}#theta^{REC} - #theta^{MC}) / #sigma_{#theta}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NewPFOS2trk_NormalizedResidualTheta = 0;
	h_NewPFOS2trk_NormalizedResidualPhi = new TH1F( "refitted Track (2 trk)" , "; (_{}#phi^{REC} - #phi^{MC}) / #sigma_{#phi}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NewPFOS2trk_NormalizedResidualPhi = 0;
	h_StdTrkProton_NormalizedResidualPx = new TH1F( "Protons (Std. Track)" , "; (_{}p_{x}^{REC} - p_{x}^{MC}) / #sigma_{p_{x}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_StdTrkProton_NormalizedResidualPx = 0;
	h_StdTrkProton_NormalizedResidualPxPy = new TH1F( "Protons (Std. Track)" , "; #sqrt{(_{}p_{x}^{REC} - p_{x}^{MC})#times(_{}p_{y}^{REC} - p_{y}^{MC}) / #sigma_{p_{x}p_{y}}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_StdTrkProton_NormalizedResidualPxPy = 0;
	h_StdTrkProton_NormalizedResidualPy = new TH1F( "Protons (Std. Track)" , "; (_{}p_{y}^{REC} - p_{y}^{MC}) / #sigma_{p_{y}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_StdTrkProton_NormalizedResidualPy = 0;
	h_StdTrkProton_NormalizedResidualPxPz = new TH1F( "Protons (Std. Track)" , "; #sqrt{(_{}p_{x}^{REC} - p_{x}^{MC})#times(_{}p_{z}^{REC} - p_{z}^{MC}) / #sigma_{p_{x}p_{z}}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_StdTrkProton_NormalizedResidualPxPz = 0;
	h_StdTrkProton_NormalizedResidualPyPz = new TH1F( "Protons (Std. Track)" , "; #sqrt{(_{}p_{y}^{REC} - p_{y}^{MC})#times(_{}p_{z}^{REC} - p_{z}^{MC}) / #sigma_{p_{y}p_{z}}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_StdTrkProton_NormalizedResidualPyPz = 0;
	h_StdTrkProton_NormalizedResidualPz = new TH1F( "Protons (Std. Track)" , "; (_{}p_{z}^{REC} - p_{z}^{MC}) / #sigma_{p_{z}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_StdTrkProton_NormalizedResidualPz = 0;
	h_StdTrkProton_NormalizedResidualPxE = new TH1F( "Protons (Std. Track)" , "; #sqrt{(_{}p_{x}^{REC} - p_{x}^{MC})#times(_{}E^{REC} - E^{MC}) / #sigma_{p_{x}E}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_StdTrkProton_NormalizedResidualPxE = 0;
	h_StdTrkProton_NormalizedResidualPyE = new TH1F( "Protons (Std. Track)" , "; #sqrt{(_{}p_{y}^{REC} - p_{y}^{MC})#times(_{}E^{REC} - E^{MC}) / #sigma_{p_{y}E}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_StdTrkProton_NormalizedResidualPyE = 0;
	h_StdTrkProton_NormalizedResidualPzE = new TH1F( "Protons (Std. Track)" , "; #sqrt{(_{}p_{z}^{REC} - p_{z}^{MC})#times(_{}E^{REC} - E^{MC}) / #sigma_{p_{z}E}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_StdTrkProton_NormalizedResidualPzE = 0;
	h_StdTrkProton_NormalizedResidualE = new TH1F( "Protons (Std. Track)" , "; (_{}E^{REC} - E^{MC}) / #sigma_{E}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_StdTrkProton_NormalizedResidualE = 0;
	h_StdTrkProton_NormalizedResidualTheta = new TH1F( "Protons (Std. Track)" , "; (_{}#theta^{REC} - #theta^{MC}) / #sigma_{#theta}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_StdTrkProton_NormalizedResidualTheta = 0;
	h_StdTrkProton_NormalizedResidualPhi = new TH1F( "Protons (Std. Track)" , "; (_{}#phi^{REC} - #phi^{MC}) / #sigma_{#phi}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_StdTrkProton_NormalizedResidualPhi = 0;
	h_RFTrkProton_NormalizedResidualPx = new TH1F( "Protons (refitted Track)" , "; (_{}p_{x}^{REC} - p_{x}^{MC}) / #sigma_{p_{x}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_RFTrkProton_NormalizedResidualPx = 0;
	h_RFTrkProton_NormalizedResidualPxPy = new TH1F( "Protons (refitted Track)" , "; #sqrt{(_{}p_{x}^{REC} - p_{x}^{MC})#times(_{}p_{y}^{REC} - p_{y}^{MC}) / #sigma_{p_{x}p_{y}}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_RFTrkProton_NormalizedResidualPxPy = 0;
	h_RFTrkProton_NormalizedResidualPy = new TH1F( "Protons (refitted Track)" , "; (_{}p_{y}^{REC} - p_{y}^{MC}) / #sigma_{p_{y}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_RFTrkProton_NormalizedResidualPy = 0;
	h_RFTrkProton_NormalizedResidualPxPz = new TH1F( "Protons (refitted Track)" , "; #sqrt{(_{}p_{x}^{REC} - p_{x}^{MC})#times(_{}p_{z}^{REC} - p_{z}^{MC}) / #sigma_{p_{x}p_{z}}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_RFTrkProton_NormalizedResidualPxPz = 0;
	h_RFTrkProton_NormalizedResidualPyPz = new TH1F( "Protons (refitted Track)" , "; #sqrt{(_{}p_{y}^{REC} - p_{y}^{MC})#times(_{}p_{z}^{REC} - p_{z}^{MC}) / #sigma_{p_{y}p_{z}}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_RFTrkProton_NormalizedResidualPyPz = 0;
	h_RFTrkProton_NormalizedResidualPz = new TH1F( "Protons (refitted Track)" , "; (_{}p_{z}^{REC} - p_{z}^{MC}) / #sigma_{p_{z}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_RFTrkProton_NormalizedResidualPz = 0;
	h_RFTrkProton_NormalizedResidualPxE = new TH1F( "Protons (refitted Track)" , "; #sqrt{(_{}p_{x}^{REC} - p_{x}^{MC})#times(_{}E^{REC} - E^{MC}) / #sigma_{p_{x}E}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_RFTrkProton_NormalizedResidualPxE = 0;
	h_RFTrkProton_NormalizedResidualPyE = new TH1F( "Protons (refitted Track)" , "; #sqrt{(_{}p_{y}^{REC} - p_{y}^{MC})#times(_{}E^{REC} - E^{MC}) / #sigma_{p_{y}E}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_RFTrkProton_NormalizedResidualPyE = 0;
	h_RFTrkProton_NormalizedResidualPzE = new TH1F( "Protons (refitted Track)" , "; #sqrt{(_{}p_{z}^{REC} - p_{z}^{MC})#times(_{}E^{REC} - E^{MC}) / #sigma_{p_{z}E}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_RFTrkProton_NormalizedResidualPzE = 0;
	h_RFTrkProton_NormalizedResidualE = new TH1F( "Protons (refitted Track)" , "; (_{}E^{REC} - E^{MC}) / #sigma_{E}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_RFTrkProton_NormalizedResidualE = 0;
	h_RFTrkProton_NormalizedResidualTheta = new TH1F( "Protons (refitted Track)" , "; (_{}#theta^{REC} - #theta^{MC}) / #sigma_{#theta}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_RFTrkProton_NormalizedResidualTheta = 0;
	h_RFTrkProton_NormalizedResidualPhi = new TH1F( "Protons (refitted Track)" , "; (_{}#phi^{REC} - #phi^{MC}) / #sigma_{#phi}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_RFTrkProton_NormalizedResidualPhi = 0;
	h_StdTrkKaon_NormalizedResidualPx = new TH1F( "Kaons (Std. Track)" , "; (_{}p_{x}^{REC} - p_{x}^{MC}) / #sigma_{p_{x}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_StdTrkKaon_NormalizedResidualPx = 0;
	h_StdTrkKaon_NormalizedResidualPxPy = new TH1F( "Kaons (Std. Track)" , "; #sqrt{(_{}p_{x}^{REC} - p_{x}^{MC})#times(_{}p_{y}^{REC} - p_{y}^{MC}) / #sigma_{p_{x}p_{y}}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_StdTrkKaon_NormalizedResidualPxPy = 0;
	h_StdTrkKaon_NormalizedResidualPy = new TH1F( "Kaons (Std. Track)" , "; (_{}p_{y}^{REC} - p_{y}^{MC}) / #sigma_{p_{y}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_StdTrkKaon_NormalizedResidualPy = 0;
	h_StdTrkKaon_NormalizedResidualPxPz = new TH1F( "Kaons (Std. Track)" , "; #sqrt{(_{}p_{x}^{REC} - p_{x}^{MC})#times(_{}p_{z}^{REC} - p_{z}^{MC}) / #sigma_{p_{x}p_{z}}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_StdTrkKaon_NormalizedResidualPxPz = 0;
	h_StdTrkKaon_NormalizedResidualPyPz = new TH1F( "Kaons (Std. Track)" , "; #sqrt{(_{}p_{y}^{REC} - p_{y}^{MC})#times(_{}p_{z}^{REC} - p_{z}^{MC}) / #sigma_{p_{y}p_{z}}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_StdTrkKaon_NormalizedResidualPyPz = 0;
	h_StdTrkKaon_NormalizedResidualPz = new TH1F( "Kaons (Std. Track)" , "; (_{}p_{z}^{REC} - p_{z}^{MC}) / #sigma_{p_{z}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_StdTrkKaon_NormalizedResidualPz = 0;
	h_StdTrkKaon_NormalizedResidualPxE = new TH1F( "Kaons (Std. Track)" , "; #sqrt{(_{}p_{x}^{REC} - p_{x}^{MC})#times(_{}E^{REC} - E^{MC}) / #sigma_{p_{x}E}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_StdTrkKaon_NormalizedResidualPxE = 0;
	h_StdTrkKaon_NormalizedResidualPyE = new TH1F( "Kaons (Std. Track)" , "; #sqrt{(_{}p_{y}^{REC} - p_{y}^{MC})#times(_{}E^{REC} - E^{MC}) / #sigma_{p_{y}E}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_StdTrkKaon_NormalizedResidualPyE = 0;
	h_StdTrkKaon_NormalizedResidualPzE = new TH1F( "Kaons (Std. Track)" , "; #sqrt{(_{}p_{z}^{REC} - p_{z}^{MC})#times(_{}E^{REC} - E^{MC}) / #sigma_{p_{z}E}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_StdTrkKaon_NormalizedResidualPzE = 0;
	h_StdTrkKaon_NormalizedResidualE = new TH1F( "Kaons (Std. Track)" , "; (_{}E^{REC} - E^{MC}) / #sigma_{E}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_StdTrkKaon_NormalizedResidualE = 0;
	h_StdTrkKaon_NormalizedResidualTheta = new TH1F( "Kaons (Std. Track)" , "; (_{}#theta^{REC} - #theta^{MC}) / #sigma_{#theta}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_StdTrkKaon_NormalizedResidualTheta = 0;
	h_StdTrkKaon_NormalizedResidualPhi = new TH1F( "Kaons (Std. Track)" , "; (_{}#phi^{REC} - #phi^{MC}) / #sigma_{#phi}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_StdTrkKaon_NormalizedResidualPhi = 0;
	h_RFTrkKaon_NormalizedResidualPx = new TH1F( "Kaons (refitted Track)" , "; (_{}p_{x}^{REC} - p_{x}^{MC}) / #sigma_{p_{x}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_RFTrkKaon_NormalizedResidualPx = 0;
	h_RFTrkKaon_NormalizedResidualPxPy = new TH1F( "Kaons (refitted Track)" , "; #sqrt{(_{}p_{x}^{REC} - p_{x}^{MC})#times(_{}p_{y}^{REC} - p_{y}^{MC}) / #sigma_{p_{x}p_{y}}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_RFTrkKaon_NormalizedResidualPxPy = 0;
	h_RFTrkKaon_NormalizedResidualPy = new TH1F( "Kaons (refitted Track)" , "; (_{}p_{y}^{REC} - p_{y}^{MC}) / #sigma_{p_{y}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_RFTrkKaon_NormalizedResidualPy = 0;
	h_RFTrkKaon_NormalizedResidualPxPz = new TH1F( "Kaons (refitted Track)" , "; #sqrt{(_{}p_{x}^{REC} - p_{x}^{MC})#times(_{}p_{z}^{REC} - p_{z}^{MC}) / #sigma_{p_{x}p_{z}}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_RFTrkKaon_NormalizedResidualPxPz = 0;
	h_RFTrkKaon_NormalizedResidualPyPz = new TH1F( "Kaons (refitted Track)" , "; #sqrt{(_{}p_{y}^{REC} - p_{y}^{MC})#times(_{}p_{z}^{REC} - p_{z}^{MC}) / #sigma_{p_{y}p_{z}}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_RFTrkKaon_NormalizedResidualPyPz = 0;
	h_RFTrkKaon_NormalizedResidualPz = new TH1F( "Kaons (refitted Track)" , "; (_{}p_{z}^{REC} - p_{z}^{MC}) / #sigma_{p_{z}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_RFTrkKaon_NormalizedResidualPz = 0;
	h_RFTrkKaon_NormalizedResidualPxE = new TH1F( "Kaons (refitted Track)" , "; #sqrt{(_{}p_{x}^{REC} - p_{x}^{MC})#times(_{}E^{REC} - E^{MC}) / #sigma_{p_{x}E}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_RFTrkKaon_NormalizedResidualPxE = 0;
	h_RFTrkKaon_NormalizedResidualPyE = new TH1F( "Kaons (refitted Track)" , "; #sqrt{(_{}p_{y}^{REC} - p_{y}^{MC})#times(_{}E^{REC} - E^{MC}) / #sigma_{p_{y}E}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_RFTrkKaon_NormalizedResidualPyE = 0;
	h_RFTrkKaon_NormalizedResidualPzE = new TH1F( "Kaons (refitted Track)" , "; #sqrt{(_{}p_{z}^{REC} - p_{z}^{MC})#times(_{}E^{REC} - E^{MC}) / #sigma_{p_{z}E}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_RFTrkKaon_NormalizedResidualPzE = 0;
	h_RFTrkKaon_NormalizedResidualE = new TH1F( "Kaons (refitted Track)" , "; (_{}E^{REC} - E^{MC}) / #sigma_{E}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_RFTrkKaon_NormalizedResidualE = 0;
	h_RFTrkKaon_NormalizedResidualTheta = new TH1F( "Kaons (refitted Track)" , "; (_{}#theta^{REC} - #theta^{MC}) / #sigma_{#theta}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_RFTrkKaon_NormalizedResidualTheta = 0;
	h_RFTrkKaon_NormalizedResidualPhi = new TH1F( "Kaons (refitted Track)" , "; (_{}#phi^{REC} - #phi^{MC}) / #sigma_{#phi}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_RFTrkKaon_NormalizedResidualPhi = 0;
	h_StdTrk_NormalizedResidualPx = new TH1F( "Kaons (Std. Track)" , "; (_{}p_{x}^{REC} - p_{x}^{MC}) / #sigma_{p_{x}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_StdTrk_NormalizedResidualPx = 0;
	h_StdTrk_NormalizedResidualPxPy = new TH1F( "Kaons (Std. Track)" , "; #sqrt{(_{}p_{x}^{REC} - p_{x}^{MC})#times(_{}p_{y}^{REC} - p_{y}^{MC}) / #sigma_{p_{x}p_{y}}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_StdTrk_NormalizedResidualPxPy = 0;
	h_StdTrk_NormalizedResidualPy = new TH1F( "Kaons (Std. Track)" , "; (_{}p_{y}^{REC} - p_{y}^{MC}) / #sigma_{p_{y}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_StdTrk_NormalizedResidualPy = 0;
	h_StdTrk_NormalizedResidualPxPz = new TH1F( "Kaons (Std. Track)" , "; #sqrt{(_{}p_{x}^{REC} - p_{x}^{MC})#times(_{}p_{z}^{REC} - p_{z}^{MC}) / #sigma_{p_{x}p_{z}}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_StdTrk_NormalizedResidualPxPz = 0;
	h_StdTrk_NormalizedResidualPyPz = new TH1F( "Kaons (Std. Track)" , "; #sqrt{(_{}p_{y}^{REC} - p_{y}^{MC})#times(_{}p_{z}^{REC} - p_{z}^{MC}) / #sigma_{p_{y}p_{z}}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_StdTrk_NormalizedResidualPyPz = 0;
	h_StdTrk_NormalizedResidualPz = new TH1F( "Kaons (Std. Track)" , "; (_{}p_{z}^{REC} - p_{z}^{MC}) / #sigma_{p_{z}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_StdTrk_NormalizedResidualPz = 0;
	h_StdTrk_NormalizedResidualPxE = new TH1F( "Kaons (Std. Track)" , "; #sqrt{(_{}p_{x}^{REC} - p_{x}^{MC})#times(_{}E^{REC} - E^{MC}) / #sigma_{p_{x}E}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_StdTrk_NormalizedResidualPxE = 0;
	h_StdTrk_NormalizedResidualPyE = new TH1F( "Kaons (Std. Track)" , "; #sqrt{(_{}p_{y}^{REC} - p_{y}^{MC})#times(_{}E^{REC} - E^{MC}) / #sigma_{p_{y}E}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_StdTrk_NormalizedResidualPyE = 0;
	h_StdTrk_NormalizedResidualPzE = new TH1F( "Kaons (Std. Track)" , "; #sqrt{(_{}p_{z}^{REC} - p_{z}^{MC})#times(_{}E^{REC} - E^{MC}) / #sigma_{p_{z}E}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_StdTrk_NormalizedResidualPzE = 0;
	h_StdTrk_NormalizedResidualE = new TH1F( "Kaons (Std. Track)" , "; (_{}E^{REC} - E^{MC}) / #sigma_{E}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_StdTrk_NormalizedResidualE = 0;
	h_StdTrk_NormalizedResidualTheta = new TH1F( "Kaons (Std. Track)" , "; (_{}#theta^{REC} - #theta^{MC}) / #sigma_{#theta}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_StdTrk_NormalizedResidualTheta = 0;
	h_StdTrk_NormalizedResidualPhi = new TH1F( "Kaons (Std. Track)" , "; (_{}#phi^{REC} - #phi^{MC}) / #sigma_{#phi}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_StdTrk_NormalizedResidualPhi = 0;
	h_RFTrk_NormalizedResidualPx = new TH1F( "Kaons (refitted Track)" , "; (_{}p_{x}^{REC} - p_{x}^{MC}) / #sigma_{p_{x}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_RFTrk_NormalizedResidualPx = 0;
	h_RFTrk_NormalizedResidualPxPy = new TH1F( "Kaons (refitted Track)" , "; #sqrt{(_{}p_{x}^{REC} - p_{x}^{MC})#times(_{}p_{y}^{REC} - p_{y}^{MC}) / #sigma_{p_{x}p_{y}}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_RFTrk_NormalizedResidualPxPy = 0;
	h_RFTrk_NormalizedResidualPy = new TH1F( "Kaons (refitted Track)" , "; (_{}p_{y}^{REC} - p_{y}^{MC}) / #sigma_{p_{y}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_RFTrk_NormalizedResidualPy = 0;
	h_RFTrk_NormalizedResidualPxPz = new TH1F( "Kaons (refitted Track)" , "; #sqrt{(_{}p_{x}^{REC} - p_{x}^{MC})#times(_{}p_{z}^{REC} - p_{z}^{MC}) / #sigma_{p_{x}p_{z}}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_RFTrk_NormalizedResidualPxPz = 0;
	h_RFTrk_NormalizedResidualPyPz = new TH1F( "Kaons (refitted Track)" , "; #sqrt{(_{}p_{y}^{REC} - p_{y}^{MC})#times(_{}p_{z}^{REC} - p_{z}^{MC}) / #sigma_{p_{y}p_{z}}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_RFTrk_NormalizedResidualPyPz = 0;
	h_RFTrk_NormalizedResidualPz = new TH1F( "Kaons (refitted Track)" , "; (_{}p_{z}^{REC} - p_{z}^{MC}) / #sigma_{p_{z}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_RFTrk_NormalizedResidualPz = 0;
	h_RFTrk_NormalizedResidualPxE = new TH1F( "Kaons (refitted Track)" , "; #sqrt{(_{}p_{x}^{REC} - p_{x}^{MC})#times(_{}E^{REC} - E^{MC}) / #sigma_{p_{x}E}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_RFTrk_NormalizedResidualPxE = 0;
	h_RFTrk_NormalizedResidualPyE = new TH1F( "Kaons (refitted Track)" , "; #sqrt{(_{}p_{y}^{REC} - p_{y}^{MC})#times(_{}E^{REC} - E^{MC}) / #sigma_{p_{y}E}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_RFTrk_NormalizedResidualPyE = 0;
	h_RFTrk_NormalizedResidualPzE = new TH1F( "Kaons (refitted Track)" , "; #sqrt{(_{}p_{z}^{REC} - p_{z}^{MC})#times(_{}E^{REC} - E^{MC}) / #sigma_{p_{z}E}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_RFTrk_NormalizedResidualPzE = 0;
	h_RFTrk_NormalizedResidualE = new TH1F( "Kaons (refitted Track)" , "; (_{}E^{REC} - E^{MC}) / #sigma_{E}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_RFTrk_NormalizedResidualE = 0;
	h_RFTrk_NormalizedResidualTheta = new TH1F( "Kaons (refitted Track)" , "; (_{}#theta^{REC} - #theta^{MC}) / #sigma_{#theta}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_RFTrk_NormalizedResidualTheta = 0;
	h_RFTrk_NormalizedResidualPhi = new TH1F( "Kaons (refitted Track)" , "; (_{}#phi^{REC} - #phi^{MC}) / #sigma_{#phi}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_RFTrk_NormalizedResidualPhi = 0;
}

void ChargedPFOCorrection::Clear()
{
	m_pfoCharge.clear();
	m_nTracksOfPFO.clear();
	m_TrkToMCPLinkWeight.clear();
	m_oldPFO_Px.clear();
	m_oldPFO_Py.clear();
	m_oldPFO_Pz.clear();
	m_oldPFO_E.clear();
	m_oldPFO_Theta.clear();
	m_oldPFO_Phi.clear();
	m_newPFO_Px.clear();
	m_newPFO_Py.clear();
	m_newPFO_Pz.clear();
	m_newPFO_E.clear();
	m_newPFO_Theta.clear();
	m_newPFO_Phi.clear();
	m_oldPFO_ResidualPx.clear();
	m_oldPFO_ResidualPy.clear();
	m_oldPFO_ResidualPz.clear();
	m_oldPFO_ResidualE.clear();
	m_oldPFO_ResidualTheta.clear();
	m_oldPFO_ResidualPhi.clear();
	m_newPFO_ResidualPx.clear();
	m_newPFO_ResidualPy.clear();
	m_newPFO_ResidualPz.clear();
	m_newPFO_ResidualE.clear();
	m_newPFO_ResidualTheta.clear();
	m_newPFO_ResidualPhi.clear();
	m_oldPFO_NormalizedResidualPx.clear();
	m_oldPFO_NormalizedResidualPy.clear();
	m_oldPFO_NormalizedResidualPz.clear();
	m_oldPFO_NormalizedResidualE.clear();
	m_oldPFO_NormalizedResidualTheta.clear();
	m_oldPFO_NormalizedResidualPhi.clear();
	m_newPFO_NormalizedResidualPx.clear();
	m_newPFO_NormalizedResidualPy.clear();
	m_newPFO_NormalizedResidualPz.clear();
	m_newPFO_NormalizedResidualE.clear();
	m_newPFO_NormalizedResidualTheta.clear();
	m_newPFO_NormalizedResidualPhi.clear();
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
			std::vector<double> TrackInfo;
			int TrackID = 0;
			double maxweightTRKtoMCP = 0.;
			double maxweightMCPtoTRK = 0.;
			Track *refittedTrack = NULL;
			bool m_updatePFO = true;
			bool foundLinkedMCP = true;
//			float pfoMass = 0.0;
			TLorentzVector pfoFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
			TLorentzVector mcpFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
			std::vector<float> oldResiduals( 6 , 0.0 );
			std::vector<float> newResiduals( 6 , 0.0 );
			std::vector<float> oldPFOCovMat = inputPFO->getCovMatrix();
			std::vector<float> newPFOCovMat( 10, 0.0 );
			std::vector<float> oldAngularUncertainties( 2 , 0.0 );
			std::vector<float> newAngularUncertainties( 2 , 0.0 );
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
			m_pfoCharge.push_back( inputPFO->getCharge() );
			m_nTracksOfPFO.push_back( nTRKsofPFO );
			h_pfoCharge_nTracks->Fill( inputPFO->getCharge() , nTRKsofPFO );

			if ( nTRKsofPFO == 0 )
			{
				m_updatePFO = false;
				foundLinkedMCP = false;
			}
			else if ( nTRKsofPFO == 1 )
			{
				Track *inputTrk = (Track*)inputPFOtrkvec.at(0);
				TrackInfo = this->getTruthTrkID( pLCEvent, inputTrk );
				TrackID = TrackInfo[ 0 ];
				maxweightTRKtoMCP = TrackInfo[ 1 ];
				maxweightMCPtoTRK = TrackInfo[ 2 ];
				int weightTRKtoMCP = maxweightTRKtoMCP * 100;
				int weightMCPtoTRK = maxweightMCPtoTRK * 100;
				int TrackIndex = this->getTrackIndex( MarlinTrkTracks , inputTrk );
				if ( TrackID != 0 )
				{
					streamlog_out(DEBUG1) << "	Track linked to a MCParticle (LinkWeight = " << weightTRKtoMCP << "%)	,	MCParticle linked to same Track (LinkWeight = " << weightMCPtoTRK << "%)" << std::endl;
					mcpFourMomentum = this->getMCPFourMomentum( pLCEvent, inputTrk );
					if ( maxweightTRKtoMCP >= m_MinWeightTrackMCTruthLink && maxweightMCPtoTRK >= m_MinWeightMCTruthTrackLink )
					{
						foundLinkedMCP = true;
					}
					else
					{
						foundLinkedMCP = false;
					}
				}
				else
				{
					streamlog_out(DEBUG1) << "	Couldn't find a MCParticle linked to Track in both direction!" << std::endl;
					foundLinkedMCP = false;
					m_updatePFO = false;
				}
				streamlog_out(DEBUG2) << "	TrackID: 	" << TrackID << std::endl;
				streamlog_out(DEBUG2) << "	TrackIndex: 	" << TrackIndex << std::endl;

				if ( abs( TrackID ) == 2212 && foundLinkedMCP )
				{
					refittedTrack = dynamic_cast<EVENT::Track*>( MarlinTrkTracksPROTON->getElementAt( TrackIndex ) );
					outputPFOtrkvec.push_back( refittedTrack );
					trackMass = m_proton_mass;
					m_updatePFO = m_updatePFOwithOneTrack;
					streamlog_out(DEBUG2) << "	Standard track is replaced with track refitted with proton mass" << std::endl;
				}
				else if ( abs( TrackID ) == 321 && foundLinkedMCP )
				{
					refittedTrack = dynamic_cast<EVENT::Track*>( MarlinTrkTracksKAON->getElementAt( TrackIndex ) );
					outputPFOtrkvec.push_back( refittedTrack );
					trackMass = m_kaon_mass;
					m_updatePFO = m_updatePFOwithOneTrack;
					streamlog_out(DEBUG2) << "	Standard track is replaced with track refitted with kaon mass" << std::endl;
				}
				else
				{
					refittedTrack = dynamic_cast<EVENT::Track*>( MarlinTrkTracks->getElementAt( TrackIndex ) );
					outputPFOtrkvec.push_back( refittedTrack );
					trackMass = m_pion_mass;
					m_updatePFO = m_updatePFOwithPionTrack;
					streamlog_out(DEBUG2) << "	Standard track is used for updating PFO" << std::endl;
				}

				pfoFourMomentum = this->getTrackFourMomentum( refittedTrack , trackMass );
				newPFOCovMat = this->UpdateChargedPFOCovMat( refittedTrack , trackMass );
			}
			else if ( nTRKsofPFO == 2 )
			{
				if ( inputPFO->getCharge() == 0.0 )
				{
					h_InnermostRadiusHit_Neutral->Fill( RadiusOfInnermostHit[ 0 ] , RadiusOfInnermostHit[ 0 ] );
					for ( int i_trk = 0 ; i_trk < nTRKsofPFO ; ++i_trk )
					{
						streamlog_out(DEBUG2) << "	---------------------" << std::endl;
						streamlog_out(DEBUG2) << "	Adding Track " << i_trk << " to PFO " << std::endl;
						streamlog_out(DEBUG2) << "	---------------------" << std::endl;
						Track *inputTrk = (Track*)inputPFOtrkvec.at( i_trk );
						TrackInfo = this->getTruthTrkID( pLCEvent, inputTrk );
						TrackID = TrackInfo[ 0 ];
						maxweightTRKtoMCP = TrackInfo[ 1 ];
						maxweightMCPtoTRK = TrackInfo[ 2 ];
						int weightTRKtoMCP = maxweightTRKtoMCP * 100;
						int weightMCPtoTRK = maxweightMCPtoTRK * 100;
						int TrackIndex = this->getTrackIndex( MarlinTrkTracks , inputTrk );
						if ( TrackID != 0 )
						{
							streamlog_out(DEBUG1) << "	Track linked to a MCParticle (LinkWeight = " << weightTRKtoMCP << "%)	,	MCParticle linked to same Track (LinkWeight = " << weightMCPtoTRK << "%)" << std::endl;
							mcpFourMomentum += this->getMCPFourMomentum( pLCEvent, inputTrk );
							if ( foundLinkedMCP && maxweightTRKtoMCP >= m_MinWeightTrackMCTruthLink && maxweightMCPtoTRK >= m_MinWeightMCTruthTrackLink )
							{
								foundLinkedMCP = true;
							}
							else
							{
								foundLinkedMCP = false;
							}
						}
						else
						{
							streamlog_out(DEBUG1) << "	Couldn't find a MCParticle linked to Track in both direction!" << std::endl;
							if ( foundLinkedMCP ) foundLinkedMCP = false;
						}
						streamlog_out(DEBUG2) << "	TrackID: 	" << TrackID << std::endl;
						streamlog_out(DEBUG2) << "	TrackIndex: 	" << TrackIndex << std::endl;

						if ( abs( TrackID ) == 2212 && foundLinkedMCP )
						{
							refittedTrack = dynamic_cast<EVENT::Track*>( MarlinTrkTracksPROTON->getElementAt( TrackIndex ) );
							outputPFOtrkvec.push_back( refittedTrack );
							trackMass = m_proton_mass;
							m_updatePFO = m_updatePFOwithTwoTrack;
							streamlog_out(DEBUG2) << "	Standard track is replaced with track refitted with proton mass" << std::endl;
						}
						else if ( abs( TrackID ) == 321 && foundLinkedMCP )
						{
							refittedTrack = dynamic_cast<EVENT::Track*>( MarlinTrkTracksKAON->getElementAt( TrackIndex ) );
							outputPFOtrkvec.push_back( refittedTrack );
							trackMass = m_kaon_mass;
							m_updatePFO = m_updatePFOwithTwoTrack;
							streamlog_out(DEBUG2) << "	Standard track is replaced with track refitted with kaon mass" << std::endl;
						}
						else
						{
							refittedTrack = dynamic_cast<EVENT::Track*>( MarlinTrkTracks->getElementAt( TrackIndex ) );
							outputPFOtrkvec.push_back( refittedTrack );
							trackMass = m_pion_mass;
							m_updatePFO = m_updatePFOwithPionTrack;
							streamlog_out(DEBUG2) << "	Standard track is used for updating PFO" << std::endl;
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
					foundLinkedMCP = false;
					m_updatePFO = false;
				}
			}
			else
			{
				foundLinkedMCP = false;
				m_updatePFO = false;
			}
			if( m_updatePFO )
			{
				this->updatePFO( inputPFO , outputPFO , outputPFOtrkvec , pfoFourMomentum , newPFOCovMat );
				m_newPFO_Px.push_back( pfoFourMomentum.Px() );		m_newPFO_Py.push_back( pfoFourMomentum.Py() );		m_newPFO_Pz.push_back( pfoFourMomentum.Pz() );
				m_newPFO_E.push_back( pfoFourMomentum.E() );		m_newPFO_Theta.push_back( pfoFourMomentum.Theta() );	m_newPFO_Phi.push_back( pfoFourMomentum.Phi() );
			}
			else
			{
				this->updatePFO( inputPFO , outputPFO , inputPFOtrkvec , oldpfoFourMomentum , oldPFOCovMat );
				m_newPFO_Px.push_back( oldpfoFourMomentum.Px() );	m_newPFO_Py.push_back( oldpfoFourMomentum.Py() );	m_newPFO_Pz.push_back( oldpfoFourMomentum.Pz() );
				m_newPFO_E.push_back( oldpfoFourMomentum.E() );		m_newPFO_Theta.push_back( oldpfoFourMomentum.Theta() );	m_newPFO_Phi.push_back( oldpfoFourMomentum.Phi() );
			}
			m_oldPFO_Px.push_back( oldpfoFourMomentum.Px() );	m_oldPFO_Py.push_back( oldpfoFourMomentum.Py() );	m_oldPFO_Pz.push_back( oldpfoFourMomentum.Pz() );
			m_oldPFO_E.push_back( oldpfoFourMomentum.E() );		m_oldPFO_Theta.push_back( oldpfoFourMomentum.Theta() );	m_oldPFO_Phi.push_back( oldpfoFourMomentum.Phi() );

			streamlog_out(DEBUG3) << "	------------------------------" << std::endl;
			streamlog_out(DEBUG3) << "	Getting Residuals for inputPFO" << std::endl;
			if ( foundLinkedMCP ) oldResiduals = this->getPFOResidual( oldpfoFourMomentum , mcpFourMomentum );
			streamlog_out(DEBUG3) << "	-------------------------------" << std::endl;
			streamlog_out(DEBUG3) << "	Getting Residuals for outputPFO" << std::endl;
			TLorentzVector newpfoFourMomentum( outputPFO->getMomentum()[ 0 ] , outputPFO->getMomentum()[ 1 ] , outputPFO->getMomentum()[ 2 ] , outputPFO->getEnergy() );
			if ( foundLinkedMCP ) newResiduals = this->getPFOResidual( newpfoFourMomentum , mcpFourMomentum );
			if ( foundLinkedMCP )
			{
				oldAngularUncertainties = this->getAngularUncertainties( oldpfoFourMomentum , oldPFOCovMat );
				newAngularUncertainties = this->getAngularUncertainties( newpfoFourMomentum , newPFOCovMat );

				m_oldPFO_ResidualPx.push_back( oldResiduals[ 0 ] );	m_oldPFO_ResidualPy.push_back( oldResiduals[ 1 ] );	m_oldPFO_ResidualPz.push_back( oldResiduals[ 2 ] );
				m_oldPFO_ResidualE.push_back( oldResiduals[ 3 ] );	m_oldPFO_ResidualTheta.push_back( oldResiduals[ 4 ] );	m_oldPFO_ResidualPhi.push_back( oldResiduals[ 5 ] );
				m_newPFO_ResidualPx.push_back( newResiduals[ 0 ] );	m_newPFO_ResidualPy.push_back( newResiduals[ 1 ] );	m_newPFO_ResidualPz.push_back( newResiduals[ 2 ] );
				m_newPFO_ResidualE.push_back( newResiduals[ 3 ] );	m_newPFO_ResidualTheta.push_back( newResiduals[ 4 ] );	m_newPFO_ResidualPhi.push_back( newResiduals[ 5 ] );
				if ( nTRKsofPFO == 1 )
				{
					m_oldPFO_NormalizedResidualPx.push_back( oldResiduals[ 0 ] / sqrt( oldPFOCovMat[ 0 ] ) );
					m_oldPFO_NormalizedResidualPy.push_back( oldResiduals[ 1 ] / sqrt( oldPFOCovMat[ 2 ] ) );
					m_oldPFO_NormalizedResidualPz.push_back( oldResiduals[ 2 ] / sqrt( oldPFOCovMat[ 5 ] ) );
					m_oldPFO_NormalizedResidualE.push_back( oldResiduals[ 3 ] / sqrt( oldPFOCovMat[ 9 ] ) );
					m_oldPFO_NormalizedResidualTheta.push_back( oldResiduals[ 4 ] / oldAngularUncertainties[ 0 ] );
					m_oldPFO_NormalizedResidualPhi.push_back( oldResiduals[ 5 ] / oldAngularUncertainties[ 1 ] );
					h_OldPFOS1trk_NormalizedResidualPx->Fill( oldResiduals[ 0 ] / sqrt( oldPFOCovMat[ 0 ] ) ); ++n_OldPFOS1trk_NormalizedResidualPx;
					int OldPFO_SignNormResi1 = ( oldResiduals[ 0 ] * oldResiduals[ 1 ] / oldPFOCovMat[ 1 ] >= 0 ? 1 : -1 );
					h_OldPFOS1trk_NormalizedResidualPxPy->Fill( OldPFO_SignNormResi1 * sqrt( fabs( oldResiduals[ 0 ] ) * fabs( oldResiduals[ 1 ] ) / fabs( oldPFOCovMat[ 1 ] ) ) ); ++n_OldPFOS1trk_NormalizedResidualPxPy;
					h_OldPFOS1trk_NormalizedResidualPy->Fill( oldResiduals[ 1 ] / sqrt( oldPFOCovMat[ 2 ] ) ); ++n_OldPFOS1trk_NormalizedResidualPy;
					int OldPFO_SignNormResi3 = ( oldResiduals[ 0 ] * oldResiduals[ 2 ] / oldPFOCovMat[ 3 ] >= 0 ? 1 : -1 );
					h_OldPFOS1trk_NormalizedResidualPxPz->Fill( OldPFO_SignNormResi3 * sqrt( fabs( oldResiduals[ 0 ] ) * fabs( oldResiduals[ 2 ] ) / fabs( oldPFOCovMat[ 3 ] ) ) ); ++n_OldPFOS1trk_NormalizedResidualPxPz;
					int OldPFO_SignNormResi4 = ( oldResiduals[ 1 ] * oldResiduals[ 2 ] / oldPFOCovMat[ 4 ] >= 0 ? 1 : -1 );
					h_OldPFOS1trk_NormalizedResidualPyPz->Fill( OldPFO_SignNormResi4 * sqrt( fabs( oldResiduals[ 1 ] ) * fabs( oldResiduals[ 2 ] ) / fabs( oldPFOCovMat[ 4 ] ) ) ); ++n_OldPFOS1trk_NormalizedResidualPyPz;
					h_OldPFOS1trk_NormalizedResidualPz->Fill( oldResiduals[ 2 ] / sqrt( oldPFOCovMat[ 5 ] ) ); ++n_OldPFOS1trk_NormalizedResidualPz;
					int OldPFO_SignNormResi6 = ( oldResiduals[ 0 ] * oldResiduals[ 3 ] / oldPFOCovMat[ 6 ] >= 0 ? 1 : -1 );
					h_OldPFOS1trk_NormalizedResidualPxE->Fill( OldPFO_SignNormResi6 * sqrt( fabs( oldResiduals[ 0 ] ) * fabs( oldResiduals[ 3 ] ) / fabs( oldPFOCovMat[ 6 ] ) ) ); ++n_OldPFOS1trk_NormalizedResidualPxE;
					int OldPFO_SignNormResi7 = ( oldResiduals[ 1 ] * oldResiduals[ 3 ] / oldPFOCovMat[ 7 ] >= 0 ? 1 : -1 );
					h_OldPFOS1trk_NormalizedResidualPyE->Fill( OldPFO_SignNormResi7 * sqrt( fabs( oldResiduals[ 1 ] ) * fabs( oldResiduals[ 3 ] ) / fabs( oldPFOCovMat[ 7 ] ) ) ); ++n_OldPFOS1trk_NormalizedResidualPyE;
					int OldPFO_SignNormResi8 = ( oldResiduals[ 2 ] * oldResiduals[ 3 ] / oldPFOCovMat[ 6 ] >= 0 ? 1 : -1 );
					h_OldPFOS1trk_NormalizedResidualPzE->Fill( OldPFO_SignNormResi8 * sqrt( fabs( oldResiduals[ 2 ] ) * fabs( oldResiduals[ 3 ] ) / fabs( oldPFOCovMat[ 8 ] ) ) ); ++n_OldPFOS1trk_NormalizedResidualPzE;
					h_OldPFOS1trk_NormalizedResidualE->Fill( oldResiduals[ 3 ] / sqrt( oldPFOCovMat[ 9 ] ) ); ++n_OldPFOS1trk_NormalizedResidualE;
					h_OldPFOS1trk_NormalizedResidualTheta->Fill( oldResiduals[ 4 ] / oldAngularUncertainties[ 0 ] ); ++n_OldPFOS1trk_NormalizedResidualTheta;
					h_OldPFOS1trk_NormalizedResidualPhi->Fill( oldResiduals[ 5 ] / oldAngularUncertainties[ 1 ] ); ++n_OldPFOS1trk_NormalizedResidualPhi;
					m_newPFO_NormalizedResidualPx.push_back( newResiduals[ 0 ] / sqrt( newPFOCovMat[ 0 ] ) );
					m_newPFO_NormalizedResidualPy.push_back( newResiduals[ 1 ] / sqrt( newPFOCovMat[ 2 ] ) );
					m_newPFO_NormalizedResidualPz.push_back( newResiduals[ 2 ] / sqrt( newPFOCovMat[ 5 ] ) );
					m_newPFO_NormalizedResidualE.push_back( newResiduals[ 3 ] / sqrt( newPFOCovMat[ 9 ] ) );
					m_newPFO_NormalizedResidualTheta.push_back( newResiduals[ 4 ] / newAngularUncertainties[ 0 ] );
					m_newPFO_NormalizedResidualPhi.push_back( newResiduals[ 5 ] / newAngularUncertainties[ 1 ] );
					h_NewPFOS1trk_NormalizedResidualPx->Fill( newResiduals[ 0 ] / sqrt( newPFOCovMat[ 0 ] ) ); ++n_NewPFOS1trk_NormalizedResidualPx;
					int NewPFO_SignNormResi1 = ( newResiduals[ 0 ] * newResiduals[ 1 ] / newPFOCovMat[ 1 ] >= 0 ? 1 : -1 );
					h_NewPFOS1trk_NormalizedResidualPxPy->Fill( NewPFO_SignNormResi1 * sqrt( fabs( newResiduals[ 0 ] ) * fabs( newResiduals[ 1 ] ) / fabs( newPFOCovMat[ 1 ] ) ) ); ++n_NewPFOS1trk_NormalizedResidualPxPy;
					h_NewPFOS1trk_NormalizedResidualPy->Fill( newResiduals[ 1 ] / sqrt( newPFOCovMat[ 2 ] ) ); ++n_NewPFOS1trk_NormalizedResidualPy;
					int NewPFO_SignNormResi3 = ( newResiduals[ 0 ] * newResiduals[ 2 ] / newPFOCovMat[ 3 ] >= 0 ? 1 : -1 );
					h_NewPFOS1trk_NormalizedResidualPxPz->Fill( NewPFO_SignNormResi3 * sqrt( fabs( newResiduals[ 0 ] ) * fabs( newResiduals[ 2 ] ) / fabs( newPFOCovMat[ 3 ] ) ) ); ++n_NewPFOS1trk_NormalizedResidualPxPz;
					int NewPFO_SignNormResi4 = ( newResiduals[ 1 ] * newResiduals[ 2 ] / newPFOCovMat[ 4 ] >= 0 ? 1 : -1 );
					h_NewPFOS1trk_NormalizedResidualPyPz->Fill( NewPFO_SignNormResi4 * sqrt( fabs( newResiduals[ 1 ] ) * fabs( newResiduals[ 2 ] ) / fabs( newPFOCovMat[ 4 ] ) ) ); ++n_NewPFOS1trk_NormalizedResidualPyPz;
					h_NewPFOS1trk_NormalizedResidualPz->Fill( newResiduals[ 2 ] / sqrt( newPFOCovMat[ 5 ] ) ); ++n_NewPFOS1trk_NormalizedResidualPz;
					int NewPFO_SignNormResi6 = ( newResiduals[ 0 ] * newResiduals[ 3 ] / newPFOCovMat[ 6 ] >= 0 ? 1 : -1 );
					h_NewPFOS1trk_NormalizedResidualPxE->Fill( NewPFO_SignNormResi6 * sqrt( fabs( newResiduals[ 0 ] ) * fabs( newResiduals[ 3 ] ) / fabs( newPFOCovMat[ 6 ] ) ) ); ++n_NewPFOS1trk_NormalizedResidualPxE;
					int NewPFO_SignNormResi7 = ( newResiduals[ 1 ] * newResiduals[ 3 ] / newPFOCovMat[ 7 ] >= 0 ? 1 : -1 );
					h_NewPFOS1trk_NormalizedResidualPyE->Fill( NewPFO_SignNormResi7 * sqrt( fabs( newResiduals[ 1 ] ) * fabs( newResiduals[ 3 ] ) / fabs( newPFOCovMat[ 7 ] ) ) ); ++n_NewPFOS1trk_NormalizedResidualPyE;
					int NewPFO_SignNormResi8 = ( newResiduals[ 2 ] * newResiduals[ 3 ] / newPFOCovMat[ 6 ] >= 0 ? 1 : -1 );
					h_NewPFOS1trk_NormalizedResidualPzE->Fill( NewPFO_SignNormResi8 * sqrt( fabs( newResiduals[ 2 ] ) * fabs( newResiduals[ 3 ] ) / fabs( newPFOCovMat[ 8 ] ) ) ); ++n_NewPFOS1trk_NormalizedResidualPzE;
					h_NewPFOS1trk_NormalizedResidualE->Fill( newResiduals[ 3 ] / sqrt( newPFOCovMat[ 9 ] ) ); ++n_NewPFOS1trk_NormalizedResidualE;
					h_NewPFOS1trk_NormalizedResidualTheta->Fill( newResiduals[ 4 ] / newAngularUncertainties[ 0 ] ); ++n_NewPFOS1trk_NormalizedResidualTheta;
					h_NewPFOS1trk_NormalizedResidualPhi->Fill( newResiduals[ 5 ] / newAngularUncertainties[ 1 ] ); ++n_NewPFOS1trk_NormalizedResidualPhi;
					if ( abs( TrackID ) == 2212 )
					{
						h_StdTrkProton_NormalizedResidualPx->Fill( oldResiduals[ 0 ] / sqrt( oldPFOCovMat[ 0 ] ) ); ++n_StdTrkProton_NormalizedResidualPx;
						h_StdTrkProton_NormalizedResidualPxPy->Fill( OldPFO_SignNormResi1 * sqrt( fabs( oldResiduals[ 0 ] ) * fabs( oldResiduals[ 1 ] ) / fabs( oldPFOCovMat[ 1 ] ) ) ); ++n_StdTrkProton_NormalizedResidualPxPy;
						h_StdTrkProton_NormalizedResidualPy->Fill( oldResiduals[ 1 ] / sqrt( oldPFOCovMat[ 2 ] ) ); ++n_StdTrkProton_NormalizedResidualPy;
						h_StdTrkProton_NormalizedResidualPxPz->Fill( OldPFO_SignNormResi3 * sqrt( fabs( oldResiduals[ 0 ] ) * fabs( oldResiduals[ 2 ] ) / fabs( oldPFOCovMat[ 3 ] ) ) ); ++n_StdTrkProton_NormalizedResidualPxPz;
						h_StdTrkProton_NormalizedResidualPyPz->Fill( OldPFO_SignNormResi4 * sqrt( fabs( oldResiduals[ 1 ] ) * fabs( oldResiduals[ 2 ] ) / fabs( oldPFOCovMat[ 4 ] ) ) ); ++n_StdTrkProton_NormalizedResidualPyPz;
						h_StdTrkProton_NormalizedResidualPz->Fill( oldResiduals[ 2 ] / sqrt( oldPFOCovMat[ 5 ] ) ); ++n_StdTrkProton_NormalizedResidualPz;
						h_StdTrkProton_NormalizedResidualPxE->Fill( OldPFO_SignNormResi6 * sqrt( fabs( oldResiduals[ 0 ] ) * fabs( oldResiduals[ 3 ] ) / fabs( oldPFOCovMat[ 6 ] ) ) ); ++n_StdTrkProton_NormalizedResidualPxE;
						h_StdTrkProton_NormalizedResidualPyE->Fill( OldPFO_SignNormResi7 * sqrt( fabs( oldResiduals[ 1 ] ) * fabs( oldResiduals[ 3 ] ) / fabs( oldPFOCovMat[ 7 ] ) ) ); ++n_StdTrkProton_NormalizedResidualPyE;
						h_StdTrkProton_NormalizedResidualPzE->Fill( OldPFO_SignNormResi8 * sqrt( fabs( oldResiduals[ 2 ] ) * fabs( oldResiduals[ 3 ] ) / fabs( oldPFOCovMat[ 8 ] ) ) ); ++n_StdTrkProton_NormalizedResidualPzE;
						h_StdTrkProton_NormalizedResidualE->Fill( oldResiduals[ 3 ] / sqrt( oldPFOCovMat[ 9 ] ) ); ++n_StdTrkProton_NormalizedResidualE;
						h_StdTrkProton_NormalizedResidualTheta->Fill( oldResiduals[ 4 ] / oldAngularUncertainties[ 0 ] ); ++n_StdTrkProton_NormalizedResidualTheta;
						h_StdTrkProton_NormalizedResidualPhi->Fill( oldResiduals[ 5 ] / oldAngularUncertainties[ 1 ] ); ++n_StdTrkProton_NormalizedResidualPhi;
						h_RFTrkProton_NormalizedResidualPx->Fill( newResiduals[ 0 ] / sqrt( newPFOCovMat[ 0 ] ) ); ++n_RFTrkProton_NormalizedResidualPx;
						h_RFTrkProton_NormalizedResidualPxPy->Fill( NewPFO_SignNormResi1 * sqrt( fabs( newResiduals[ 0 ] ) * fabs( newResiduals[ 1 ] ) / fabs( newPFOCovMat[ 1 ] ) ) ); ++n_RFTrkProton_NormalizedResidualPxPy;
						h_RFTrkProton_NormalizedResidualPy->Fill( newResiduals[ 1 ] / sqrt( newPFOCovMat[ 2 ] ) ); ++n_RFTrkProton_NormalizedResidualPy;
						h_RFTrkProton_NormalizedResidualPxPz->Fill( NewPFO_SignNormResi3 * sqrt( fabs( newResiduals[ 0 ] ) * fabs( newResiduals[ 2 ] ) / fabs( newPFOCovMat[ 3 ] ) ) ); ++n_RFTrkProton_NormalizedResidualPxPz;
						h_RFTrkProton_NormalizedResidualPyPz->Fill( NewPFO_SignNormResi4 * sqrt( fabs( newResiduals[ 1 ] ) * fabs( newResiduals[ 2 ] ) / fabs( newPFOCovMat[ 4 ] ) ) ); ++n_RFTrkProton_NormalizedResidualPyPz;
						h_RFTrkProton_NormalizedResidualPz->Fill( newResiduals[ 2 ] / sqrt( newPFOCovMat[ 5 ] ) ); ++n_RFTrkProton_NormalizedResidualPz;
						h_RFTrkProton_NormalizedResidualPxE->Fill( NewPFO_SignNormResi6 * sqrt( fabs( newResiduals[ 0 ] ) * fabs( newResiduals[ 3 ] ) / fabs( newPFOCovMat[ 6 ] ) ) ); ++n_RFTrkProton_NormalizedResidualPxE;
						h_RFTrkProton_NormalizedResidualPyE->Fill( NewPFO_SignNormResi7 * sqrt( fabs( newResiduals[ 1 ] ) * fabs( newResiduals[ 3 ] ) / fabs( newPFOCovMat[ 7 ] ) ) ); ++n_RFTrkProton_NormalizedResidualPyE;
						h_RFTrkProton_NormalizedResidualPzE->Fill( NewPFO_SignNormResi8 * sqrt( fabs( newResiduals[ 2 ] ) * fabs( newResiduals[ 3 ] ) / fabs( newPFOCovMat[ 8 ] ) ) ); ++n_RFTrkProton_NormalizedResidualPzE;
						h_RFTrkProton_NormalizedResidualE->Fill( newResiduals[ 3 ] / sqrt( newPFOCovMat[ 9 ] ) ); ++n_RFTrkProton_NormalizedResidualE;
						h_RFTrkProton_NormalizedResidualTheta->Fill( newResiduals[ 4 ] / newAngularUncertainties[ 0 ] ); ++n_RFTrkProton_NormalizedResidualTheta;
						h_RFTrkProton_NormalizedResidualPhi->Fill( newResiduals[ 5 ] / newAngularUncertainties[ 1 ] ); ++n_RFTrkProton_NormalizedResidualPhi;
					}
					else if ( abs( TrackID ) == 321 )
					{
						h_StdTrkKaon_NormalizedResidualPx->Fill( oldResiduals[ 0 ] / sqrt( oldPFOCovMat[ 0 ] ) ); ++n_StdTrkKaon_NormalizedResidualPx;
						h_StdTrkKaon_NormalizedResidualPxPy->Fill( OldPFO_SignNormResi1 * sqrt( fabs( oldResiduals[ 0 ] ) * fabs( oldResiduals[ 1 ] ) / fabs( oldPFOCovMat[ 1 ] ) ) ); ++n_StdTrkKaon_NormalizedResidualPxPy;
						h_StdTrkKaon_NormalizedResidualPy->Fill( oldResiduals[ 1 ] / sqrt( oldPFOCovMat[ 2 ] ) ); ++n_StdTrkKaon_NormalizedResidualPy;
						h_StdTrkKaon_NormalizedResidualPxPz->Fill( OldPFO_SignNormResi3 * sqrt( fabs( oldResiduals[ 0 ] ) * fabs( oldResiduals[ 2 ] ) / fabs( oldPFOCovMat[ 3 ] ) ) ); ++n_StdTrkKaon_NormalizedResidualPxPz;
						h_StdTrkKaon_NormalizedResidualPyPz->Fill( OldPFO_SignNormResi4 * sqrt( fabs( oldResiduals[ 1 ] ) * fabs( oldResiduals[ 2 ] ) / fabs( oldPFOCovMat[ 4 ] ) ) ); ++n_StdTrkKaon_NormalizedResidualPyPz;
						h_StdTrkKaon_NormalizedResidualPz->Fill( oldResiduals[ 2 ] / sqrt( oldPFOCovMat[ 5 ] ) ); ++n_StdTrkKaon_NormalizedResidualPz;
						h_StdTrkKaon_NormalizedResidualPxE->Fill( OldPFO_SignNormResi6 * sqrt( fabs( oldResiduals[ 0 ] ) * fabs( oldResiduals[ 3 ] ) / fabs( oldPFOCovMat[ 6 ] ) ) ); ++n_StdTrkKaon_NormalizedResidualPxE;
						h_StdTrkKaon_NormalizedResidualPyE->Fill( OldPFO_SignNormResi7 * sqrt( fabs( oldResiduals[ 1 ] ) * fabs( oldResiduals[ 3 ] ) / fabs( oldPFOCovMat[ 7 ] ) ) ); ++n_StdTrkKaon_NormalizedResidualPyE;
						h_StdTrkKaon_NormalizedResidualPzE->Fill( OldPFO_SignNormResi8 * sqrt( fabs( oldResiduals[ 2 ] ) * fabs( oldResiduals[ 3 ] ) / fabs( oldPFOCovMat[ 8 ] ) ) ); ++n_StdTrkKaon_NormalizedResidualPzE;
						h_StdTrkKaon_NormalizedResidualE->Fill( oldResiduals[ 3 ] / sqrt( oldPFOCovMat[ 9 ] ) ); ++n_StdTrkKaon_NormalizedResidualE;
						h_StdTrkKaon_NormalizedResidualTheta->Fill( oldResiduals[ 4 ] / oldAngularUncertainties[ 0 ] ); ++n_StdTrkKaon_NormalizedResidualTheta;
						h_StdTrkKaon_NormalizedResidualPhi->Fill( oldResiduals[ 5 ] / oldAngularUncertainties[ 1 ] ); ++n_StdTrkKaon_NormalizedResidualPhi;
						h_RFTrkKaon_NormalizedResidualPx->Fill( newResiduals[ 0 ] / sqrt( newPFOCovMat[ 0 ] ) ); ++n_RFTrkKaon_NormalizedResidualPx;
						h_RFTrkKaon_NormalizedResidualPxPy->Fill( NewPFO_SignNormResi1 * sqrt( fabs( newResiduals[ 0 ] ) * fabs( newResiduals[ 1 ] ) / fabs( newPFOCovMat[ 1 ] ) ) ); ++n_RFTrkKaon_NormalizedResidualPxPy;
						h_RFTrkKaon_NormalizedResidualPy->Fill( newResiduals[ 1 ] / sqrt( newPFOCovMat[ 2 ] ) ); ++n_RFTrkKaon_NormalizedResidualPy;
						h_RFTrkKaon_NormalizedResidualPxPz->Fill( NewPFO_SignNormResi3 * sqrt( fabs( newResiduals[ 0 ] ) * fabs( newResiduals[ 2 ] ) / fabs( newPFOCovMat[ 3 ] ) ) ); ++n_RFTrkKaon_NormalizedResidualPxPz;
						h_RFTrkKaon_NormalizedResidualPyPz->Fill( NewPFO_SignNormResi4 * sqrt( fabs( newResiduals[ 1 ] ) * fabs( newResiduals[ 2 ] ) / fabs( newPFOCovMat[ 4 ] ) ) ); ++n_RFTrkKaon_NormalizedResidualPyPz;
						h_RFTrkKaon_NormalizedResidualPz->Fill( newResiduals[ 2 ] / sqrt( newPFOCovMat[ 5 ] ) ); ++n_RFTrkKaon_NormalizedResidualPz;
						h_RFTrkKaon_NormalizedResidualPxE->Fill( NewPFO_SignNormResi6 * sqrt( fabs( newResiduals[ 0 ] ) * fabs( newResiduals[ 3 ] ) / fabs( newPFOCovMat[ 6 ] ) ) ); ++n_RFTrkKaon_NormalizedResidualPxE;
						h_RFTrkKaon_NormalizedResidualPyE->Fill( NewPFO_SignNormResi7 * sqrt( fabs( newResiduals[ 1 ] ) * fabs( newResiduals[ 3 ] ) / fabs( newPFOCovMat[ 7 ] ) ) ); ++n_RFTrkKaon_NormalizedResidualPyE;
						h_RFTrkKaon_NormalizedResidualPzE->Fill( NewPFO_SignNormResi8 * sqrt( fabs( newResiduals[ 2 ] ) * fabs( newResiduals[ 3 ] ) / fabs( newPFOCovMat[ 8 ] ) ) ); ++n_RFTrkKaon_NormalizedResidualPzE;
						h_RFTrkKaon_NormalizedResidualE->Fill( newResiduals[ 3 ] / sqrt( newPFOCovMat[ 9 ] ) ); ++n_RFTrkKaon_NormalizedResidualE;
						h_RFTrkKaon_NormalizedResidualTheta->Fill( newResiduals[ 4 ] / newAngularUncertainties[ 0 ] ); ++n_RFTrkKaon_NormalizedResidualTheta;
						h_RFTrkKaon_NormalizedResidualPhi->Fill( newResiduals[ 5 ] / newAngularUncertainties[ 1 ] ); ++n_RFTrkKaon_NormalizedResidualPhi;
					}
				}
				else if ( nTRKsofPFO == 2 )
				{
					m_oldPFO_NormalizedResidualPx.push_back( 0 );
					m_oldPFO_NormalizedResidualPy.push_back( 0 );
					m_oldPFO_NormalizedResidualPz.push_back( 0 );
					m_oldPFO_NormalizedResidualE.push_back( 0 );
					m_oldPFO_NormalizedResidualTheta.push_back( 0 );
					m_oldPFO_NormalizedResidualPhi.push_back( 0 );
					m_newPFO_NormalizedResidualPx.push_back( newResiduals[ 0 ] / sqrt( newPFOCovMat[ 0 ] ) );
					m_newPFO_NormalizedResidualPy.push_back( newResiduals[ 1 ] / sqrt( newPFOCovMat[ 2 ] ) );
					m_newPFO_NormalizedResidualPz.push_back( newResiduals[ 2 ] / sqrt( newPFOCovMat[ 5 ] ) );
					m_newPFO_NormalizedResidualE.push_back( newResiduals[ 3 ] / sqrt( newPFOCovMat[ 9 ] ) );
					m_newPFO_NormalizedResidualTheta.push_back( newResiduals[ 4 ] / newAngularUncertainties[ 0 ] );
					m_newPFO_NormalizedResidualPhi.push_back( newResiduals[ 5 ] / newAngularUncertainties[ 1 ] );
					h_NewPFOS2trk_NormalizedResidualPx->Fill( newResiduals[ 0 ] / sqrt( newPFOCovMat[ 0 ] ) ); ++n_NewPFOS2trk_NormalizedResidualPx;
					int NewPFO_SignNormResi1 = ( newResiduals[ 0 ] * newResiduals[ 1 ] / newPFOCovMat[ 1 ] >= 0 ? 1 : -1 );
					h_NewPFOS2trk_NormalizedResidualPxPy->Fill( NewPFO_SignNormResi1 * sqrt( fabs( newResiduals[ 0 ] ) * fabs( newResiduals[ 1 ] ) / fabs( newPFOCovMat[ 1 ] ) ) ); ++n_NewPFOS2trk_NormalizedResidualPxPy;
					h_NewPFOS2trk_NormalizedResidualPy->Fill( newResiduals[ 1 ] / sqrt( newPFOCovMat[ 2 ] ) ); ++n_NewPFOS2trk_NormalizedResidualPy;
					int NewPFO_SignNormResi3 = ( newResiduals[ 0 ] * newResiduals[ 2 ] / newPFOCovMat[ 3 ] >= 0 ? 1 : -1 );
					h_NewPFOS2trk_NormalizedResidualPxPz->Fill( NewPFO_SignNormResi3 * sqrt( fabs( newResiduals[ 0 ] ) * fabs( newResiduals[ 2 ] ) / fabs( newPFOCovMat[ 3 ] ) ) ); ++n_NewPFOS2trk_NormalizedResidualPxPz;
					int NewPFO_SignNormResi4 = ( newResiduals[ 1 ] * newResiduals[ 2 ] / newPFOCovMat[ 4 ] >= 0 ? 1 : -1 );
					h_NewPFOS2trk_NormalizedResidualPyPz->Fill( NewPFO_SignNormResi4 * sqrt( fabs( newResiduals[ 1 ] ) * fabs( newResiduals[ 2 ] ) / fabs( newPFOCovMat[ 4 ] ) ) ); ++n_NewPFOS2trk_NormalizedResidualPyPz;
					h_NewPFOS2trk_NormalizedResidualPz->Fill( newResiduals[ 2 ] / sqrt( newPFOCovMat[ 5 ] ) ); ++n_NewPFOS2trk_NormalizedResidualPz;
					int NewPFO_SignNormResi6 = ( newResiduals[ 0 ] * newResiduals[ 3 ] / newPFOCovMat[ 6 ] >= 0 ? 1 : -1 );
					h_NewPFOS2trk_NormalizedResidualPxE->Fill( NewPFO_SignNormResi6 * sqrt( fabs( newResiduals[ 0 ] ) * fabs( newResiduals[ 3 ] ) / fabs( newPFOCovMat[ 6 ] ) ) ); ++n_NewPFOS2trk_NormalizedResidualPxE;
					int NewPFO_SignNormResi7 = ( newResiduals[ 1 ] * newResiduals[ 3 ] / newPFOCovMat[ 7 ] >= 0 ? 1 : -1 );
					h_NewPFOS2trk_NormalizedResidualPyE->Fill( NewPFO_SignNormResi7 * sqrt( fabs( newResiduals[ 1 ] ) * fabs( newResiduals[ 3 ] ) / fabs( newPFOCovMat[ 7 ] ) ) ); ++n_NewPFOS2trk_NormalizedResidualPyE;
					int NewPFO_SignNormResi8 = ( newResiduals[ 2 ] * newResiduals[ 3 ] / newPFOCovMat[ 6 ] >= 0 ? 1 : -1 );
					h_NewPFOS2trk_NormalizedResidualPzE->Fill( NewPFO_SignNormResi8 * sqrt( fabs( newResiduals[ 2 ] ) * fabs( newResiduals[ 3 ] ) / fabs( newPFOCovMat[ 8 ] ) ) ); ++n_NewPFOS2trk_NormalizedResidualPzE;
					h_NewPFOS2trk_NormalizedResidualE->Fill( newResiduals[ 3 ] / sqrt( newPFOCovMat[ 9 ] ) ); ++n_NewPFOS2trk_NormalizedResidualE;
					h_NewPFOS2trk_NormalizedResidualTheta->Fill( newResiduals[ 4 ] / newAngularUncertainties[ 0 ] ); ++n_NewPFOS2trk_NormalizedResidualTheta;
					h_NewPFOS2trk_NormalizedResidualPhi->Fill( newResiduals[ 5 ] / newAngularUncertainties[ 1 ] ); ++n_NewPFOS2trk_NormalizedResidualPhi;
					for ( int i_trk = 0 ; i_trk < nTRKsofPFO ; ++i_trk)
					{
						Track *inputTrk = (Track*)inputPFOtrkvec.at( i_trk );
						Track *outputTrk = (Track*)outputPFOtrkvec.at( i_trk );
						TLorentzVector inputTrack4Momentum( 0.0 , 0.0 , 0.0 , 0.0 );
						TLorentzVector outputTrack4Momentum( 0.0 , 0.0 , 0.0 , 0.0 );
						TLorentzVector mcp4Momentum( 0.0 , 0.0 , 0.0 , 0.0 );
						std::vector<float> inputTrackResiduals( 6 , 0.0 );
						std::vector<float> outputTrackResiduals( 6 , 0.0 );
						std::vector<float> outputTrackUncertainties( 2 , 0.0 );
						float outputTrackMass = m_pion_mass;
						float inputTrackMass = m_pion_mass;
						if ( getTrackIndex( MarlinTrkTracksPROTON , outputTrk ) != -1 ) outputTrackMass = m_proton_mass;
						if ( getTrackIndex( MarlinTrkTracksKAON , outputTrk ) != -1 ) outputTrackMass = m_kaon_mass;
						if ( getTrackIndex( MarlinTrkTracks , outputTrk ) != -1  ) outputTrackMass = m_pion_mass;
						mcp4Momentum = this->getMCPFourMomentum( pLCEvent, inputTrk );
						inputTrack4Momentum = this->getTrackFourMomentum( inputTrk , inputTrackMass );
						outputTrack4Momentum = this->getTrackFourMomentum( outputTrk , outputTrackMass );
						inputTrackResiduals = this->getPFOResidual( inputTrack4Momentum , mcp4Momentum );
						outputTrackResiduals = this->getPFOResidual( outputTrack4Momentum , mcp4Momentum );
						std::vector<float> inputTrkCovMat = this->UpdateChargedPFOCovMat( inputTrk , inputTrackMass );
						std::vector<float> inputTrackAngularUncertainties = this->getAngularUncertainties( inputTrack4Momentum , inputTrkCovMat );
						std::vector<float> outputTrkCovMat = this->UpdateChargedPFOCovMat( outputTrk , outputTrackMass );
						std::vector<float> outputTrackAngularUncertainties = this->getAngularUncertainties( outputTrack4Momentum , outputTrkCovMat );
						h_StdTrk_NormalizedResidualPx->Fill( inputTrackResiduals[ 0 ] / inputTrkCovMat[ 0 ] ); ++n_StdTrk_NormalizedResidualPx;
						int StdTrk_SignNormResi1 = ( inputTrackResiduals[ 0 ] * inputTrackResiduals[ 1 ] / inputTrkCovMat[ 1 ] >= 0 ? 1 : -1 );
						h_StdTrk_NormalizedResidualPxPy->Fill( StdTrk_SignNormResi1 * sqrt( fabs( inputTrackResiduals[ 0 ] ) * fabs( inputTrackResiduals[ 1 ] ) / fabs( inputTrkCovMat[ 1 ] ) ) ); ++n_StdTrk_NormalizedResidualPxPy;
						h_StdTrk_NormalizedResidualPy->Fill( inputTrackResiduals[ 1 ] / inputTrkCovMat[ 2 ] ); ++n_StdTrk_NormalizedResidualPy;
						int StdTrk_SignNormResi3 = ( inputTrackResiduals[ 0 ] * inputTrackResiduals[ 2 ] / inputTrkCovMat[ 3 ] >= 0 ? 1 : -1 );
						h_StdTrk_NormalizedResidualPxPz->Fill( StdTrk_SignNormResi3 * sqrt( fabs( inputTrackResiduals[ 0 ] ) * fabs( inputTrackResiduals[ 2 ] ) / fabs( inputTrkCovMat[ 3 ] ) ) ); ++n_StdTrk_NormalizedResidualPxPz;
						int StdTrk_SignNormResi4 = ( inputTrackResiduals[ 1 ] * inputTrackResiduals[ 2 ] / inputTrkCovMat[ 4 ] >= 0 ? 1 : -1 );
						h_StdTrk_NormalizedResidualPyPz->Fill( StdTrk_SignNormResi4 * sqrt( fabs( inputTrackResiduals[ 1 ] ) * fabs( inputTrackResiduals[ 2 ] ) / fabs( inputTrkCovMat[ 4 ] ) ) ); ++n_StdTrk_NormalizedResidualPyPz;
						h_StdTrk_NormalizedResidualPz->Fill( inputTrackResiduals[ 2 ] / inputTrkCovMat[ 5 ] ); ++n_StdTrk_NormalizedResidualPz;
						int StdTrk_SignNormResi6 = ( inputTrackResiduals[ 0 ] * inputTrackResiduals[ 3 ] / inputTrkCovMat[ 6 ] >= 0 ? 1 : -1 );
						h_StdTrk_NormalizedResidualPxE->Fill( StdTrk_SignNormResi6 * sqrt( fabs( inputTrackResiduals[ 0 ] ) * fabs( inputTrackResiduals[ 3 ] ) / fabs( inputTrkCovMat[ 6 ] ) ) ); ++n_StdTrk_NormalizedResidualPxE;
						int StdTrk_SignNormResi7 = ( inputTrackResiduals[ 1 ] * inputTrackResiduals[ 3 ] / inputTrkCovMat[ 7 ] >= 0 ? 1 : -1 );
						h_StdTrk_NormalizedResidualPyE->Fill( StdTrk_SignNormResi7 * sqrt( fabs( inputTrackResiduals[ 1 ] ) * fabs( inputTrackResiduals[ 3 ] ) / fabs( inputTrkCovMat[ 7 ] ) ) ); ++n_StdTrk_NormalizedResidualPyE;
						int StdTrk_SignNormResi8 = ( inputTrackResiduals[ 2 ] * inputTrackResiduals[ 3 ] / inputTrkCovMat[ 8 ] >= 0 ? 1 : -1 );
						h_StdTrk_NormalizedResidualPzE->Fill( StdTrk_SignNormResi8 * sqrt( fabs( inputTrackResiduals[ 2 ] ) * fabs( inputTrackResiduals[ 3 ] ) / fabs( inputTrkCovMat[ 8 ] ) ) ); ++n_StdTrk_NormalizedResidualPzE;
						h_StdTrk_NormalizedResidualE->Fill( inputTrackResiduals[ 3 ] / inputTrkCovMat[ 9 ] ); ++n_StdTrk_NormalizedResidualE;
						h_StdTrk_NormalizedResidualTheta->Fill( inputTrackResiduals[ 4 ] / inputTrackAngularUncertainties[ 0 ] ); ++n_StdTrk_NormalizedResidualTheta;
						h_StdTrk_NormalizedResidualPhi->Fill( inputTrackResiduals[ 5 ] / inputTrackAngularUncertainties[ 1 ] ); ++n_StdTrk_NormalizedResidualPhi;
						h_RFTrk_NormalizedResidualPx->Fill( outputTrackResiduals[ 0 ] / outputTrkCovMat[ 0 ] ); ++n_RFTrk_NormalizedResidualPx;
						int RFTrk_SignNormResi1 = ( outputTrackResiduals[ 0 ] * outputTrackResiduals[ 1 ] / outputTrkCovMat[ 1 ] >= 0 ? 1 : -1 );
						h_RFTrk_NormalizedResidualPxPy->Fill( RFTrk_SignNormResi1 * sqrt( fabs( outputTrackResiduals[ 0 ] ) * fabs( outputTrackResiduals[ 1 ] ) / fabs( outputTrkCovMat[ 1 ] ) ) ); ++n_RFTrk_NormalizedResidualPxPy;
						h_RFTrk_NormalizedResidualPy->Fill( outputTrackResiduals[ 1 ] / outputTrkCovMat[ 2 ] ); ++n_RFTrk_NormalizedResidualPy;
						int RFTrk_SignNormResi3 = ( outputTrackResiduals[ 0 ] * outputTrackResiduals[ 2 ] / outputTrkCovMat[ 3 ] >= 0 ? 1 : -1 );
						h_RFTrk_NormalizedResidualPxPz->Fill( RFTrk_SignNormResi3 * sqrt( fabs( outputTrackResiduals[ 0 ] ) * fabs( outputTrackResiduals[ 2 ] ) / fabs( outputTrkCovMat[ 3 ] ) ) ); ++n_RFTrk_NormalizedResidualPxPz;
						int RFTrk_SignNormResi4 = ( outputTrackResiduals[ 1 ] * outputTrackResiduals[ 2 ] / outputTrkCovMat[ 4 ] >= 0 ? 1 : -1 );
						h_RFTrk_NormalizedResidualPyPz->Fill( RFTrk_SignNormResi4 * sqrt( fabs( outputTrackResiduals[ 1 ] ) * fabs( outputTrackResiduals[ 2 ] ) / fabs( outputTrkCovMat[ 4 ] ) ) ); ++n_RFTrk_NormalizedResidualPyPz;
						h_RFTrk_NormalizedResidualPz->Fill( outputTrackResiduals[ 2 ] / outputTrkCovMat[ 5 ] ); ++n_RFTrk_NormalizedResidualPz;
						int RFTrk_SignNormResi6 = ( outputTrackResiduals[ 0 ] * outputTrackResiduals[ 3 ] / outputTrkCovMat[ 6 ] >= 0 ? 1 : -1 );
						h_RFTrk_NormalizedResidualPxE->Fill( RFTrk_SignNormResi6 * sqrt( fabs( outputTrackResiduals[ 0 ] ) * fabs( outputTrackResiduals[ 3 ] ) / fabs( outputTrkCovMat[ 6 ] ) ) ); ++n_RFTrk_NormalizedResidualPxE;
						int RFTrk_SignNormResi7 = ( outputTrackResiduals[ 1 ] * outputTrackResiduals[ 3 ] / outputTrkCovMat[ 7 ] >= 0 ? 1 : -1 );
						h_RFTrk_NormalizedResidualPyE->Fill( RFTrk_SignNormResi7 * sqrt( fabs( outputTrackResiduals[ 1 ] ) * fabs( outputTrackResiduals[ 3 ] ) / fabs( outputTrkCovMat[ 7 ] ) ) ); ++n_RFTrk_NormalizedResidualPyE;
						int RFTrk_SignNormResi8 = ( outputTrackResiduals[ 2 ] * outputTrackResiduals[ 3 ] / outputTrkCovMat[ 8 ] >= 0 ? 1 : -1 );
						h_RFTrk_NormalizedResidualPzE->Fill( RFTrk_SignNormResi8 * sqrt( fabs( outputTrackResiduals[ 2 ] ) * fabs( outputTrackResiduals[ 3 ] ) / fabs( outputTrkCovMat[ 8 ] ) ) ); ++n_RFTrk_NormalizedResidualPzE;
						h_RFTrk_NormalizedResidualE->Fill( outputTrackResiduals[ 3 ] / outputTrkCovMat[ 9 ] ); ++n_RFTrk_NormalizedResidualE;
						h_RFTrk_NormalizedResidualTheta->Fill( outputTrackResiduals[ 4 ] / outputTrackAngularUncertainties[ 0 ] ); ++n_RFTrk_NormalizedResidualTheta;
						h_RFTrk_NormalizedResidualPhi->Fill( outputTrackResiduals[ 5 ] / outputTrackAngularUncertainties[ 1 ] ); ++n_RFTrk_NormalizedResidualPhi;
					}
				}
			}
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

std::vector<double> ChargedPFOCorrection::getTruthTrkID( EVENT::LCEvent *pLCEvent , EVENT::Track* inputTrk )
{
	std::vector<float> TRKMCPlink;
	std::vector<double> TrackInfo;
	TrackInfo.push_back( 0 );
	TrackInfo.push_back( 0 );
	TrackInfo.push_back( 0 );
	MCParticle *linkedMCP{};
	try
	{
		LCRelationNavigator TrackMCParticleNav( pLCEvent->getCollection( m_TrackMCTruthLinkCollection ) );
		LCRelationNavigator MCParticleTrackNav( pLCEvent->getCollection( m_MCTruthTrackLinkCollection ) );
	}
	catch (DataNotAvailableException &e)
	{
		streamlog_out(WARNING) << "Could not find the LCRelationNavigator for Track <-> MCParticle" << std::endl;
		return TrackInfo;
	}
	LCRelationNavigator TrackMCParticleNav( pLCEvent->getCollection( m_TrackMCTruthLinkCollection ) );
	const EVENT::LCObjectVec& mcpvec = TrackMCParticleNav.getRelatedToObjects(inputTrk);
	const EVENT::FloatVec&  mcpweightvec = TrackMCParticleNav.getRelatedToWeights(inputTrk);
	int TrackID = 0;
	double maxweightTRKtoMCP = 0.;
	double maxweightMCPtoTRK = 0.;
	int iTRKtoMCPmax = -1;
	for ( unsigned int i_mcp = 0; i_mcp < mcpvec.size(); i_mcp++ )
	{
		double mcp_weight = mcpweightvec.at(i_mcp);
		MCParticle *testMCP = (MCParticle *) mcpvec.at(i_mcp);
		if ( mcp_weight > maxweightTRKtoMCP )//&& mcp_weight >= m_MinWeightTrackMCTruthLink )
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
		int iMCPtoTRKmax = -1;
		for ( unsigned int i_trk = 0; i_trk < trkvec.size(); i_trk++ )
		{
			double track_weight = trkweightvec.at( i_trk );
			if ( track_weight > maxweightMCPtoTRK )//&& track_weight >= m_MinWeightTrackMCTruthLink )
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
//				TrackInfo.push_back( TrackID );
//				TrackInfo.push_back( maxweightTRKtoMCP );
//				TrackInfo.push_back( maxweightMCPtoTRK );
				TrackInfo[ 0 ] = TrackID;
				TrackInfo[ 1 ] = maxweightTRKtoMCP;
				TrackInfo[ 2 ] = maxweightMCPtoTRK;
			}
		}
	}
	if ( TrackID != 0 ) m_TrkToMCPLinkWeight.push_back( maxweightTRKtoMCP );
	return TrackInfo;
}

int ChargedPFOCorrection::getTrackIndex( EVENT::LCCollection *TrackCollection , EVENT::Track* inputTrk )
{
	LCCollection *MarlinTrkTracks{};
	try
	{
		MarlinTrkTracks = TrackCollection;//pLCEvent->getCollection( TrackCollection );
	}
	catch (DataNotAvailableException &e)
	{
		streamlog_out(WARNING) << "	Could not find the " << TrackCollection << " Collection" << std::endl;
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
			if ( track_weight > maxweightMCPtoTRK && track_weight >= m_MinWeightMCTruthTrackLink )
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

std::vector<float> ChargedPFOCorrection::getAngularUncertainties( TLorentzVector pfoFourMomentum , std::vector<float> pfoCovMat )
{
	std::vector<float> AngularUncertainties{};
	float Px	= pfoFourMomentum.Px();
	float Py	= pfoFourMomentum.Py();
	float Pz	= pfoFourMomentum.Pz();
//	float E		= pfoFourMomentum.E();
	float P2	= pow( Px , 2 ) + pow( Py , 2 ) + pow( Pz , 2 );
	float Pt2	= pow( Px , 2 ) + pow( Py , 2 );
	float Pt	= sqrt( pow( Px , 2 ) + pow( Py , 2 ) );
	float SigPx2	= pfoCovMat[0];
	float SigPxPy	= pfoCovMat[1];
	float SigPy2	= pfoCovMat[2];
	float SigPxPz	= pfoCovMat[3];
	float SigPyPz	= pfoCovMat[4];
	float SigPz2	= pfoCovMat[5];
//	float SigPxE	= pfoCovMat[6];
//	float SigPyE	= pfoCovMat[7];
//	float SigPzE	= pfoCovMat[8];
//	float SigE2	= pfoCovMat[9];

	float dth_dpx	= Px * Pz / ( P2 * Pt );
	float dth_dpy	= Py * Pz / ( P2 * Pt );
	float dth_dpz	= -Pt / P2;

	float dphi_dpx	= -Py / Pt2;
	float dphi_dpy	= Px / Pt2;

	float SigmaTheta= std::sqrt( std::fabs( SigPx2 * std::pow( dth_dpx , 2 ) + SigPy2 * std::pow( dth_dpy , 2 ) + SigPz2 * std::pow( dth_dpz , 2 ) + 2 * ( SigPxPy * dth_dpx * dth_dpy ) + 2 * ( SigPyPz * dth_dpy * dth_dpz ) + 2 * ( SigPxPz * dth_dpx * dth_dpz ) ) );
	float SigmaPhi	= std::sqrt( std::fabs( SigPx2 * std::pow( dphi_dpx , 2 ) + SigPy2 * std::pow( dphi_dpy , 2 ) + 2 * ( SigPxPy * dphi_dpx * dphi_dpy ) ) );

	AngularUncertainties.push_back( SigmaTheta );
	AngularUncertainties.push_back( SigmaPhi );
	return AngularUncertainties;
}

void ChargedPFOCorrection::InitializeHistogram( TH1F *histogram , int scale , int color , int lineWidth , int markerSize , int markerStyle )
{
	histogram->Scale( 1.0 / scale );
	histogram->SetLineColor( color );
	histogram->SetLineWidth( lineWidth );
	histogram->SetMarkerSize( markerSize );
	histogram->SetMarkerStyle( markerStyle );
	histogram->SetMarkerColor( color );
	float fit_range = 2.0;
	float fit_min = -2.0;
	float fit_max = 2.0;
	for ( int i_fit = 0 ; i_fit < 5 ; ++i_fit )
	{
		histogram->Fit( "gaus" , "" , "" , fit_min , fit_max );
		TF1 *fitFunction = (TF1 *)histogram->GetFunction("gaus");
		double fitMean = fitFunction->GetParameter( 1 );
		double fitSigma = fitFunction->GetParameter( 2 );
		fit_min = fitMean - fit_range * fitSigma;
		fit_max = fitMean + fit_range * fitSigma;
	}
	histogram->GetFunction("gaus")->SetLineColor( color );
	float y_max = 1.2 * histogram->GetMaximum();
	histogram->GetYaxis()->SetRangeUser(0.0, y_max);
	histogram->Write();
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

	if ( m_fillRootTree )
	{
		m_pTFile->cd();
		m_pTTree->Write();
		m_Histograms->cd();
		h_nClusters_nTracks->Write();
		h_pfoCharge_nTracks->Write();
		h_InnermostRadiusHit_Neutral->Write();
		h_InnermostRadiusHit_Charged->Write();
		h_FirstSubDet_Charged->Write();
		m_OldPFOs_1Trk->cd();
		streamlog_out(MESSAGE) << "	Fill Histograms for PFOs with ONE track using standard tracks" << std::endl;
		InitializeHistogram( h_OldPFOS1trk_NormalizedResidualPx , n_OldPFOS1trk_NormalizedResidualPx , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_OldPFOS1trk_NormalizedResidualPxPy , n_OldPFOS1trk_NormalizedResidualPxPy , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_OldPFOS1trk_NormalizedResidualPy , n_OldPFOS1trk_NormalizedResidualPy , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_OldPFOS1trk_NormalizedResidualPxPz , n_OldPFOS1trk_NormalizedResidualPxPz , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_OldPFOS1trk_NormalizedResidualPyPz , n_OldPFOS1trk_NormalizedResidualPyPz , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_OldPFOS1trk_NormalizedResidualPz , n_OldPFOS1trk_NormalizedResidualPz , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_OldPFOS1trk_NormalizedResidualPxE , n_OldPFOS1trk_NormalizedResidualPxE , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_OldPFOS1trk_NormalizedResidualPyE , n_OldPFOS1trk_NormalizedResidualPyE , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_OldPFOS1trk_NormalizedResidualPzE , n_OldPFOS1trk_NormalizedResidualPzE , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_OldPFOS1trk_NormalizedResidualE , n_OldPFOS1trk_NormalizedResidualE , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_OldPFOS1trk_NormalizedResidualTheta , n_OldPFOS1trk_NormalizedResidualTheta , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_OldPFOS1trk_NormalizedResidualPhi , n_OldPFOS1trk_NormalizedResidualPhi , 4 , 1 , 1.0 , 1 );
		streamlog_out(MESSAGE) << "	Histograms for PFOs with ONE track have been filled using standard tracks" << std::endl;
		m_NewPFOs_1Trk->cd();
		streamlog_out(MESSAGE) << "	Fill Histograms for PFOs with ONE track using tracks refitted with true mass" << std::endl;
		InitializeHistogram( h_NewPFOS1trk_NormalizedResidualPx , n_NewPFOS1trk_NormalizedResidualPx , 2 , 1 , 1.0 , 1 );
		InitializeHistogram( h_NewPFOS1trk_NormalizedResidualPxPy , n_NewPFOS1trk_NormalizedResidualPxPy , 2 , 1 , 1.0 , 1 );
		InitializeHistogram( h_NewPFOS1trk_NormalizedResidualPy , n_NewPFOS1trk_NormalizedResidualPy , 2 , 1 , 1.0 , 1 );
		InitializeHistogram( h_NewPFOS1trk_NormalizedResidualPxPz , n_NewPFOS1trk_NormalizedResidualPxPz , 2 , 1 , 1.0 , 1 );
		InitializeHistogram( h_NewPFOS1trk_NormalizedResidualPyPz , n_NewPFOS1trk_NormalizedResidualPyPz , 2 , 1 , 1.0 , 1 );
		InitializeHistogram( h_NewPFOS1trk_NormalizedResidualPz , n_NewPFOS1trk_NormalizedResidualPz , 2 , 1 , 1.0 , 1 );
		InitializeHistogram( h_NewPFOS1trk_NormalizedResidualPxE , n_NewPFOS1trk_NormalizedResidualPxE , 2 , 1 , 1.0 , 1 );
		InitializeHistogram( h_NewPFOS1trk_NormalizedResidualPyE , n_NewPFOS1trk_NormalizedResidualPyE , 2 , 1 , 1.0 , 1 );
		InitializeHistogram( h_NewPFOS1trk_NormalizedResidualPzE , n_NewPFOS1trk_NormalizedResidualPzE , 2 , 1 , 1.0 , 1 );
		InitializeHistogram( h_NewPFOS1trk_NormalizedResidualE , n_NewPFOS1trk_NormalizedResidualE , 2 , 1 , 1.0 , 1 );
		InitializeHistogram( h_NewPFOS1trk_NormalizedResidualTheta , n_NewPFOS1trk_NormalizedResidualTheta , 2 , 1 , 1.0 , 1 );
		InitializeHistogram( h_NewPFOS1trk_NormalizedResidualPhi , n_NewPFOS1trk_NormalizedResidualPhi , 2 , 1 , 1.0 , 1 );
		streamlog_out(MESSAGE) << "	Histograms for PFOs with ONE track have been filled using tracks refitted with true mass" << std::endl;
		m_NewPFOs_2Trk->cd();
		streamlog_out(MESSAGE) << "	Fill Histograms for PFOs with TWO tracks using tracks refitted with true mass" << std::endl;
		InitializeHistogram( h_NewPFOS2trk_NormalizedResidualPx , n_NewPFOS2trk_NormalizedResidualPx , 6 , 1 , 1.0 , 1 );
		InitializeHistogram( h_NewPFOS2trk_NormalizedResidualPxPy , n_NewPFOS2trk_NormalizedResidualPxPy , 6 , 1 , 1.0 , 1 );
		InitializeHistogram( h_NewPFOS2trk_NormalizedResidualPy , n_NewPFOS2trk_NormalizedResidualPy , 6 , 1 , 1.0 , 1 );
		InitializeHistogram( h_NewPFOS2trk_NormalizedResidualPxPz , n_NewPFOS2trk_NormalizedResidualPxPz , 6 , 1 , 1.0 , 1 );
		InitializeHistogram( h_NewPFOS2trk_NormalizedResidualPyPz , n_NewPFOS2trk_NormalizedResidualPyPz , 6 , 1 , 1.0 , 1 );
		InitializeHistogram( h_NewPFOS2trk_NormalizedResidualPz , n_NewPFOS2trk_NormalizedResidualPz , 6 , 1 , 1.0 , 1 );
		InitializeHistogram( h_NewPFOS2trk_NormalizedResidualPxE , n_NewPFOS2trk_NormalizedResidualPxE , 6 , 1 , 1.0 , 1 );
		InitializeHistogram( h_NewPFOS2trk_NormalizedResidualPyE , n_NewPFOS2trk_NormalizedResidualPyE , 6 , 1 , 1.0 , 1 );
		InitializeHistogram( h_NewPFOS2trk_NormalizedResidualPzE , n_NewPFOS2trk_NormalizedResidualPzE , 6 , 1 , 1.0 , 1 );
		InitializeHistogram( h_NewPFOS2trk_NormalizedResidualE , n_NewPFOS2trk_NormalizedResidualE , 6 , 1 , 1.0 , 1 );
		InitializeHistogram( h_NewPFOS2trk_NormalizedResidualTheta , n_NewPFOS2trk_NormalizedResidualTheta , 6 , 1 , 1.0 , 1 );
		InitializeHistogram( h_NewPFOS2trk_NormalizedResidualPhi , n_NewPFOS2trk_NormalizedResidualPhi , 6 , 1 , 1.0 , 1 );
		streamlog_out(MESSAGE) << "	Histograms for PFOs with ONE track have been filled using tracks refitted with true mass" << std::endl;
		m_TrueProtons_1Trk->cd();
		streamlog_out(MESSAGE) << "	Fill Histograms for protons in PFOs with ONE tracks" << std::endl;
		InitializeHistogram( h_StdTrkProton_NormalizedResidualPx , n_StdTrkProton_NormalizedResidualPx , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_StdTrkProton_NormalizedResidualPxPy , n_StdTrkProton_NormalizedResidualPxPy , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_StdTrkProton_NormalizedResidualPy , n_StdTrkProton_NormalizedResidualPy , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_StdTrkProton_NormalizedResidualPxPz , n_StdTrkProton_NormalizedResidualPxPz , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_StdTrkProton_NormalizedResidualPyPz , n_StdTrkProton_NormalizedResidualPyPz , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_StdTrkProton_NormalizedResidualPz , n_StdTrkProton_NormalizedResidualPz , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_StdTrkProton_NormalizedResidualPxE , n_StdTrkProton_NormalizedResidualPxE , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_StdTrkProton_NormalizedResidualPyE , n_StdTrkProton_NormalizedResidualPyE , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_StdTrkProton_NormalizedResidualPzE , n_StdTrkProton_NormalizedResidualPzE , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_StdTrkProton_NormalizedResidualE , n_StdTrkProton_NormalizedResidualE , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_StdTrkProton_NormalizedResidualTheta , n_StdTrkProton_NormalizedResidualTheta , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_StdTrkProton_NormalizedResidualPhi , n_StdTrkProton_NormalizedResidualPhi , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_RFTrkProton_NormalizedResidualPx , n_RFTrkProton_NormalizedResidualPx , 2 , 1 , 1.0 , 1 );
		InitializeHistogram( h_RFTrkProton_NormalizedResidualPxPy , n_RFTrkProton_NormalizedResidualPxPy , 2 , 1 , 1.0 , 1 );
		InitializeHistogram( h_RFTrkProton_NormalizedResidualPy , n_RFTrkProton_NormalizedResidualPy , 2 , 1 , 1.0 , 1 );
		InitializeHistogram( h_RFTrkProton_NormalizedResidualPxPz , n_RFTrkProton_NormalizedResidualPxPz , 2 , 1 , 1.0 , 1 );
		InitializeHistogram( h_RFTrkProton_NormalizedResidualPyPz , n_RFTrkProton_NormalizedResidualPyPz , 2 , 1 , 1.0 , 1 );
		InitializeHistogram( h_RFTrkProton_NormalizedResidualPz , n_RFTrkProton_NormalizedResidualPz , 2 , 1 , 1.0 , 1 );
		InitializeHistogram( h_RFTrkProton_NormalizedResidualPxE , n_RFTrkProton_NormalizedResidualPxE , 2 , 1 , 1.0 , 1 );
		InitializeHistogram( h_RFTrkProton_NormalizedResidualPyE , n_RFTrkProton_NormalizedResidualPyE , 2 , 1 , 1.0 , 1 );
		InitializeHistogram( h_RFTrkProton_NormalizedResidualPzE , n_RFTrkProton_NormalizedResidualPzE , 2 , 1 , 1.0 , 1 );
		InitializeHistogram( h_RFTrkProton_NormalizedResidualE , n_RFTrkProton_NormalizedResidualE , 2 , 1 , 1.0 , 1 );
		InitializeHistogram( h_RFTrkProton_NormalizedResidualTheta , n_RFTrkProton_NormalizedResidualTheta , 2 , 1 , 1.0 , 1 );
		InitializeHistogram( h_RFTrkProton_NormalizedResidualPhi , n_RFTrkProton_NormalizedResidualPhi , 2 , 1 , 1.0 , 1 );
		streamlog_out(MESSAGE) << "	Histograms for protons in PFOs with ONE track have been filled" << std::endl;
		m_TrueKaons_1Trk->cd();
		streamlog_out(MESSAGE) << "	Fill Histograms for kaons in PFOs with ONE tracks" << std::endl;
		InitializeHistogram( h_StdTrkKaon_NormalizedResidualPx , n_StdTrkKaon_NormalizedResidualPx , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_StdTrkKaon_NormalizedResidualPxPy , n_StdTrkKaon_NormalizedResidualPxPy , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_StdTrkKaon_NormalizedResidualPy , n_StdTrkKaon_NormalizedResidualPy , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_StdTrkKaon_NormalizedResidualPxPz , n_StdTrkKaon_NormalizedResidualPxPz , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_StdTrkKaon_NormalizedResidualPyPz , n_StdTrkKaon_NormalizedResidualPyPz , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_StdTrkKaon_NormalizedResidualPz , n_StdTrkKaon_NormalizedResidualPz , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_StdTrkKaon_NormalizedResidualPxE , n_StdTrkKaon_NormalizedResidualPxE , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_StdTrkKaon_NormalizedResidualPyE , n_StdTrkKaon_NormalizedResidualPyE , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_StdTrkKaon_NormalizedResidualPzE , n_StdTrkKaon_NormalizedResidualPzE , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_StdTrkKaon_NormalizedResidualE , n_StdTrkKaon_NormalizedResidualE , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_StdTrkKaon_NormalizedResidualTheta , n_StdTrkKaon_NormalizedResidualTheta , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_StdTrkKaon_NormalizedResidualPhi , n_StdTrkKaon_NormalizedResidualPhi , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_RFTrkKaon_NormalizedResidualPx , n_RFTrkKaon_NormalizedResidualPx , 2 , 1 , 1.0 , 1 );
		InitializeHistogram( h_RFTrkKaon_NormalizedResidualPxPy , n_RFTrkKaon_NormalizedResidualPxPy , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_RFTrkKaon_NormalizedResidualPy , n_RFTrkKaon_NormalizedResidualPy , 2 , 1 , 1.0 , 1 );
		InitializeHistogram( h_RFTrkKaon_NormalizedResidualPxPz , n_RFTrkKaon_NormalizedResidualPxPz , 2 , 1 , 1.0 , 1 );
		InitializeHistogram( h_RFTrkKaon_NormalizedResidualPyPz , n_RFTrkKaon_NormalizedResidualPyPz , 2 , 1 , 1.0 , 1 );
		InitializeHistogram( h_RFTrkKaon_NormalizedResidualPz , n_RFTrkKaon_NormalizedResidualPz , 2 , 1 , 1.0 , 1 );
		InitializeHistogram( h_RFTrkKaon_NormalizedResidualPxE , n_RFTrkKaon_NormalizedResidualPxE , 2 , 1 , 1.0 , 1 );
		InitializeHistogram( h_RFTrkKaon_NormalizedResidualPyE , n_RFTrkKaon_NormalizedResidualPyE , 2 , 1 , 1.0 , 1 );
		InitializeHistogram( h_RFTrkKaon_NormalizedResidualPzE , n_RFTrkKaon_NormalizedResidualPzE , 2 , 1 , 1.0 , 1 );
		InitializeHistogram( h_RFTrkKaon_NormalizedResidualE , n_RFTrkKaon_NormalizedResidualE , 2 , 1 , 1.0 , 1 );
		InitializeHistogram( h_RFTrkKaon_NormalizedResidualTheta , n_RFTrkKaon_NormalizedResidualTheta , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_RFTrkKaon_NormalizedResidualPhi , n_RFTrkKaon_NormalizedResidualPhi , 2 , 1 , 1.0 , 1 );
		streamlog_out(MESSAGE) << "	Histograms for kaons in PFOs with ONE track have been filled" << std::endl;
		m_IndividualTracks_2Trk->cd();
		streamlog_out(MESSAGE) << "	Fill Histograms for individual tracks in PFOs with TWO tracks" << std::endl;
		InitializeHistogram( h_StdTrk_NormalizedResidualPx , n_StdTrk_NormalizedResidualPx , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_StdTrk_NormalizedResidualPxPy , n_StdTrk_NormalizedResidualPxPy , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_StdTrk_NormalizedResidualPx , n_StdTrk_NormalizedResidualPx , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_StdTrk_NormalizedResidualPxPz , n_StdTrk_NormalizedResidualPxPz , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_StdTrk_NormalizedResidualPyPz , n_StdTrk_NormalizedResidualPyPz , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_StdTrk_NormalizedResidualPx , n_StdTrk_NormalizedResidualPz , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_StdTrk_NormalizedResidualPxE , n_StdTrk_NormalizedResidualPxE , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_StdTrk_NormalizedResidualPyE , n_StdTrk_NormalizedResidualPyE , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_StdTrk_NormalizedResidualPzE , n_StdTrk_NormalizedResidualPzE , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_StdTrk_NormalizedResidualE , n_StdTrk_NormalizedResidualE , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_StdTrk_NormalizedResidualTheta , n_StdTrk_NormalizedResidualTheta , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_StdTrk_NormalizedResidualPhi , n_StdTrk_NormalizedResidualPhi , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_RFTrk_NormalizedResidualPx , n_RFTrk_NormalizedResidualPx , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_RFTrk_NormalizedResidualPxPy , n_RFTrk_NormalizedResidualPxPy , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_RFTrk_NormalizedResidualPx , n_RFTrk_NormalizedResidualPx , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_RFTrk_NormalizedResidualPxPz , n_RFTrk_NormalizedResidualPxPz , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_RFTrk_NormalizedResidualPyPz , n_RFTrk_NormalizedResidualPyPz , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_RFTrk_NormalizedResidualPx , n_RFTrk_NormalizedResidualPz , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_RFTrk_NormalizedResidualPxE , n_RFTrk_NormalizedResidualPxE , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_RFTrk_NormalizedResidualPyE , n_RFTrk_NormalizedResidualPyE , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_RFTrk_NormalizedResidualPzE , n_RFTrk_NormalizedResidualPzE , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_RFTrk_NormalizedResidualE , n_RFTrk_NormalizedResidualE , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_RFTrk_NormalizedResidualTheta , n_RFTrk_NormalizedResidualTheta , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_RFTrk_NormalizedResidualPhi , n_RFTrk_NormalizedResidualPhi , 4 , 1 , 1.0 , 1 );
		streamlog_out(MESSAGE) << "	Histograms for individual tracks in PFOs with TWO tracks have been filled" << std::endl;

		m_pTFile->Close();
		delete m_pTFile;
	}

	streamlog_out(MESSAGE) << " " << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////	processed events: 	" << m_nEvtSum << "	////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	streamlog_out(MESSAGE) << " " << std::endl;

}
