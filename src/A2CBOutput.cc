#include "A2CBOutput.hh"
#include "A2FileGenerator.hh"
#include "G4RunManager.hh"
#include "CLHEP/Units/SystemOfUnits.h"
#include "A2UserTrackInformation.hh"
#include "G4String.hh"

using namespace CLHEP;

A2CBOutput::A2CBOutput()
{

  fFile=NULL;
  fTree=NULL;
  fPGA=const_cast<A2PrimaryGeneratorAction*>(static_cast<const A2PrimaryGeneratorAction*>(G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction()));

  fDET=const_cast<A2DetectorConstruction*>(static_cast<const A2DetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction()));

  //Need to get the number of initial particles
  fnpart=fPGA->GetNGenMaxParticles();
  if(fnpart>100){
    G4cerr<<"A2CBOutput::A2CBOutput() can only store 100 particles. Need to edit fdircos[100][3] in A2CBOutput.hh if you want more particles."<<G4endl;
    exit(1);
  }

  //Get the LorentzVectors of the initial particles
  fGenLorentzVec=(fPGA->GetGenLorentzVecs()); //Only exists if ntuple input
  fBeamLorentzVec=fPGA->GetBeamLorentzVec();//Will take the default beam if no ntuple
  fGenPartType=fPGA->GetGenPartType();
  fvertex=fPGA->GetVertex();
  Int_t* fpart = fPGA->GetParticlesToBeTracked();
  A2UserTrackInformation *tracks;
  //G4Track* ftrack=tracks->GetTrack(fpart);
  //Initialise arrays dependent on number of particles
  //Note can't initialize 2 dimesional arrays like this when writing 
  //to a ROOT TTree, you only get the correct info for the first particle
  //fdircos now declared as a 2D array with 100 and 3 dimensions 
  //fdircos=new Float_t*[fnpart];
  //for(Int_t i=0;i<fnpart;i++) fdircos[i]=new Float_t[3];
  felab=new Float_t[fnpart]; 
  fklab=new Float_t[fnpart];
  fplab=new Float_t[fnpart]; 
  fidpart=new Int_t[fnpart]; 

  //create Tree
  fTree=new TTree("h12","Crystals");
  fTree->SetAutoSave();

  //Collection IDs for decoding hit collections
  fCBCollID = -1;
  fTAPSCollID = -1;

  // store IDs of primary particles
  fStorePrimaries = true;

  //TOF stuff
  if(fDET->GetNToFbars()>0){
    fToFTot=fDET->GetNToFbars();
    ftofi=new Int_t[fToFTot];
    ftofe=new Float_t[fToFTot];
    ftoft=new Float_t[fToFTot];
    ftofx=new Float_t[fToFTot];
    ftofy=new Float_t[fToFTot];
    ftofz=new Float_t[fToFTot];
  }

  fIsGiBUU = false;
  if (fPGA->GetFileGen())
    fIsGiBUU = (fPGA->GetFileGen()->GetType() == A2FileGenerator::kGiBUU);
  fweight = 1;
}

A2CBOutput::~A2CBOutput()
{
	delete fidpart;
	delete fplab;
	delete felab;
	delete fklab;
	if(fTree)delete fTree;

}

void A2CBOutput::SetBranches()
{

	if(!fTree)
	{
		G4cout<<"A2CBOutput::SetBranches() Can't set branches have to set fTree first!"<<G4endl;
		return;
	}
	Int_t basket =64000;

	// Beam, TAPS, and CB
	fTree->Branch("nhits",&fnhits,"fnhits/I",basket);
	fTree->Branch("npart",&fnpart,"fnpart/I",basket);
	fTree->Branch("ntaps",&fntaps,"fntaps/I",basket);
	fTree->Branch("nvtaps",&fnvtaps,"fnvtaps/I",basket);
	fTree->Branch("vhits",&fvhits,"fvhits/I",basket);
	fTree->Branch("plab",fplab,"fplab[fnpart]/F",basket);
	fTree->Branch("tctaps",ftctaps,"ftctaps[fntaps]/F",basket);
	fTree->Branch("vertex",fvertex,"fvertex[3]/F",basket);
	fTree->Branch("beam",fbeam,"fbeam[5]/F",basket);
	fTree->Branch("dircos",fdircos,"fdircos[fnpart][3]/F",basket);
	fTree->Branch("ecryst",fecryst,"fecryst[fnhits]/F",basket);
	fTree->Branch("tcryst",ftcryst,"ftcryst[fnhits]/F",basket);
	fTree->Branch("ectapfs",fectapfs,"fectapfs[fntaps]/F",basket);
	fTree->Branch("ectapsl",fectapsl,"fectapsl[fntaps]/F",basket);
	fTree->Branch("elab",felab,"felab[fnpart]/F",basket);
	fTree->Branch("klab",fklab,"fklab[fnpart]/F",basket);
	//fTree->Branch("theta",ftheta,"ftheta[fnpart]/F",basket);
	fTree->Branch("eleak",&feleak,"feleak/F",basket);
	fTree->Branch("enai",&fenai,"fenai/F",basket);
	fTree->Branch("etot",&fetot,"fetot/F",basket);
	fTree->Branch("eveto",feveto,"feveto[fvhits]/F",basket);
	fTree->Branch("tveto",ftveto,"ftveto[fvhits]/F",basket);
	fTree->Branch("evtaps",fevtaps,"fevtaps[fnvtaps]/F",basket);
	fTree->Branch("icryst",ficryst,"ficryst[fnhits]/I",basket);
	fTree->Branch("ictaps",fictaps,"fictaps[fntaps]/I",basket);

	// TAPS Veto
	fTree->Branch("ivtaps",fivtaps,"fictaps[fnvtaps]/I",basket);
	fTree->Branch("idpart",fidpart,"fidpart[fnpart]/I",basket);
	fTree->Branch("iveto",fiveto,"fiveto[fvhits]/I",basket);

	// MWPC
	if (fDET->GetUseMWPC() && fDET->GetUseMWPC() / 10 == 0)
	{
		fTree->Branch("nmwpc",&fnmwpc,"fnmwpc/I",basket);
		fTree->Branch("imwpc",fimwpc,"fimwpc[fnmwpc]/I",basket);
		fTree->Branch("mposx",fmposx,"fmposx[fnmwpc]/F",basket);
		fTree->Branch("mposy",fmposy,"fmposy[fnmwpc]/F",basket);
		fTree->Branch("mposz",fmposz,"fmposz[fnmwpc]/F",basket);
		fTree->Branch("emwpc",femwpc,"femwpc[fnmwpc]/F",basket);
	}
	else G4cout<<"A2CBOutput::SetBranches() Disabling MWPC readout"<<G4endl;
 
	// TOF
	if(fToFTot>0)
	{
		fTree->Branch("ntof",&fntof,"fntof/I",basket);
		fTree->Branch("tofi",ftofi,"ftofi[fntof]/I",basket);
		fTree->Branch("tofe",ftofe,"ftofe[fntof]/F",basket);
		fTree->Branch("toft",ftoft,"ftoft[fntof]/F",basket);
		fTree->Branch("tofx",ftofx,"ftofx[fntof]/F",basket);
		fTree->Branch("tofy",ftofy,"ftofy[fntof]/F",basket);
		fTree->Branch("tofz",ftofz,"ftofz[fntof]/F",basket);
	}

	// Active He3
	fTree->Branch("nhe3",&fnhe3,"fnhe3/I",basket);
	fTree->Branch("ihe3",fihe3,"fihe3[fnhe3]/I",basket);
	fTree->Branch("ehe3",fehe3,"fehe3[fnhe3]/F",basket);
	fTree->Branch("the3",fthe3,"fthe3[fnhe3]/F",basket);

	// TPC
	fTree->Branch("ntpc",&fntpc,"fntpc/I",basket);
	fTree->Branch("itpc",fitpc,"fitpc[fntpc]/I",basket);
	fTree->Branch("qtpc",fqtpc,"fqtpc[fntpc]/F",basket);
	fTree->Branch("ttpc",fttpc,"fttpc[fntpc]/F",basket);

	// Pizza
	fTree->Branch("npiz",&fnpiz,"fnpiz/I",basket);
	fTree->Branch("ipiz",fipiz,"fipiz[fnpiz]/I",basket);
	fTree->Branch("epiz",fepiz,"fepiz[fnpiz]/F",basket);
	fTree->Branch("tpiz",ftpiz,"ftpiz[fnpiz]/F",basket);

	// CATS
	if ( fDET->GetUseCATS()) 
	{
		// Core
		fTree->Branch( "nCATSCore", &fnCATSCore, "fnCATSCore/I", basket);
		fTree->Branch( "iCATSCore", fiCATSCore, "fiCATSCore[fnCATSCore]/I", basket);
		fTree->Branch( "eCATSCore", feCATSCore, "feCATSCore[fnCATSCore]/F", basket);
		fTree->Branch( "tCATSCore", ftCATSCore, "ftCATSCore[fnCATSCore]/F", basket);

		// Annulus
		fTree->Branch( "nCATSAnn", &fnCATSAnn, "fnCATSAnn/I", basket);
		fTree->Branch( "iCATSAnn", fiCATSAnn, "fiCATSAnn[fnCATSAnn]/I", basket);
		fTree->Branch( "eCATSAnn", feCATSAnn, "feCATSAnn[fnCATSAnn]/F", basket);
		fTree->Branch( "tCATSAnn", ftCATSAnn, "ftCATSAnn[fnCATSAnn]/F", basket);

		// Cosmic Shield
		fTree->Branch( "nCATSShield", &fnCATSShield, "fnCATSShield/I", basket);
		fTree->Branch( "iCATSShield", fiCATSShield, "fiCATSShield[fnCATSShield]/I", basket);
		fTree->Branch( "eCATSShield", feCATSShield, "feCATSShield[fnCATSShield]/F", basket);
		fTree->Branch( "tCATSShield", ftCATSShield, "ftCATSShield[fnCATSShield]/F", basket);

		// Rear Shield
		fTree->Branch( "nCATSRShield", &fnCATSRShield, "fnCATSRShield/I", basket);
		fTree->Branch( "iCATSRShield", fiCATSRShield, "fiCATSRShield[fnCATSRShield]/I", basket);
		fTree->Branch( "eCATSRShield", feCATSRShield, "feCATSRShield[fnCATSRShield]/F", basket);
		fTree->Branch( "tCATSRShield", ftCATSRShield, "ftCATSRShield[fnCATSRShield]/F", basket);

		// Veto
		fTree->Branch( "nCATSVeto", &fnCATSVeto, "fnCATSVeto/I", basket);
		fTree->Branch( "iCATSVeto", fiCATSVeto, "fiCATSVeto[fnCATSVeto]/I", basket);
		fTree->Branch( "eCATSVeto", feCATSVeto, "feCATSVeto[fnCATSVeto]/F", basket);
		fTree->Branch( "tCATSVeto", ftCATSVeto, "ftCATSVeto[fnCATSVeto]/F", basket);

	}

	// Store Primaries
	if (fStorePrimaries)
	{
		G4cout << "Storing IDs of primary particles" << G4endl;
		fTree->Branch("pcryst",fpcryst,"fpcryst[fnhits]/I",basket);
		fTree->Branch("pctaps",fpctaps,"fpctaps[fntaps]/I",basket);
		fTree->Branch("pveto",fpveto,"fpveto[fvhits]/I",basket);
		fTree->Branch("pvtaps",fpvtaps,"fpvtaps[fnvtaps]/I",basket);
	}

	if (fIsGiBUU) fTree->Branch("weight",&fweight,"fweight/F",basket);

}

void A2CBOutput::WriteHit(G4HCofThisEvent* HitsColl){
  G4int CollSize=HitsColl->GetNumberOfCollections();
  //G4cout<<"Collection size "<<CollSize<<" "<<HitsColl->GetHC(0)->GetName()<<" "<<HitsColl->GetHC(1)->GetName()<<G4endl;
  //G4cout<<"Collection size "<<CollSize<<G4endl;
  fnhits=fntaps=fnvtaps=fvhits=fntof=fnpiz=fnmwpc=fnhe3=fntpc=0;
  fnCATSCore = fnCATSAnn = fnCATSShield = fnCATSRShield = fnCATSVeto= 0;  
  fetot=0;
  G4int hci=0;
  for(G4int i=0;i<CollSize;i++){
    //Get the ball hit info to be written to output
    A2HitsCollection* hc;
    //Look for the hits collections with entries 
    while(!(hc=static_cast<A2HitsCollection*>(HitsColl->GetHC(hci++))));
    //if(!hc)continue; //no hits in that detector
    G4int hc_nhits=hc->entries();
    //G4cout<<i<<" "<<hc->GetName()<< " "<<hc_nhits<<G4endl;
    if(hc->GetName()=="A2SDHitsCBSD"||hc->GetName()=="A2SDHitsVisCBSD"){
      fnhits=hc_nhits;
      for(Int_t ii=0;ii<fnhits;ii++){
	A2Hit* hit=static_cast<A2Hit*>(hc->GetHit(ii));
	fecryst[ii]=hit->GetEdep()/GeV;
	ftcryst[ii]=hit->GetTime()/ns;
	ficryst[ii]=hit->GetID();
	fpcryst[ii]=hit->GetParticle();
	fetot+=fecryst[ii];
      }
    }
    if(hc->GetName()=="A2SDHitsTAPSSD"||hc->GetName()=="A2SDHitsTAPSVisSD"){
      fntaps=hc_nhits;
      for(Int_t ii=0;ii<fntaps;ii++){
	A2Hit* hit=static_cast<A2Hit*>(hc->GetHit(ii));
	fectapsl[ii]=hit->GetEdep()/GeV;
	fictaps[ii]=hit->GetID();
	ftctaps[ii]=hit->GetTime()/ns;
	fpctaps[ii]=hit->GetParticle();
	//fetot+=fectapsl[i];//***TEMP!!!!
      }
    }
    if(hc->GetName()=="A2SDHitsTAPSVSD"||hc->GetName()=="A2SDHitsTAPSVVisSD"){
      fnvtaps=hc_nhits;
      for(Int_t ii=0;ii<fnvtaps;ii++){
	A2Hit* hit=static_cast<A2Hit*>(hc->GetHit(ii));
	fevtaps[ii]=hit->GetEdep()/GeV;
	fivtaps[ii]=hit->GetID();
	fpvtaps[ii]=hit->GetParticle();
      }
    }
    if(hc->GetName()=="A2SDHitsPIDSD"){
      fvhits=hc_nhits;
      for(Int_t ii=0;ii<fvhits;ii++){
	A2Hit* hit=static_cast<A2Hit*>(hc->GetHit(ii));
	feveto[ii]=hit->GetEdep()/GeV;
	ftveto[ii]=hit->GetTime()/ns;
	fiveto[ii]=hit->GetID();
	fpveto[ii]=hit->GetParticle();
      }
    }
		if(G4StrUtil::contains(hc->GetName(),"A2MWPCSD"))
		{
//    if(hc->GetName().contains("A2MWPCSD")){
      // fnmwpc=hc_nhits;
      for(Int_t ii=0;ii<hc_nhits;ii++){
	A2Hit* hit=static_cast<A2Hit*>(hc->GetHit(ii));
	fimwpc[fnmwpc+ii]   = hit->GetID();
	fmposx[fnmwpc+ii] = (Float_t)hit->GetPos().getX();
	fmposy[fnmwpc+ii] = (Float_t)hit->GetPos().getY();
	fmposz[fnmwpc+ii] = (Float_t)hit->GetPos().getZ();
	femwpc[fnmwpc+ii] = (Float_t)hit->GetEdep()/GeV;
	//	G4cout<<fnmwpc+i<<" "<<fimwpc[fnmwpc+i]<<" "<<fmpos[fnmwpc+i][0]<<" "<<fmpos[fnmwpc+i][1]<<" "<<fmpos[fnmwpc+i][2]<<G4endl;
      }
      fnmwpc+=hc_nhits;
    }
    if(hc->GetName()=="A2SDHitsTOFSD"){
      fntof=hc_nhits;
      for(Int_t ii=0;ii<fntof;ii++){
	A2Hit* hit=static_cast<A2Hit*>(hc->GetHit(ii));
	ftofe[ii]=hit->GetEdep()/GeV;
	ftoft[ii]=hit->GetTime()/ns;
	ftofx[ii]=hit->GetPos().x()/cm;
	ftofy[ii]=hit->GetPos().y()/cm;
	ftofz[ii]=hit->GetPos().z()/cm;
	ftofi[ii]=hit->GetID();
      }
    }
    if(hc->GetName()=="A2SDHitsPizzaSD" || hc->GetName()=="A2SDHitsPizzaVisSD"){
      fnpiz=hc_nhits;
      for(Int_t ii=0;ii<fnpiz;ii++){
	A2Hit* hit=static_cast<A2Hit*>(hc->GetHit(ii));
	fepiz[ii]=hit->GetEdep()/GeV;
	ftpiz[ii]=hit->GetTime()/ns;
	fipiz[ii]=hit->GetID();
      }
    }

		// CATS Core
		if ( hc->GetName() == "A2SDHitsCATSCoreSD" || hc->GetName() == "A2SDHitsCATSCoreVisSD")
		{
	      fnCATSCore = hc_nhits;
	      for ( Int_t ii = 0; ii < fnCATSCore; ii++)
			{
				A2Hit* hit=static_cast<A2Hit*>(hc->GetHit(ii));
				feCATSCore[ii]=hit->GetEdep()/GeV;
				ftCATSCore[ii]=hit->GetTime()/ns;
				fiCATSCore[ii]=hit->GetID();
//				G4cout << " fnCATSCore = " << fnCATSCore;
//				G4cout << " ii = " << ii;
//				G4cout << " hit = " << hit->GetID();
//				G4cout << G4endl;
			}
		}

		// CATS Annulus
		if ( hc->GetName() == "A2SDHitsCATSAnnSD" || hc->GetName() == "A2SDHitsCATSAnnVisSD")
		{
	      fnCATSAnn = hc_nhits;
	      for ( Int_t ii = 0; ii < fnCATSAnn; ii++)
			{
				A2Hit* hit=static_cast<A2Hit*>(hc->GetHit(ii));
				feCATSAnn[ii]=hit->GetEdep()/GeV;
				ftCATSAnn[ii]=hit->GetTime()/ns;
				fiCATSAnn[ii]=hit->GetID();
//				G4cout << " fnCATSAnn = " << fnCATSAnn;
//				G4cout << " ii = " << ii;
//				G4cout << " hit = " << hit->GetID();
//				G4cout << G4endl;
			}
		}

		// CATS Cosmic Shield
		if ( hc->GetName() == "A2SDHitsCATSShieldSD" || hc->GetName() == "A2SDHitsCATSShieldVisSD")
		{
	      fnCATSShield = hc_nhits;
	      for ( Int_t ii = 0; ii < fnCATSShield; ii++)
			{
				A2Hit* hit=static_cast<A2Hit*>(hc->GetHit(ii));
				feCATSShield[ii]=hit->GetEdep()/GeV;
				ftCATSShield[ii]=hit->GetTime()/ns;
				fiCATSShield[ii]=hit->GetID();
//				G4cout << " fnCATSShield = " << fnCATSShield;
//				G4cout << " ii = " << ii;
//				G4cout << " hit = " << hit->GetID();
//				G4cout << G4endl;
			}
		}

		// CATS Rear Cosmic Shield
		if ( hc->GetName() == "A2SDHitsCATSRShieldSD" || hc->GetName() == "A2SDHitsCATSRShieldVisSD")
		{
	      fnCATSRShield = hc_nhits;
	      for ( Int_t ii = 0; ii < fnCATSRShield; ii++)
			{
				A2Hit* hit=static_cast<A2Hit*>(hc->GetHit(ii));
				feCATSRShield[ii]=hit->GetEdep()/GeV;
				ftCATSRShield[ii]=hit->GetTime()/ns;
				fiCATSRShield[ii]=hit->GetID();
//				G4cout << " fnCATSRShield = " << fnCATSRShield;
//				G4cout << " ii = " << ii;
//				G4cout << " hit = " << hit->GetID();
//				G4cout << G4endl;
			}
		}

		// CATS Veto
		if ( hc->GetName() == "A2SDHitsCATSVetoSD" || hc->GetName() == "A2SDHitsCATSVetoVisSD")
		{
	      fnCATSVeto = hc_nhits;
	      for ( Int_t ii = 0; ii < fnCATSVeto; ii++)
			{
				A2Hit* hit=static_cast<A2Hit*>(hc->GetHit(ii));
				feCATSVeto[ii]=hit->GetEdep()/GeV;
				ftCATSVeto[ii]=hit->GetTime()/ns;
				fiCATSVeto[ii]=hit->GetID();
//				G4cout << " fnCATSVeto = " << fnCATSVeto;
//				G4cout << " ii = " << ii;
//				G4cout << " hit = " << hit->GetID();
//				G4cout << G4endl;
			}
		}


    if(hc->GetName()=="A2SDHitsAHe3SD"){
      fnhe3=hc_nhits;
      for(Int_t ii=0;ii<fnhe3;ii++){
	A2Hit* hit=static_cast<A2Hit*>(hc->GetHit(ii));
    fehe3[ii]=hit->GetEdep()/GeV;
    fthe3[ii]=hit->GetTime()/ms;
	//fhe3x[ii]=hit->GetPos().x()/cm;
	//fhe3y[ii]=hit->GetPos().y()/cm;
	//fhe3z[ii]=hit->GetPos().z()/cm;
    fihe3[ii]=hit->GetID();
      }
    }
    //new 2021: adding TPC anode hits
    if(hc->GetName()=="A2SDHitsAnodeSD"){
      fntpc=hc_nhits;
      for(Int_t ii=0;ii<fntpc;ii++){
        A2Hit* hit=static_cast<A2Hit*>(hc->GetHit(ii));
	//hit->Print();
    fqtpc[ii]=hit->GetQdep()/eplus; //units???
    fttpc[ii]=hit->GetTime()/ms;
        //ftpcx[ii]=hit->GetPos().x()/cm;
        //ftpcy[ii]=hit->GetPos().y()/cm;
        //ftpcz[ii]=hit->GetPos().z()/cm;
    fitpc[ii]=hit->GetID();
      }
    }
  }
  
}
void A2CBOutput::WriteGenInput(){
  //Note fvertex is already the pointer to fPGA::fGenPosition 
  //Get the generated input info to be written to output 
  TVector3 vec=fBeamLorentzVec->Vect().Unit();
  fbeam[0]=vec.X();
  fbeam[1]=vec.Y();
  fbeam[2]=vec.Z();
  fbeam[3]=fBeamLorentzVec->Rho()/GeV;
  fbeam[4]=fBeamLorentzVec->E()/GeV;

  //Loop over the input particles and write their real kinematics
  fnpart=fPGA->GetNGenParticles();
  for(Int_t i=0;i<fnpart;i++){
    vec=fGenLorentzVec[i]->Vect().Unit();
    fdircos[i][0]=static_cast<Float_t>(vec.X());
    fdircos[i][1]=static_cast<Float_t>(vec.Y());
    fdircos[i][2]=static_cast<Float_t>(vec.Z());
    felab[i]=fGenLorentzVec[i]->E()/GeV;
    fklab[i]=(fGenLorentzVec[i]->E()-fGenLorentzVec[i]->M())/MeV;
    fplab[i]=fGenLorentzVec[i]->Rho()/GeV;
    //ftheta[i]=fGenLorentzVec[i]->Vect().Theta();
    fidpart[i]=fGenPartType[i];
  }
  if (fIsGiBUU) fweight = fPGA->GetFileGen()->GetWeight();
}
