//
// CATS Simulation
//
// Vincent Bruening
// Mount Allison University

#include "A2DetCATS.hh"
#include "CLHEP/Units/SystemOfUnits.h"
#include "G4PVPlacement.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4SDManager.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4MultiUnion.hh"
#include "G4RotationMatrix.hh"
#include "G4SubtractionSolid.hh"
#include "G4Cons.hh"
#include "A2Utils.hh"

#include "TString.h"

using namespace CLHEP;

// half-length of CATS
const G4double z_cryst = 31.75*cm; // For crystal
const G4double w_scint = 40.75*cm; // For scintillators

const G4int nAnn = 6;	// Number of Annulus segments
const G4int nShi = 5;	// Number of Cosmic Shield segments

// Constructor
A2DetCATS::A2DetCATS()
{ 

	fIsCheckOverlap = true;

	fregionCATS = new G4Region( "CATS");

	fCATSCoreSD = 0;
	fCATSCoreVisSD = 0;
	fCATSAnnSD = 0;
	fCATSAnnVisSD = 0;
	fCATSShieldSD = 0;
	fCATSShieldVisSD = 0;
	fCATSRShieldSD = 0;
	fCATSRShieldVisSD = 0;
	fCATSVetoSD = 0;
	fCATSVetoVisSD = 0;

	fNistManager = G4NistManager::Instance();

}

// Destructor
A2DetCATS::~A2DetCATS()
{

	if (fregionCATS) delete fregionCATS;
	if (fCATSCoreSD) delete fCATSCoreSD;
	if (fCATSAnnSD) delete fCATSAnnSD;
	if (fCATSShieldSD) delete fCATSShieldSD;
	if (fCATSRShieldSD) delete fCATSRShieldSD;
	if (fCATSVetoSD) delete fCATSVetoSD;
	if (fCATSCoreVisSD) delete fCATSCoreVisSD;
	if (fCATSAnnVisSD) delete fCATSAnnVisSD;
	if (fCATSShieldVisSD) delete fCATSShieldVisSD;
	if (fCATSRShieldVisSD) delete fCATSRShieldVisSD;
	if (fCATSVetoVisSD) delete fCATSVetoVisSD;

}

//Main Construction instructions
G4VPhysicalVolume* A2DetCATS::Construct( G4LogicalVolume* motherLogic)
{
	// Set the mother volume
	fMotherLogic = motherLogic;

	// Construct the Air Box
	G4Box* airbox = new G4Box( "litterbox", 70*cm, 70*cm, 100*cm);
	fMyLogic = new G4LogicalVolume( airbox, fNistManager->FindOrBuildMaterial( "G4_AIR"), "litterbox");
	fMyLogic->SetVisAttributes( G4VisAttributes::GetInvisible); 

	G4RotationMatrix* rot = new G4RotationMatrix();
	rot->rotateY( 40.*deg);

//	fMyPhysi = new G4PVPlacement( rot, G4ThreeVector( -96*cm, 0, 315*cm), fMyLogic, "A2DetCATS", fMotherLogic, false, 0);

//	fMyPhysi = new G4PVPlacement( rot, G4ThreeVector( 0*cm, 0*cm, 150*cm), fMyLogic, "litterbox", fMotherLogic, false, 0);  
	fMyPhysi = new G4PVPlacement( rot, G4ThreeVector( -96*cm, 0*cm, 115*cm), fMyLogic, "litterbox", fMotherLogic, false, 0);  

	if ( fIsCheckOverlap) A2Utils::CheckOverlapAndAbort( fMyPhysi, "A2DetCATS::Construct()");

// Construct the various parts of CATS

	MakeCore();
	MakeAnnulus();
	MakeShield();
	MakeRearShield();
	MakeVeto();
	MakeRing();
	MakeLeadShield();

	return fMyPhysi;

}

// Function to build NaI CATS core
void A2DetCATS::MakeCore()
{

// Define Core colour
	G4VisAttributes* col1 = new G4VisAttributes( G4Colour( 0.4, 0.5, 0.1));
	col1->SetVisibility( true);

// Make the Core
	G4Tubs *fCore = new G4Tubs( "Core", 0*cm, 13.35*cm, z_cryst, 0.*deg, 360.*deg);
	G4LogicalVolume* fCoreLogic = new G4LogicalVolume( fCore, fNistManager->FindOrBuildMaterial( "G4_SODIUM_IODIDE"), "CoreLogic");
	fCoreLogic->SetVisAttributes( col1);
	G4VPhysicalVolume* fCorePhysi = new G4PVPlacement( 0, G4ThreeVector( 0, 0, 0), fCoreLogic, "CorePlacement", fMyLogic, false, 1);

// Check for overlap
	if ( fIsCheckOverlap) A2Utils::CheckOverlapAndAbort( fCorePhysi, "A2DetCATS::Construct()");

// Make Sensitive
	G4SDManager* SDman = G4SDManager::GetSDMpointer();
	if ( fIsInteractive == 1)
	{
		if ( !fCATSCoreVisSD) fCATSCoreVisSD = new A2VisSD( "CATSCoreVisSD", 1);
		SDman->AddNewDetector( fCATSCoreVisSD );
		fCoreLogic->SetSensitiveDetector( fCATSCoreVisSD);
		fregionCATS->AddRootLogicalVolume( fCoreLogic);
	}
	else
	{
		// Core
		if ( !fCATSCoreVisSD) fCATSCoreVisSD = new A2VisSD( "CATSCoreVisSD", 1);
		SDman->AddNewDetector( fCATSCoreVisSD );
		fCoreLogic->SetSensitiveDetector( fCATSCoreVisSD);
		fregionCATS->AddRootLogicalVolume( fCoreLogic);
	}

}

//
// Annulus
//
void A2DetCATS::MakeAnnulus()
{

//	Define Annulus colours
	G4VisAttributes* col2 = new G4VisAttributes( G4Colour( 0.0, 1.0, 1.0));
	col2->SetVisibility( true);
	G4VisAttributes* col3 = new G4VisAttributes( G4Colour( 1.0, 0.0, 0.0));
	col3->SetVisibility( true);

	G4Tubs *fAnnulus = new G4Tubs( "Annulus", 13.35*cm, 24.15*cm, z_cryst, 30.*deg, 60.*deg);
	G4LogicalVolume* fAnnulusLogic = new G4LogicalVolume( fAnnulus, fNistManager->FindOrBuildMaterial( "G4_SODIUM_IODIDE"), "AnnulusLogic");

// Make the 6 Annulus segments
	for ( G4int ii = 0; ii < 6; ii++)
	{

		G4RotationMatrix* rotate = new G4RotationMatrix();
		G4double rot_ang;

		rot_ang = 30.*deg*(1. + 2.*ii);

		rotate->rotateZ( rot_ang);

		fAnnulusLogic->SetVisAttributes( col2);
		G4VPhysicalVolume* fAnnulusPhysi = new G4PVPlacement( rotate, G4ThreeVector( 0, 0, 0), fAnnulusLogic,
				TString::Format( "AnnulusPlacement%d", ii+1).Data(), fMyLogic, false, ii+1);

		// Check for overlap
		if ( fIsCheckOverlap) A2Utils::CheckOverlapAndAbort( fAnnulusPhysi, "A2DetCATS::Construct()");

	}

// Make Sensitive
	G4SDManager* SDman = G4SDManager::GetSDMpointer();

	if ( fIsInteractive == 1)
	{ 

		if ( !fCATSAnnVisSD) fCATSAnnVisSD = new A2VisSD( "CATSAnnVisSD", nAnn);
		SDman->AddNewDetector( fCATSAnnVisSD );
		fAnnulusLogic->SetSensitiveDetector( fCATSAnnVisSD);
		fregionCATS->AddRootLogicalVolume( fAnnulusLogic);

	}
	else
	{

		if ( !fCATSAnnSD) fCATSAnnSD = new A2SD( "CATSAnnSD", nAnn);
		SDman->AddNewDetector( fCATSAnnSD );
		fAnnulusLogic->SetSensitiveDetector( fCATSAnnSD);
		fregionCATS->AddRootLogicalVolume( fAnnulusLogic);

	}
}

//
// Cosmic Veto Shield
//
void A2DetCATS::MakeShield()
{

	// Define Shield colour
	G4VisAttributes* col4 = new G4VisAttributes( G4Colour( 1.0, 1.0, 0.0));
	col4->SetVisibility( true);

	G4Tubs *fScint = new G4Tubs( "Scint", 25.15*cm, 35.15*cm, w_scint, 36.*deg, 72*deg);
	G4LogicalVolume* fScintLogic = new G4LogicalVolume( fScint, fNistManager->FindOrBuildMaterial( "G4_PLASTIC_SC_VINYLTOLUENE"), "ScintLogic");
	fScintLogic->SetVisAttributes( col4);

	// Make the 5 segments around the NaI
	for ( G4int ii = 0; ii < nShi; ii++)
	{

		G4RotationMatrix* rotate = new G4RotationMatrix();
		G4double rot_ang;

		rot_ang = 36.*deg*(1.+2.*ii);

		rotate->rotateZ( rot_ang);

		G4VPhysicalVolume* fScintPhysi = new G4PVPlacement( rotate, G4ThreeVector( 0, 0, 0), fScintLogic,
				TString::Format( "ScintPlacement%d", ii+1).Data(), fMyLogic, false, ii+1);

		if ( fIsCheckOverlap) A2Utils::CheckOverlapAndAbort( fScintPhysi, "A2DetCATS::Construct()");

	}

// Make Sensitive
	G4SDManager* SDman = G4SDManager::GetSDMpointer();
	if ( fIsInteractive == 1)
	{ 
		if ( !fCATSShieldVisSD) fCATSShieldVisSD = new A2VisSD( "CATSShieldVisSD", nShi);
		SDman->AddNewDetector( fCATSShieldVisSD );
		fScintLogic->SetSensitiveDetector( fCATSShieldVisSD);
		fregionCATS->AddRootLogicalVolume( fScintLogic);
	}
	else
	{
		if ( !fCATSShieldSD) fCATSShieldSD = new A2SD( "CATSShieldSD", nShi);
		SDman->AddNewDetector( fCATSShieldSD );
		fScintLogic->SetSensitiveDetector( fCATSShieldSD);
		fregionCATS->AddRootLogicalVolume( fScintLogic);
	}
}

//
// Rear Cosmic Shield
//
void A2DetCATS::MakeRearShield()
{

	// Define Shield colour
	G4VisAttributes* col4 = new G4VisAttributes( G4Colour( 1.0, 1.0, 0.0));
	col4->SetVisibility( true);

// Make Rear Shield 
	G4Tubs *fScintBack = new G4Tubs( "ScintBack", 0*cm, 35.0*cm, 4.7*cm, 0.*deg, 360*deg);
	G4LogicalVolume* fScintBackLogic = new G4LogicalVolume( fScintBack, fNistManager->FindOrBuildMaterial( "G4_PLASTIC_SC_VINYLTOLUENE"), "ScintBackLogic");
	fScintBackLogic->SetVisAttributes( col4);
	G4VPhysicalVolume* fScintBackPhysi = new G4PVPlacement( 0, G4ThreeVector(0,0,55.45*cm), fScintBackLogic, "BackScintPlacement", fMyLogic, false, 1);

	if ( fIsCheckOverlap) A2Utils::CheckOverlapAndAbort( fScintBackPhysi, "A2DetCATS::Construct()");

// Make Sensitive
	G4SDManager* SDman = G4SDManager::GetSDMpointer();
	if ( fIsInteractive == 1)
	{ 
		if ( !fCATSRShieldVisSD) fCATSRShieldVisSD = new A2VisSD( "CATSRShieldVisSD", 1);
		SDman->AddNewDetector( fCATSRShieldVisSD );
		fScintBackLogic->SetSensitiveDetector( fCATSRShieldVisSD);
		fregionCATS->AddRootLogicalVolume( fScintBackLogic);
	}
	else
	{
		if ( !fCATSRShieldSD) fCATSRShieldSD = new A2SD( "CATSRShieldSD", 1);
		SDman->AddNewDetector( fCATSRShieldSD );
		fScintBackLogic->SetSensitiveDetector( fCATSRShieldSD);
		fregionCATS->AddRootLogicalVolume( fScintBackLogic);
	}
}

//
// Front Veto covering entrance
//
void A2DetCATS::MakeVeto()
{

	G4VisAttributes* col7 = new G4VisAttributes( G4Colour(0,0.3,1.0));
	col7->SetVisibility( true);

	G4Box* fVeto = new G4Box( "Veto", 7*cm, 7*cm, 2.5*mm);
	G4LogicalVolume* fVetoLogic = new G4LogicalVolume(fVeto, fNistManager->FindOrBuildMaterial( "G4_PLASTIC_SC_VINYLTOLUENE"), "VetoLogic");
	fVetoLogic->SetVisAttributes( col7);
	G4VPhysicalVolume* fVetoPhysi = new G4PVPlacement( 0, G4ThreeVector(0,5*mm,-805.5*mm), fVetoLogic, "VetoPlacement", fMyLogic, false, 1);

	if ( fIsCheckOverlap) A2Utils::CheckOverlapAndAbort( fVetoPhysi, "A2DetCATS::Construct()");

	// Make Sensitive
	G4SDManager* SDman = G4SDManager::GetSDMpointer();
	if ( fIsInteractive==1)
	{ 
		if ( !fCATSVetoVisSD) fCATSVetoVisSD = new A2VisSD( "CATSVetoVisSD", 1);
		SDman->AddNewDetector( fCATSVetoVisSD);
		fVetoLogic->SetSensitiveDetector( fCATSVetoVisSD);
		fregionCATS->AddRootLogicalVolume( fVetoLogic);
	}
	else
	{
		if ( !fCATSVetoSD) fCATSVetoSD = new A2SD( "CATSVetoSD", 1);
		SDman->AddNewDetector( fCATSVetoSD );
		fVetoLogic->SetSensitiveDetector( fCATSVetoSD);
		fregionCATS->AddRootLogicalVolume( fVetoLogic);
	}
}

//
// Lithium-Carbonate Ring around NaI
//
void A2DetCATS::MakeRing()
{

	G4VisAttributes* col4 = new G4VisAttributes( G4Colour( 1.0, 0.0, 1.0));
	col4->SetVisibility( true);

	G4Tubs *fRing = new G4Tubs( "Ring", 24.15*cm, 25.15*cm, z_cryst, 0.*deg, 360.*deg);
	G4LogicalVolume* fRingLogic = new G4LogicalVolume( fRing, fNistManager->FindOrBuildMaterial( "G4_LITHIUM_CARBONATE"), "RingLogic");
	fRingLogic->SetVisAttributes( col4);
	G4VPhysicalVolume* fRingPhysi = new G4PVPlacement( 0, G4ThreeVector( 0, 0, 0), fRingLogic, "RingPlacement", fMyLogic, false, 1);

	if ( fIsCheckOverlap) A2Utils::CheckOverlapAndAbort( fRingPhysi, "A2DetCATS::Construct()");

}

//
// Lead Shielding
//
void A2DetCATS::MakeLeadShield()
{

	G4double Leadx = 351.5*mm; 

	G4double Longy = 50*mm;
	G4double Longzz = 650.5*mm;

//Two of the "Long" parts of the box
	G4Box* fLongBit1 = new G4Box("LongBit1", Leadx+111*mm, Longy, Longzz);
	G4Box* fLongBit2 = new G4Box("LongBit2", Leadx+111*mm, Longy, Longzz); 

//The other two "Long" parts
	G4Box* fOtherLongBit1 = new G4Box("OtherLongBit1", Longy, Leadx+10*mm, Longzz);
	G4Box* fOtherLongBit2 = new G4Box("OtherLongBit2", Longy, Leadx+10*mm, Longzz);

	G4Box* fEndBit = new G4Box("EndBit", Leadx+110*mm, Leadx+105*mm, Longy);

//subtraction solid here!
	G4Box* fFrontBit = new G4Box("FrontBit", 460.5*mm, 460*mm, 105*mm);
	G4Tubs* fFrontBitHole = new G4Tubs("FrontBitHole", 0*cm, 6.9*cm, 105.1*mm, 0*deg, 360*deg);
	//had to make the hole slightly longer than the box- otherwise you could not see the hole on both sides of the front bit

	G4SubtractionSolid *subtraction = new G4SubtractionSolid("subtraction", fFrontBit, fFrontBitHole);

//making the inner cone! 
	G4Cons* fLeadCone = new G4Cons("LeadCone", 138*mm, 260*mm, 187*mm, 260*mm, 73*mm, 0*deg, 360*deg);
	G4LogicalVolume* fLeadConeLogic = new G4LogicalVolume(fLeadCone, fNistManager->FindOrBuildMaterial( "G4_Pb"), "LeadConeLogic");
	G4VPhysicalVolume* fLeadConePhysi = new G4PVPlacement( 0, G4ThreeVector(0,0,-480.5*mm), fLeadConeLogic, "PLeadCone", fMyLogic, false, 1);

	if ( fIsCheckOverlap) A2Utils::CheckOverlapAndAbort( fLeadConePhysi, "A2DetCATS::Construct()");

	G4VisAttributes* col5 = new G4VisAttributes( G4Colour(0,0,1.0));
	col5->SetVisibility(true);
//col5->SetVisibility(false);
	fLeadConeLogic->SetVisAttributes(col5);

//endbit transform
	G4ThreeVector EndDisplacement = G4ThreeVector(0,0,66.15*cm);
	G4RotationMatrix none = G4RotationMatrix(); //for no rotation
	G4Transform3D EndBitTransform = G4Transform3D(none, EndDisplacement);

	G4double Placement = 401.5*mm;
	G4double LongBitZ= -1.3*cm;

//LongBit1 transform
	G4ThreeVector HighLongBit1 = G4ThreeVector(0,Placement+1*cm,LongBitZ); 
	G4Transform3D LongBit1Transform = G4Transform3D(none, HighLongBit1);

//LongBit2 transform
	G4ThreeVector LowLongBit2 = G4ThreeVector(0,-Placement,LongBitZ);
	G4Transform3D LongBit2Transform = G4Transform3D(none, LowLongBit2);

//OtherLongBit1 transform
	G4ThreeVector OLB1Displacement = G4ThreeVector(Placement+1*cm,0,LongBitZ);
	G4Transform3D OLB1Transform = G4Transform3D(none, OLB1Displacement);

//OtherLongBit2 transform
	G4ThreeVector OLB2Displacement = G4ThreeVector(-Placement-1*cm,0,LongBitZ);
	G4Transform3D OLB2Transform = G4Transform3D(none, OLB2Displacement);
//front
	G4ThreeVector FrontDisplacement = G4ThreeVector(1*mm,5*mm,-695.5*mm);//1,0,-495.5
	G4Transform3D FrontBitTransform = G4Transform3D(none, FrontDisplacement);

	G4MultiUnion* LeadBox = new G4MultiUnion("LeadBox");

	LeadBox->AddNode(*fEndBit, EndBitTransform);
	LeadBox->AddNode(*fLongBit1, LongBit1Transform);
	LeadBox->AddNode(*fLongBit2, LongBit2Transform);
	LeadBox->AddNode(*fOtherLongBit1, OLB1Transform);
	LeadBox->AddNode(*fOtherLongBit2, OLB2Transform);
	LeadBox->AddNode(*subtraction, FrontBitTransform);

	LeadBox->Voxelize();

	G4LogicalVolume* fLeadBoxLogic = new G4LogicalVolume(LeadBox, fNistManager->FindOrBuildMaterial( "G4_Pb"), "LeadLogic");

	G4VisAttributes* col6 = new G4VisAttributes( G4Colour(0.3,1.0,0.0));
	col6->SetVisibility(true);
//col6->SetVisibility(false);

	fLeadBoxLogic->SetVisAttributes(col6);

	G4VPhysicalVolume* fLeadBoxPhysi = new G4PVPlacement( 0, G4ThreeVector( 0, 0, 0), fLeadBoxLogic, "PLeadBox", fMyLogic, false, 1);

	if ( fIsCheckOverlap) A2Utils::CheckOverlapAndAbort( fLeadBoxPhysi, "A2DetCATS::Construct()");

}
