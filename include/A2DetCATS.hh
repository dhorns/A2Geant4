#ifndef A2DetCATS_h
#define A2DetCATS_h 1

#include "A2Detector.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4PVPlacement.hh"
#include "G4Box.hh"
#include "A2SD.hh"
#include "A2VisSD.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Region;
class A2SD;
class A2VisSD;

class A2DetCATS : public A2Detector 
{
	protected:
		G4bool fIsCheckOverlap;

		G4Region* fregionCATS;

		A2VisSD* fCATSCoreVisSD;
		A2VisSD* fCATSAnnVisSD;
		A2VisSD* fCATSShieldVisSD;
		A2VisSD* fCATSRShieldVisSD;
		A2VisSD* fCATSVetoVisSD;

		A2SD* fCATSCoreSD;
		A2SD* fCATSAnnSD;
		A2SD* fCATSShieldSD;
		A2SD* fCATSRShieldSD;
		A2SD* fCATSVetoSD;

	public: 
		A2DetCATS();
		~A2DetCATS();

		virtual G4VPhysicalVolume* Construct( G4LogicalVolume* motherLogic); 

		void MakeCore();
		void MakeAnnulus();
		void MakeShield();
		void MakeRearShield();
		void MakeVeto();

		void MakeSensitiveDetectors();

		void MakeRing();
		void MakeLeadShield();

};

#endif
