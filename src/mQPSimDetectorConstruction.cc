//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
//
/// \file mQPSimDetectorConstruction.cc
/// \brief Implementation of the mQPSimDetectorConstruction class

#include "mQPSimDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"

// for the ConstructSDandField method
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"

// reflection stuff
#include "G4OpBoundaryProcess.hh"

// #include "G4Scintillation.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mQPSimDetectorConstruction::mQPSimDetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mQPSimDetectorConstruction::~mQPSimDetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* mQPSimDetectorConstruction::Construct()
{

  // get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // switch on checking of volumes overlaps
  G4bool checkOverlaps = true;

  //
  // Create the World
  //-----------------

  // define the world size
  G4double world_sizeXY = 30*cm;
  G4double world_sizeZ  = 300*cm;

  // make the world a box
  G4Box* solidWorld =
    new G4Box("World",                                           //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size

  // define the world is full of air
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

  // give air a refractive index (needed for optical photons)
  G4MaterialPropertiesTable *air_mt = new G4MaterialPropertiesTable();
  G4double air_Energy[3]={2.0*eV,7.0*eV,7.14*eV};
  G4double air_RIND[3]={1.,1.,1.};
  air_mt->AddProperty("RINDEX", air_Energy, air_RIND,3);
  world_mat->SetMaterialPropertiesTable(air_mt);


  // create the world logical volume
  G4LogicalVolume* logicWorld =
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name

  // create the world physical volume
  G4VPhysicalVolume* physWorld =
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking


  //
  // Create the scintillator
  // -----------------------

  // define scintillator size & shape
  G4double scintillator_dx =  10*cm;
  G4double scintillator_dy =  10*cm;
  G4double scintillator_dz = 120*cm;

  G4Box* solidScintillator =
    new G4Box("ScintillatorBox",          //its name
              0.5*scintillator_dx,
              0.5*scintillator_dy,
              0.5*scintillator_dz);       //its size

  // define scintillator material
  G4Material* scintillator_mat = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

  // specify scintillator proporties (this is default from Listing 5.5 on page 219 of Geant4 manual)
  const G4int NUMENTRIES = 9;
  G4double Scnt_PP[NUMENTRIES] = { 6.6*eV, 6.7*eV, 6.8*eV, 6.9*eV,
                                   7.0*eV, 7.1*eV, 7.2*eV, 7.3*eV, 7.4*eV };
  G4double Scnt_FAST[NUMENTRIES] = { 0.000134, 0.004432, 0.053991, 0.241971,
                                     0.398942, 0.000134, 0.004432, 0.053991,
                                     0.241971 };
  G4double Scnt_SLOW[NUMENTRIES] = { 0.000010, 0.000020, 0.000030, 0.004000,
                                     0.008000, 0.005000, 0.020000, 0.001000,
                                     0.000010 };

  G4MaterialPropertiesTable* Scnt_MPT = new G4MaterialPropertiesTable();
  Scnt_MPT->AddProperty("FASTCOMPONENT", Scnt_PP, Scnt_FAST, NUMENTRIES);
  Scnt_MPT->AddProperty("SLOWCOMPONENT", Scnt_PP, Scnt_SLOW, NUMENTRIES);
  //Scnt_MPT->AddConstProperty("SCINTILLATIONYIELD", 5000./MeV);
  Scnt_MPT->AddConstProperty("SCINTILLATIONYIELD", 0.1/MeV); //reduce yield for testing
  Scnt_MPT->AddConstProperty("RESOLUTIONSCALE", 2.0);
  Scnt_MPT->AddConstProperty("FASTTIMECONSTANT",  1.*ns);
  Scnt_MPT->AddConstProperty("SLOWTIMECONSTANT", 10.*ns);
  Scnt_MPT->AddConstProperty("YIELDRATIO", 0.8);
  scintillator_mat->SetMaterialPropertiesTable(Scnt_MPT);

  // set Birk's constant
  scintillator_mat->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

  // create scintillator logical volume
  G4LogicalVolume* logicScintillator =
    new G4LogicalVolume(solidScintillator,         //its solid
                        scintillator_mat,          //its material
                        "ScintillatorLV");         //its name

  // define scintillator position
  G4ThreeVector scintillator_pos = G4ThreeVector(0, 0, 0);

  // create scintillator physical volume and place it in the world volume
  G4VPhysicalVolume* physScintillator =
                    new G4PVPlacement(0,     //no rotation
                    scintillator_pos,        //at position
                    logicScintillator,       //its logical volume
                    "ScintillatorPV",        //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  // define scintillator visualization attributes
  G4VisAttributes* visAttributes = new G4VisAttributes(G4Colour(0.0,0.3,0.9,0.5));
  visAttributes->SetVisibility(true);
  logicScintillator->SetVisAttributes(visAttributes);


  // following https://www-zeuthen.desy.de/geant4/g4course2011/day2/6_opticalphysics/index.html
  // and Listing 5.8 of the manual, make the scintillator surface reflective

  // define the optical surface to be the interesection between the scintillator and the world
  G4OpticalSurface *scintWrap = new G4OpticalSurface("ScintWrap");
  new G4LogicalBorderSurface("ScintWrap",        // name
                            physScintillator,    // from volume
                            physWorld,           // to volume
                            scintWrap);          // optical surface

  // set the optical boundary model (see manual page 222 'Boundary Process')
  scintWrap->SetModel(unified);
  scintWrap->SetType(dielectric_dielectric);
  scintWrap->SetFinish(polishedfrontpainted);

  // define the optical boundary properties
  const G4int NUM = 2;
  G4double pp[NUM] = {2.0*eV, 4.0*eV};
  G4double rindex[NUM] = {1.35, 1.40};
  G4double reflectivity[NUM] = {1.0, 1.0};
  G4double efficiency[NUM] = {0.0, 0.0};

  // set the optical boundary properties
  G4MaterialPropertiesTable* scintWrapProperty = new G4MaterialPropertiesTable();
  scintWrapProperty->AddProperty("RINDEX",pp,rindex,NUM);
  scintWrapProperty->AddProperty("REFLECTIVITY",pp,reflectivity,NUM);
  scintWrapProperty->AddProperty("EFFICIENCY",pp,efficiency,NUM);
  scintWrap->SetMaterialPropertiesTable(scintWrapProperty);


  // set scintillator as scoring volume
  fScoringVolume = logicScintillator;

  // set the scintillator as a sensitive detector
  //fSensitiveDetector->SetSensitiveDetector(logicScintillator)


/* LEAVE THE LIGHT GUIDE OUT FOR NOW...
  //
  // Create the light guide
  // ----------------------

  // light guide material and position
  G4Material* lightguide_mat = nist->FindOrBuildMaterial("G4_PLEXIGLASS");
  G4ThreeVector lightguide_pos = G4ThreeVector(0, 0, -80*cm);

  // light guide conical section shape
  G4double lightguide_rmina =  0.*cm, lightguide_rmaxa = 2.*cm;
  G4double lightguide_rminb =  0.*cm, lightguide_rmaxb = 5.*cm;
  G4double lightguide_hz = 20.*cm;
  G4double lightguide_phimin = 0.*deg, lightguide_phimax = 360.*deg;
  G4Cons* solidLightguide =
    new G4Cons("Lightguide",
    lightguide_rmina, lightguide_rmaxa, lightguide_rminb, lightguide_rmaxb, lightguide_hz,
    lightguide_phimin, lightguide_phimax);

  //light guide logical volume
  G4LogicalVolume* logicLightguide =
    new G4LogicalVolume(solidLightguide,         //its solid
                        lightguide_mat,          //its material
                        "Lightguide");           //its name
  // light guide placement
  new G4PVPlacement(0,                       //no rotation
                    lightguide_pos,          //at position
                    logicLightguide,         //its logical volume
                    "Lightguide",            //its name
                    //logicEnv,                //its mother  volume
                    logicWorld,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  //
  // lightguide visualization attributes
  logicLightguide->SetVisAttributes(visAttributes);
*/


  // always return the physical World
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void mQPSimDetectorConstruction::ConstructSDandField()
{

  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  // declare Absorber as a MultiFunctionalDetector scorer
  //
  auto scintillatorDetector = new G4MultiFunctionalDetector("ScintillatorMFD");
  G4SDManager::GetSDMpointer()->AddNewDetector(scintillatorDetector);

  G4VPrimitiveScorer* primitive;
  primitive = new G4PSEnergyDeposit("Edep");
  scintillatorDetector->RegisterPrimitive(primitive);

  // primitive = new G4PSTrackLength("TrackLength");
  // auto charged = new G4SDChargedFilter("chargedFilter");
  // primitive ->SetFilter(charged);
  // absDetector->RegisterPrimitive(primitive);

  SetSensitiveDetector("ScintillatorLV",scintillatorDetector);

  //
  // // sensitive detectors
  // auto sdManager = G4SDManager::GetSDMpointer();
  // G4String SDname;

  // auto scintillator = new mQPSimScintillatorSD(SDname="/scintillator");
  // sdManager->AddNewDetector(scintillator);
  // fScintillatorLogical->SetSensitiveDetector(scintillator);


}
