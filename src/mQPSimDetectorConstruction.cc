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

// for the ConstructSDandField method
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"


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

  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // Envelope parameters
  //
  G4double env_sizeXY = 20*cm, env_sizeZ = 200*cm;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_AIR");

  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //
  // World
  //
  G4double world_sizeXY = 1.2*env_sizeXY;
  G4double world_sizeZ  = 1.2*env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

  G4Box* solidWorld =
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size

  G4LogicalVolume* logicWorld =
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name

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
  // Envelope
  //
  G4Box* solidEnv =
    new G4Box("Envelope",                    //its name
        0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ); //its size

  G4LogicalVolume* logicEnv =
    new G4LogicalVolume(solidEnv,            //its solid
                        env_mat,             //its material
                        "Envelope");         //its name

  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicEnv,                //its logical volume
                    "Envelope",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking


  //
  // Create the scintillator
  //

  // scintillator material and position
  G4Material* scintillator_mat = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  G4ThreeVector scintillator_pos = G4ThreeVector(0, 0, 0);

  // scintillator shape
  G4double scintillator_dx =  10*cm;
  G4double scintillator_dy =  10*cm;
  G4double scintillator_dz = 120*cm;
  G4Box* solidScintillator =
    new G4Box("ScintillatorBox",          //its name
              0.5*scintillator_dx,
              0.5*scintillator_dy,
              0.5*scintillator_dz);    //its size

  // scintillator logical volume
  G4LogicalVolume* logicScintillator =
    new G4LogicalVolume(solidScintillator,         //its solid
                        scintillator_mat,          //its material
                        "ScintillatorLV");           //its name

  // scintillator placement
  new G4PVPlacement(0,                       //no rotation
                    scintillator_pos,        //at position
                    logicScintillator,       //its logical volume
                    "ScintillatorPV",          //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  // scintillator visualization attributes
  // visAttributes = new G4VisAttributes(G4Colour(0.0,0.0,0.9));
  // visAttributes->SetVisibility(false);
  // logicScintillator->SetVisAttributes(visAttributes);
  // fVisAttributes.push_back(visAttributes);

  // set scintillator as scoring volume
  fScoringVolume = logicScintillator;

  // set the scintillator as a sensitive detector
  //fSensitiveDetector->SetSensitiveDetector(logicScintillator)

  //
  // Create the light guide
  //
  // light guide material and position
  G4Material* lightguide_mat = nist->FindOrBuildMaterial("G4_PLEXIGLASS");
  G4ThreeVector lightguide_pos = G4ThreeVector(0, 0, -70*cm);

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
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking


  //
  //always return the physical World
  //
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
