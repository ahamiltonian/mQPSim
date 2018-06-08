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
/// \file mQPSimEventAction.cc
/// \brief Implementation of the mQPSimEventAction class

#include "mQPSimEventAction.hh"
#include "mQPSimRunAction.hh"
#include "mQPSimAnalysis.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mQPSimEventAction::mQPSimEventAction(mQPSimRunAction* runAction)
: G4UserEventAction(),
  fRunAction(runAction),
  fEdep(0.),
  fScintEdepHCID(-1)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mQPSimEventAction::~mQPSimEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4THitsMap<G4double>* mQPSimEventAction::GetHitsCollection(G4int hcID,
                                                            const G4Event* event) const
{
  auto hitsCollection
      = static_cast<G4THitsMap<G4double>*>(
          event->GetHCofThisEvent()->GetHC(hcID));

  if ( ! hitsCollection ) {
      G4ExceptionDescription msg;
      msg << "Cannot access hitsCollection ID " << hcID;
      G4Exception("mQPSimEventAction::GetHitsCollection()",
        "MyCode000X", FatalException, msg);
  }

  return hitsCollection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double mQPSimEventAction::GetSum(G4THitsMap<G4double>* hitsMap) const
{
  G4double sumValue = 0.;
  for ( auto it : *hitsMap->GetMap() ) {
    // hitsMap->GetMap() returns the map of std::map<G4int, G4double*>
    sumValue += *(it.second);
  }
  return sumValue;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void mQPSimEventAction::BeginOfEventAction(const G4Event*)
{
  fEdep = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void mQPSimEventAction::EndOfEventAction(const G4Event* event)
{


  // accumulate statistics in run action (this is for the accumulator business)
  fRunAction->AddEdep(fEdep);

  // below is for the hits collection stuff...

  // Get hist collections IDs
 if ( fScintEdepHCID == -1 ) {
   fScintEdepHCID
     = G4SDManager::GetSDMpointer()->GetCollectionID("ScintillatorMFD/Edep");
 }

 // Get sum values from hits collections
 //
 auto scintEdep = GetSum(GetHitsCollection(fScintEdepHCID, event));


  // get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // fill histograms
  //
  analysisManager->FillH1(0, fEdep);
  analysisManager->FillH1(1, scintEdep);
//  analysisManager->FillH1(2, fEdep);
//  analysisManager->FillH1(3, fEdep);
  // analysisManager->FillH1(0, absoEdep);
  // analysisManager->FillH1(1, gapEdep);
  // analysisManager->FillH1(2, absoTrackLength);
  // analysisManager->FillH1(3, gapTrackLength);

  // fill ntuple
  //
  analysisManager->FillNtupleDColumn(0, fEdep);
  analysisManager->FillNtupleDColumn(1, scintEdep);
  //analysisManager->FillNtupleDColumn(2, fEdep);
  //analysisManager->FillNtupleDColumn(3, fEdep);
  // analysisManager->FillNtupleDColumn(0, absoEdep);
  // analysisManager->FillNtupleDColumn(1, gapEdep);
  // analysisManager->FillNtupleDColumn(2, absoTrackLength);
  // analysisManager->FillNtupleDColumn(3, gapTrackLength);
  analysisManager->AddNtupleRow();

  //print per event (modulo n)
  //
  auto eventID = event->GetEventID();
  //auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  //if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {
    G4cout << "---> End of event output: " << eventID << G4endl;
    G4cout
       << "accumulator energy: "
       << std::setw(7) << G4BestUnit(fEdep, "Energy")
       << "       hits collection energy: "
       << std::setw(7) << G4BestUnit(scintEdep, "Energy")
       << G4endl;
  //}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
