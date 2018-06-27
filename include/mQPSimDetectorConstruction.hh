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
/// \file mQPSimDetectorConstruction.hh
/// \brief Definition of the mQPSimDetectorConstruction class

#ifndef mQPSimDetectorConstruction_h
#define mQPSimDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include "G4Material.hh"
#include "G4Element.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;

/// Detector construction class to define materials and geometry.

class mQPSimDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    mQPSimDetectorConstruction();
    virtual ~mQPSimDetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }


  protected:
    G4LogicalVolume*  fScoringVolume;
    G4bool            fCheckOverlaps;

  // private:
  //   void DefineMaterials();
  //
  //   //Materials & Elements
  //   G4Material* fLXe;
  //   G4Material* fAl;
  //   G4Element* fN;
  //   G4Element* fO;
  //   G4Material* fAir;
  //   G4Material* fVacuum;
  //   G4Element* fC;
  //   G4Element* fH;
  //   G4Material* fGlass;
  //   G4Material* fPstyrene;
  //   G4Material* fPMMA;
  //   G4Material* fPethylene1;
  //   G4Material* fPethylene2;
  //
  //   G4MaterialPropertiesTable* fLXe_mt;
  //   G4MaterialPropertiesTable* fMPTPStyrene;


};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
