#!/usr/bin/env python
import sys
import ROOT
from collections import OrderedDict
from ROOT import TLorentzVector
from array import array
import numpy as np
import argparse


class TreeProducer:
    def __init__(self, debug):

         # flat tree branches
         self.debug = debug

         self.t = ROOT.TTree( "mytree","TestTree" )
         self.maxn = 9999

         # declare arrays
         self.evt_size = array( 'i', [ 0 ] )

         self.vtx_size         = array( 'i', [ 0 ] )
         self.vtx_pt2          = array( 'f', self.maxn*[ 0. ] )

         self.genpart_size     = array( 'i', [ 0 ] )
         self.genpart_pid      = array( 'i', self.maxn*[ 0 ] )
         self.genpart_status   = array( 'i', self.maxn*[ 0 ] )
         self.genpart_pt       = array( 'f', self.maxn*[ 0. ] )
         self.genpart_eta      = array( 'f', self.maxn*[ 0. ] )
         self.genpart_phi      = array( 'f', self.maxn*[ 0. ] )
         self.genpart_mass     = array( 'f', self.maxn*[ 0. ] )
         self.genpart_m1       = array( 'i', self.maxn*[ 0 ] )
         self.genpart_m2       = array( 'i', self.maxn*[ 0 ] )
         self.genpart_d1       = array( 'i', self.maxn*[ 0 ] )
         self.genpart_d2       = array( 'i', self.maxn*[ 0 ] )

         self.genjet_size      = array( 'i', [ 0 ] )
         self.genjet_pt        = array( 'f', self.maxn*[ 0. ] )
         self.genjet_eta       = array( 'f', self.maxn*[ 0. ] )
         self.genjet_phi       = array( 'f', self.maxn*[ 0. ] )
         self.genjet_mass      = array( 'f', self.maxn*[ 0. ] )

         self.gamma_size       = array( 'i', [ 0 ] )
         self.gamma_pt         = array( 'f', self.maxn*[ 0. ] )
         self.gamma_eta        = array( 'f', self.maxn*[ 0. ] )
         self.gamma_phi        = array( 'f', self.maxn*[ 0. ] )
         self.gamma_mass       = array( 'f', self.maxn*[ 0. ] )
         self.gamma_idvar      = array( 'f', self.maxn*[ 0. ] )
         self.gamma_reliso     = array( 'f', self.maxn*[ 0. ] )
         self.gamma_idpass     = array( 'i', self.maxn*[ 0 ] )
         self.gamma_isopass    = array( 'i', self.maxn*[ 0 ] )

         self.elec_size        = array( 'i', [ 0 ] )
         self.elec_pt          = array( 'f', self.maxn*[ 0. ] )
         self.elec_eta         = array( 'f', self.maxn*[ 0. ] )
         self.elec_phi         = array( 'f', self.maxn*[ 0. ] )
         self.elec_mass        = array( 'f', self.maxn*[ 0. ] )
         self.elec_charge      = array( 'i', self.maxn*[ 0 ] )
         self.elec_idvar       = array( 'f', self.maxn*[ 0. ] )
         self.elec_reliso      = array( 'f', self.maxn*[ 0. ] )
         self.elec_idpass      = array( 'i', self.maxn*[ 0 ] )
         self.elec_isopass     = array( 'i', self.maxn*[ 0 ] )

         self.muon_size        = array( 'i', [ 0 ] )
         self.muon_pt          = array( 'f', self.maxn*[ 0. ] )
         self.muon_eta         = array( 'f', self.maxn*[ 0. ] )
         self.muon_phi         = array( 'f', self.maxn*[ 0. ] )
         self.muon_mass        = array( 'f', self.maxn*[ 0. ] )
         self.muon_charge      = array( 'i', self.maxn*[ 0 ] )
         self.muon_idvar       = array( 'f', self.maxn*[ 0. ] )
         self.muon_reliso      = array( 'f', self.maxn*[ 0. ] )
         self.muon_idpass      = array( 'i', self.maxn*[ 0 ] )
         self.muon_isopass     = array( 'i', self.maxn*[ 0 ] )

         self.tau_size         = array( 'i', [ 0 ] )
         self.tau_pt           = array( 'f', self.maxn*[ 0. ] )
         self.tau_eta          = array( 'f', self.maxn*[ 0. ] )
         self.tau_phi          = array( 'f', self.maxn*[ 0. ] )
         self.tau_mass         = array( 'f', self.maxn*[ 0. ] )
         self.tau_charge       = array( 'i', self.maxn*[ 0 ] )
         self.tau_decaymode    = array( 'f', self.maxn*[ 0. ] )
         self.tau_reliso       = array( 'f', self.maxn*[ 0. ] )
         self.tau_isofunction  = array( 'f', self.maxn*[ 0. ] )
         self.tau_isopass      = array( 'i', self.maxn*[ 0 ] )

         self.jet_size         = array( 'i', [ 0 ] )
         self.jet_pt           = array( 'f', self.maxn*[ 0. ] )
         self.jet_eta          = array( 'f', self.maxn*[ 0. ] )
         self.jet_phi          = array( 'f', self.maxn*[ 0. ] )
         self.jet_mass         = array( 'f', self.maxn*[ 0. ] )
         self.jet_idpass       = array( 'i', self.maxn*[ 0 ] )
         self.jet_DeepJET      = array( 'f', self.maxn*[ 0. ] )
         self.jet_btag         = array( 'i', self.maxn*[ 0 ] )

         self.met_size         = array( 'i', [ 0 ] )
         self.met_pt           = array( 'f', self.maxn*[ 0. ] )
         self.met_phi          = array( 'f', self.maxn*[ 0. ] )

         # declare tree branches
         self.t.Branch( "evt_size",self.evt_size, "evt_size/I")

         self.t.Branch( "vtx_size",self.vtx_size, "vtx_size/I")
         self.t.Branch( "vtx_pt2",self.vtx_pt2, "vtx_pt2[vtx_size]/F")

         self.t.Branch( "genpart_size",self.genpart_size, "genpart_size/I")
         self.t.Branch( "genpart_pid",self. genpart_pid, "genpart_pid[genpart_size]/I")
         self.t.Branch( "genpart_status",self.genpart_status, "genpart_status[genpart_size]/I")
         self.t.Branch( "genpart_pt",self.genpart_pt, "genpart_pt[genpart_size]/F")
         self.t.Branch( "genpart_eta",self.genpart_eta, "genpart_eta[genpart_size]/F")
         self.t.Branch( "genpart_phi",self.genpart_phi, "genpart_phi[genpart_size]/F")
         self.t.Branch( "genpart_mass",self.genpart_mass, "genpart_mass[genpart_size]/F")
         self.t.Branch( "genpart_m1",self.genpart_m1, "genpart_m1[genpart_size]/I")
         self.t.Branch( "genpart_m2",self.genpart_m2, "genpart_m2[genpart_size]/I")
         self.t.Branch( "genpart_d1",self.genpart_d1, "genpart_d1[genpart_size]/I")
         self.t.Branch( "genpart_d2",self.genpart_d2, "genpart_d2[genpart_size]/I")

         self.t.Branch( "genjet_size",self.genjet_size, "genjet_size/I")
         self.t.Branch( "genjet_pt",self.genjet_pt, "genjet_pt[genjet_size]/F")
         self.t.Branch( "genjet_eta",self.genjet_eta, "genjet_eta[genjet_size]/F")
         self.t.Branch( "genjet_phi",self.genjet_phi, "genjet_phi[genjet_size]/F")
         self.t.Branch( "genjet_mass",self.genjet_mass, "genjet_mass[genjet_size]/F")

         self.t.Branch( "gamma_size",self.gamma_size, "gamma_size/I")
         self.t.Branch( "gamma_pt",self.gamma_pt, "gamma_pt[gamma_size]/F")
         self.t.Branch( "gamma_eta",self.gamma_eta, "gamma_eta[gamma_size]/F")
         self.t.Branch( "gamma_phi",self.gamma_phi, "gamma_phi[gamma_size]/F")
         self.t.Branch( "gamma_mass",self.gamma_mass, "gamma_mass[gamma_size]/F")
         self.t.Branch( "gamma_idvar",self. gamma_idvar, "gamma_idvar[gamma_size]/F")
         self.t.Branch( "gamma_reliso",self.gamma_reliso, "gamma_reliso[gamma_size]/F")
         self.t.Branch( "gamma_idpass",self. gamma_idpass, "gamma_idpass[gamma_size]/i")
         self.t.Branch( "gamma_isopass",self. gamma_isopass, "gamma_isopass[gamma_size]/i")

         self.t.Branch( "elec_size",self.elec_size, "elec_size/I")
         self.t.Branch( "elec_pt",self.elec_pt, "elec_pt[elec_size]/F")
         self.t.Branch( "elec_eta",self.elec_eta, "elec_eta[elec_size]/F")
         self.t.Branch( "elec_phi",self.elec_phi, "elec_phi[elec_size]/F")
         self.t.Branch( "elec_mass",self.elec_mass, "elec_mass[elec_size]/F")
         self.t.Branch( "elec_charge",self.elec_charge, "elec_charge[elec_size]/I")
         self.t.Branch( "elec_idvar",self. elec_idvar, "elec_idvar[elec_size]/F")
         self.t.Branch( "elec_reliso",self.elec_reliso, "elec_reliso[elec_size]/F")
         self.t.Branch( "elec_idpass",self.elec_idpass, "elec_idpass[elec_size]/i")
         self.t.Branch( "elec_isopass",self. elec_isopass, "elec_isopass[elec_size]/i")

         self.t.Branch( "muon_size",self.muon_size, "muon_size/I")
         self.t.Branch( "muon_pt",self.muon_pt, "muon_pt[muon_size]/F")
         self.t.Branch( "muon_eta",self.muon_eta, "muon_eta[muon_size]/F")
         self.t.Branch( "muon_phi",self.muon_phi, "muon_phi[muon_size]/F")
         self.t.Branch( "muon_mass",self.muon_mass, "muon_mass[muon_size]/F")
         self.t.Branch( "muon_charge",self.muon_charge, "muon_charge[muon_size]/I")
         self.t.Branch( "muon_idvar",self. muon_idvar, "muon_idvar[muon_size]/F")
         self.t.Branch( "muon_reliso",self.muon_reliso, "muon_reliso[muon_size]/F")
         self.t.Branch( "muon_idpass",self. muon_idpass, "muon_idpass[muon_size]/i")
         self.t.Branch( "muon_isopass",self. muon_isopass, "muon_isopass[muon_size]/i")

         self.t.Branch( "tau_size",self.tau_size, "tau_size/I")
         self.t.Branch( "tau_pt",self.tau_pt, "tau_pt[tau_size]/F")
         self.t.Branch( "tau_eta",self.tau_eta, "tau_eta[tau_size]/F")
         self.t.Branch( "tau_phi",self.tau_phi, "tau_phi[tau_size]/F")
         self.t.Branch( "tau_mass",self.tau_mass, "tau_mass[tau_size]/F")
         self.t.Branch( "tau_charge",self.tau_charge, "tau_charge[tau_size]/I")
         self.t.Branch( "tau_decaymode",self.tau_decaymode, "tau_decaymode[tau_size]/F")
         self.t.Branch( "tau_reliso",self.tau_reliso, "tau_reliso[tau_size]/F")
         self.t.Branch( "tau_isofunction",self.tau_isofunction, "tau_isofunction[tau_size]/F")
         self.t.Branch( "tau_isopass",self. tau_isopass, "tau_isopass[tau_size]/i")

         self.t.Branch( "jet_size",self.jet_size, "jet_size/I")
         self.t.Branch( "jet_pt",self.jet_pt, "jet_pt[jet_size]/F")
         self.t.Branch( "jet_eta",self.jet_eta, "jet_eta[jet_size]/F")
         self.t.Branch( "jet_phi",self.jet_phi, "jet_phi[jet_size]/F")
         self.t.Branch( "jet_mass",self.jet_mass, "jet_mass[jet_size]/F")
         self.t.Branch( "jet_idpass",self. jet_idpass, "jet_idpass[jet_size]/i")
         self.t.Branch( "jet_DeepJET",self.jet_DeepJET,"jet_DeepJET[jet_size]/F")
         self.t.Branch( "jet_btag",self.jet_btag,"jet_btag[jet_size]/I")

         self.t.Branch( "met_size",self.met_size, "met_size/I")
         self.t.Branch( "met_pt",self. met_pt, "met_pt[met_size]/F")
         self.t.Branch( "met_phi",self.met_phi, "met_phi[met_size]/F")

    #___________________________________________
    def processEvent(self, entry):
        self.evt_size[0] = entry

    #___________________________________________
    def processVertices(self, vertices):
        i = 0
        for item in vertices:
            self.vtx_pt2[i] = item.SumPT2
            i += 1
        self.vtx_size[0] = i

    #___________________________________________
    def processGenParticles(self, particles):
        i = 0
        for item in particles:
            self.genpart_pid    [i] = item.PID
            self.genpart_status [i] = item.Status
            self.genpart_pt     [i] = item.PT
            self.genpart_eta    [i] = item.Eta
            self.genpart_phi    [i] = item.Phi
            self.genpart_mass   [i] = item.Mass
            self.genpart_m1     [i] = item.M1
            self.genpart_m2     [i] = item.M2
            self.genpart_d1     [i] = item.D1
            self.genpart_d2     [i] = item.D2
            i += 1
        self.genpart_size[0] = i


    #___________________________________________
    def processGenJets(self, genjets):
        i = 0
        for item in genjets:
            self.genjet_pt     [i] = item.PT
            self.genjet_eta    [i] = item.Eta
            self.genjet_phi    [i] = item.Phi
            self.genjet_mass   [i] = item.Mass
            i += 1
        self.genjet_size[0] = i


    #___________________________________________
    def processPhotons(self, photons, photons_loose, photons_medium, photons_tight):
        i = 0
        for item in photons:
            self.gamma_pt      [i] = item.PT
            self.gamma_eta     [i] = item.Eta
            self.gamma_phi     [i] = item.Phi
            self.gamma_mass    [i] = 0.
            self.gamma_idvar   [i] = 0.            # DUMMY
            self.gamma_reliso  [i] = item.IsolationVar
            self.gamma_idpass  [i] = 0             # DUMMY 
            self.gamma_isopass [i] = 0             # DUMMY

            if self.gamma_reliso[i] < 0.1:
               self.gamma_isopass[i] |= 1 << 0

            if self.gamma_reliso[i] < 0.2:
               self.gamma_isopass[i] |= 1 << 1

            if self.gamma_reliso[i] < 0.3:
               self.gamma_isopass[i] |= 1 << 2

            if self.gamma_reliso[i] < 0.4:
               self.gamma_isopass[i] |= 1 << 3

            # loop over ID collections
            for loose in photons_loose:
                if dr_match(item,loose,0.2):
                    self.gamma_idpass[i] |= 1 << 0

            for medium in photons_medium:
                if dr_match(item,medium,0.2):
                    self.gamma_idpass[i] |= 1 << 1

            for tight in photons_tight:
                if dr_match(item,tight,0.2):
                    self.gamma_idpass[i] |= 1 << 2
            i += 1

        self.gamma_size[0] = i


    #___________________________________________
    def processElectrons(self, electrons, electrons_loose, electrons_medium, electrons_tight):
        i = 0
        for item in electrons:
            self.elec_pt      [i] = item.PT
            self.elec_eta     [i] = item.Eta
            self.elec_phi     [i] = item.Phi
            self.elec_mass    [i] = item.P4().M()
            self.elec_charge  [i] = item.Charge
            self.elec_idvar   [i] = 0.            # DUMMY
            self.elec_reliso  [i] = item.IsolationVar
            self.elec_idpass  [i] = 0             # DUMMY 
            self.elec_isopass [i] = 0             # DUMMY

            if self.elec_reliso[i] < 0.1:
               self.elec_isopass[i] |= 1 << 0

            if self.elec_reliso[i] < 0.2:
               self.elec_isopass[i] |= 1 << 1

            if self.elec_reliso[i] < 0.3:
               self.elec_isopass[i] |= 1 << 2

            if self.elec_reliso[i] < 0.4:
               self.elec_isopass[i] |= 1 << 3

            # loop over ID collections
            for loose in electrons_loose:
                if dr_match(item,loose,0.2):
                    self.elec_idpass[i] |= 1 << 0

            for medium in electrons_medium:
                if dr_match(item,medium,0.2):
                    self.elec_idpass[i] |= 1 << 1

            for tight in electrons_tight:
                if dr_match(item,tight,0.2):
                    self.elec_idpass[i] |= 1 << 2

            i += 1
        self.elec_size[0] = i


    #___________________________________________
    def processMuons(self, muons, muons_loose, muons_medium, muons_tight):
        i = 0
        for item in muons:
            self.muon_pt      [i] = item.PT
            self.muon_eta     [i] = item.Eta
            self.muon_phi     [i] = item.Phi
            self.muon_mass    [i] = item.P4().M()
            self.muon_charge  [i] = item.Charge
            self.muon_idvar   [i] = 0.            # DUMMY
            self.muon_reliso  [i] = item.IsolationVar
            self.muon_idpass  [i] = 0             # DUMMY 
            self.muon_isopass [i] = 0             # DUMMY

            if self.muon_reliso[i] < 0.1:
               self.muon_isopass[i] |= 1 << 0

            if self.muon_reliso[i] < 0.2:
               self.muon_isopass[i] |= 1 << 1

            if self.muon_reliso[i] < 0.3:
              self.muon_isopass[i] |= 1 << 2

            if self.muon_reliso[i] < 0.4:
               self.muon_isopass[i] |= 1 << 3

            # loop over ID collections
            for loose in muons_loose:
                if dr_match(item,loose,0.2):
                    self.muon_idpass[i] |= 1 << 0

            for medium in muons_medium:
                if dr_match(item,medium,0.2):
                    self.muon_idpass[i] |= 1 << 1

            for tight in muons_tight:
                if dr_match(item,tight,0.2):
                    self.muon_idpass[i] |= 1 << 2

            i += 1
        self.muon_size[0] = i

    #___________________________________________
    def processJets(self, jets):

        i = 0
        for item in jets:
            jetp4 = item.P4()
            self.jet_pt      [i] = jetp4.Pt()
            self.jet_eta     [i] = jetp4.Eta()
            self.jet_phi     [i] = jetp4.Phi()
            self.jet_mass    [i] = jetp4.M()
            self.jet_idpass  [i] = 0             # DUMMY 

            ### JETID: Jet constituents seem to broken!! For now set all Jet ID to True TO BE FIXED ######

            # compute jet id by looping over jet constituents
            if self.debug : print '   new jet: ', item.PT, item.Eta, item.Phi, item.Mass

            p4tot = ROOT.TLorentzVector(0., 0., 0., 0.)

            if self.debug: print '   -> Nconst: ', len(item.Constituents)
            for j in xrange(len(item.Constituents)):
                const = item.Constituents.At(j)
                p4 = ROOT.TLorentzVector(0., 0., 0., 0.)
                if isinstance(const, ROOT.Track):
                    p4 = ROOT.Track(const).P4()
                    if self.debug: print '       Track: ', p4.Pt(), p4.Eta(), p4.Phi(), p4.M()

                if isinstance(const, ROOT.Tower):
                    p4 = ROOT.Tower(const).P4()
                    if self.debug: print '       Tower: ', p4.Pt(), p4.Eta(), p4.Phi(), p4.M()           

                if isinstance(const, ROOT.Muon):
                    p4 = ROOT.Muon(const).P4()
                    if self.debug: print '       Muon: ', p4.Pt(), p4.Eta(), p4.Phi(), p4.M()            
                p4tot += p4

            if self.debug : print '   jet const sum: ', p4tot.Pt(), p4tot.Eta(), p4tot.Phi(), p4tot.M()
            if self.debug : print '   jet          : ', jetp4.Pt(), jetp4.Eta(), jetp4.Phi(), jetp4.M()

            self.jet_idpass[i] |= 1 << 0
            self.jet_idpass[i] |= 1 << 1
            self.jet_idpass[i] |= 1 << 2

            #### BTagging

            self.jet_DeepJET [i] = 0.  ## some dummy value

            for j in range(3):
                if ( item.BTag & (1 << j) ):
                    self.jet_btag[i] |= 1 << j

            i += 1
        self.jet_size[0] = i



    #___________________________________________
    def processTaus(self, jets):
        ##### For now un-tagged tau's are simply jets for now (TO BE FIXED)

        i = 0
        for item in jets:

            jetp4 = item.P4()
            self.tau_pt         [i]  = jetp4.Pt()
            self.tau_eta        [i]  = jetp4.Eta()
            self.tau_phi        [i]  = jetp4.Phi()
            self.tau_mass       [i]  = jetp4.M()
            self.tau_charge     [i]  = item.Charge
            self.tau_decaymode  [i]  = 0. # dummy for now, has to be implemented in Delphes
            self.tau_reliso     [i]  = 0. # dummy for now, has to be implemented in Delphes
            self.tau_isofunction[i]  = 0. # dummy for now, has to be implemented in Delphes

            ## Tau-Tagging (for now only cut based ID, corresponds to 1st BIT in Delphes) ### TO BE FIXED !!!
            if ( item.TauTag & (1 << 0) ):
                self.tau_isopass    [i] |= 1 << 0

            i += 1
        self.tau_size[0] = i



    #___________________________________________
    def processMissingET(self, met):
        i = 0
        for item in met:

            self.met_pt    [i] = item.MET
            self.met_phi   [i] = item.Phi
            i += 1
        self.met_size  [0] = i


    def fill(self):
        self.t.Fill()

    def write(self):
        self.t.Write()

#_______________________________________________________
def dr_match(p1, p2, drmin):
    dr = p1.P4().DeltaR(p2.P4())
    return dr < drmin


#_____________________________________________________________________________________________________________
def main():

    ROOT.gSystem.Load("libDelphes")
    try:
      ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
      ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')
    except:
      pass

    parser = argparse.ArgumentParser()
    parser.add_argument ('-i', '--input', help='input Delphes file',  default='delphes.root')
    parser.add_argument ('-o', '--output', help='output flat tree',  default='tree.root')
    parser.add_argument ('-n', '--nev', help='number of events', type=int, default=-1)
    parser.add_argument ('-d', '--debug', help='debug flag',  action='store_true',  default=False)

    args = parser.parse_args()

    inputFile = args.input
    outputFile = args.output
    nevents = args.nev
    debug = args.debug

    chain = ROOT.TChain("Delphes")
    chain.Add(inputFile)
    

    # Create object of class ExRootTreeReader
    treeReader = ROOT.ExRootTreeReader(chain)
    numberOfEntries = treeReader.GetEntries()

    ## for now only M for electrons, LT for muons and LT for photons are defined !!
    ## should dervie new parameterisations for other working points

    branchVertex          = treeReader.UseBranch('Vertex')   
    branchParticle        = treeReader.UseBranch('Particle') 
    branchGenJet          = treeReader.UseBranch('GenJet')   

    branchPhoton          = treeReader.UseBranch('Photon')
    branchPhotonLoose     = treeReader.UseBranch('PhotonLoose')
    # TO BE FIXED (replace by Medium when available)!!!
    branchPhotonMedium    = branchPhotonLoose
    branchPhotonTight     = treeReader.UseBranch('PhotonTight')

    branchElectron        = treeReader.UseBranch('Electron')
    # TO BE FIXED (replace by Loose when available)!!!
    branchElectronMedium  = treeReader.UseBranch('ElectronMedium')
    branchElectronLoose   = branchElectronMedium
    branchElectronTight   = branchElectronMedium
    # TO BE FIXED (replace by Tight when available)!!!


    branchMuon            = treeReader.UseBranch('Muon')
    branchMuonLoose       = treeReader.UseBranch('MuonLoose')
    # TO BE FIXED (replace by Medium when available)!!!
    branchMuonMedium      = branchMuonLoose
    branchMuonTight       = treeReader.UseBranch('MuonTight')

    branchJet            = treeReader.UseBranch('JetPUPPI')
    #branchJet            = treeReader.UseBranch('Jet')
    branchMissingET       = treeReader.UseBranch('PuppiMissingET')

    # NEED these branches to access jet constituents
    branchEFlowTrack = treeReader.UseBranch('EFlowTrack');
    branchEFlowPhoton = treeReader.UseBranch('EFlowPhoton');
    branchEFlowNeutralHadron = treeReader.UseBranch('EFlowNeutralHadron');

    treeProducer = TreeProducer(debug)

    if nevents > 0:
        numberOfEntries = nevents

    ################ Start event loop #######################
    for entry in range(0, numberOfEntries):

        # Load selected branches with data from specified event
        treeReader.ReadEntry(entry)

        if (entry+1)%1 == 0:
            print ' ... processed {} events ...'.format(entry+1)

        treeProducer.processEvent(entry)
        treeProducer.processVertices(branchVertex)
        treeProducer.processGenParticles(branchParticle)
        treeProducer.processGenJets(branchGenJet)
        treeProducer.processElectrons(branchElectron, branchElectronLoose, branchElectronMedium, branchElectronTight)
        treeProducer.processMuons(branchMuon, branchMuonLoose, branchMuonMedium, branchMuonTight)
        treeProducer.processPhotons(branchPhoton, branchPhotonLoose, branchPhotonMedium, branchPhotonTight)
        treeProducer.processJets(branchJet)
        treeProducer.processTaus(branchJet)
        treeProducer.processMissingET(branchMissingET)

        ## fill tree 
        treeProducer.fill()

    out_root = ROOT.TFile(outputFile,"RECREATE")
    out_root.mkdir("myana")
    out_root.cd("myana")
    treeProducer.write()
 

#_______________________________________________________________________________________
if __name__ == "__main__":
    main()

