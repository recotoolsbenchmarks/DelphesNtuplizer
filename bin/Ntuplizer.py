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
         
         ## put dummy value for now
         self.true_int         = array( 'i', [ -1 ] )

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

         self.genmet_size      = array( 'i', [ 0 ] )
         self.genmet_pt        = array( 'f', self.maxn*[ 0. ] )
         self.genmet_phi       = array( 'f', self.maxn*[ 0. ] )

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
         self.tau_neutraliso   = array( 'f', self.maxn*[ 0. ] )
         self.tau_chargediso   = array( 'f', self.maxn*[ 0. ] )
         self.tau_combinediso  = array( 'f', self.maxn*[ 0. ] )
         self.tau_isopass      = array( 'i', self.maxn*[ 0 ] )

         self.jetpuppi_size         = array( 'i', [ 0 ] )
         self.jetpuppi_pt           = array( 'f', self.maxn*[ 0. ] )
         self.jetpuppi_eta          = array( 'f', self.maxn*[ 0. ] )
         self.jetpuppi_phi          = array( 'f', self.maxn*[ 0. ] )
         self.jetpuppi_mass         = array( 'f', self.maxn*[ 0. ] )
         self.jetpuppi_idpass       = array( 'i', self.maxn*[ 0 ] )
         self.jetpuppi_DeepJET      = array( 'f', self.maxn*[ 0. ] )
         self.jetpuppi_btag         = array( 'i', self.maxn*[ 0 ] )

         self.jetchs_size         = array( 'i', [ 0 ] )
         self.jetchs_pt           = array( 'f', self.maxn*[ 0. ] )
         self.jetchs_eta          = array( 'f', self.maxn*[ 0. ] )
         self.jetchs_phi          = array( 'f', self.maxn*[ 0. ] )
         self.jetchs_mass         = array( 'f', self.maxn*[ 0. ] )
         self.jetchs_idpass       = array( 'i', self.maxn*[ 0 ] )
         self.jetchs_DeepJET      = array( 'f', self.maxn*[ 0. ] )
         self.jetchs_btag         = array( 'i', self.maxn*[ 0 ] )

         self.metpuppi_size         = array( 'i', [ 0 ] )
         self.metpuppi_pt           = array( 'f', self.maxn*[ 0. ] )
         self.metpuppi_phi          = array( 'f', self.maxn*[ 0. ] )

         self.metpf_size         = array( 'i', [ 0 ] )
         self.metpf_pt           = array( 'f', self.maxn*[ 0. ] )
         self.metpf_phi          = array( 'f', self.maxn*[ 0. ] )

         # declare tree branches
         self.t.Branch( "evt_size",self.evt_size, "evt_size/I")

         self.t.Branch( "vtx_size",self.vtx_size, "vtx_size/I")
         self.t.Branch( "trueInteractions",self.true_int, "trueInteractions/I")
         self.t.Branch( "npuVertices",self.vtx_size, "npuVertices/I")
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

         self.t.Branch( "genmet_size",self.genmet_size, "genmet_size/I")
         self.t.Branch( "genmet_pt",self.genmet_pt, "genmet_pt[genmet_size]/F")
         self.t.Branch( "genmet_phi",self.genmet_phi, "genmet_phi[genmet_size]/F")

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
         self.t.Branch( "tau_neutraliso",self.tau_neutraliso, "tau_neutraliso[tau_size]/F")
         self.t.Branch( "tau_chargediso",self.tau_chargediso, "tau_chargediso[tau_size]/F")
         self.t.Branch( "tau_combinediso",self.tau_combinediso, "tau_combinediso[tau_size]/F")
         self.t.Branch( "tau_isopass",self. tau_isopass, "tau_isopass[tau_size]/i")

         self.t.Branch( "jetpuppi_size",self.jetpuppi_size, "jetpuppi_size/I")
         self.t.Branch( "jetpuppi_pt",self.jetpuppi_pt, "jetpuppi_pt[jetpuppi_size]/F")
         self.t.Branch( "jetpuppi_eta",self.jetpuppi_eta, "jetpuppi_eta[jetpuppi_size]/F")
         self.t.Branch( "jetpuppi_phi",self.jetpuppi_phi, "jetpuppi_phi[jetpuppi_size]/F")
         self.t.Branch( "jetpuppi_mass",self.jetpuppi_mass, "jetpuppi_mass[jetpuppi_size]/F")
         self.t.Branch( "jetpuppi_idpass",self. jetpuppi_idpass, "jetpuppi_idpass[jetpuppi_size]/i")
         self.t.Branch( "jetpuppi_DeepJET",self.jetpuppi_DeepJET,"jetpuppi_DeepJET[jetpuppi_size]/F")
         self.t.Branch( "jetpuppi_btag",self.jetpuppi_btag,"jetpuppi_btag[jetpuppi_size]/I")

         self.t.Branch( "jetchs_size",self.jetchs_size, "jetchs_size/I")
         self.t.Branch( "jetchs_pt",self.jetchs_pt, "jetchs_pt[jetchs_size]/F")
         self.t.Branch( "jetchs_eta",self.jetchs_eta, "jetchs_eta[jetchs_size]/F")
         self.t.Branch( "jetchs_phi",self.jetchs_phi, "jetchs_phi[jetchs_size]/F")
         self.t.Branch( "jetchs_mass",self.jetchs_mass, "jetchs_mass[jetchs_size]/F")
         self.t.Branch( "jetchs_idpass",self. jetchs_idpass, "jetchs_idpass[jetchs_size]/i")
         self.t.Branch( "jetchs_DeepJET",self.jetchs_DeepJET,"jetchs_DeepJET[jetchs_size]/F")
         self.t.Branch( "jetchs_btag",self.jetchs_btag,"jetchs_btag[jetchs_size]/I")

         self.t.Branch( "metpuppi_size",self.metpuppi_size, "metpuppi_size/I")
         self.t.Branch( "metpuppi_pt",self. metpuppi_pt, "metpuppi_pt[metpuppi_size]/F")
         self.t.Branch( "metpuppi_phi",self.metpuppi_phi, "metpuppi_phi[metpuppi_size]/F")

         self.t.Branch( "metpf_size",self.metpf_size, "metpf_size/I")
         self.t.Branch( "metpf_pt",self. metpf_pt, "metpf_pt[metpf_size]/F")
         self.t.Branch( "metpf_phi",self.metpf_phi, "metpf_phi[metpf_size]/F")

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
    def processGenMissingET(self, met):
        i = 0
        for item in met:

            self.genmet_pt    [i] = item.MET
            self.genmet_phi   [i] = item.Phi
            i += 1
        self.genmet_size  [0] = i

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
               self.gamma_isopass[i] |= 1 << 2

            if self.gamma_reliso[i] < 0.2:
               self.gamma_isopass[i] |= 1 << 1

            if self.gamma_reliso[i] < 0.3:
               self.gamma_isopass[i] |= 1 << 0

            #if self.gamma_reliso[i] < 0.4:
            #   self.gamma_isopass[i] |= 1 << 3

            # loop over ID collections
            for loose in photons_loose:
                if dr_match(item,loose,0.005):
                    self.gamma_idpass[i] |= 1 << 0

            for medium in photons_medium:
                if dr_match(item,medium,0.005):
                    self.gamma_idpass[i] |= 1 << 1

            for tight in photons_tight:
                if dr_match(item,tight,0.005):
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
               self.elec_isopass[i] |= 1 << 2

            if self.elec_reliso[i] < 0.2:
               self.elec_isopass[i] |= 1 << 1

            if self.elec_reliso[i] < 0.3:
               self.elec_isopass[i] |= 1 << 0

            #if self.elec_reliso[i] < 0.4:
            #   self.elec_isopass[i] |= 1 << 3

            # loop over ID collections
            for loose in electrons_loose:
                if dr_match(item,loose,0.005):
                    self.elec_idpass[i] |= 1 << 0

            for medium in electrons_medium:
                if dr_match(item,medium,0.005):
                    self.elec_idpass[i] |= 1 << 1

            for tight in electrons_tight:
                if dr_match(item,tight,0.005):
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

            if self.muon_reliso[i] < 0.15:
               self.muon_isopass[i] |= 1 << 2

            if self.muon_reliso[i] < 0.20:
               self.muon_isopass[i] |= 1 << 1

            if self.muon_reliso[i] < 0.25:
              self.muon_isopass[i] |= 1 << 0

            #if self.muon_reliso[i] < 0.4:
            #   self.muon_isopass[i] |= 1 << 0

            # loop over ID collections
            for loose in muons_loose:
                if dr_match(item,loose,0.005):
                    self.muon_idpass[i] |= 1 << 0

            for medium in muons_medium:
                if dr_match(item,medium,0.005):
                    self.muon_idpass[i] |= 1 << 1

            for tight in muons_tight:
                if dr_match(item,tight,0.005):
                    self.muon_idpass[i] |= 1 << 2

            i += 1
        self.muon_size[0] = i

    #___________________________________________
    def processPuppiJets(self, jets):

        i = 0

        for item in jets:
            jetp4 = item.P4()
            self.jetpuppi_pt      [i] = jetp4.Pt()
            self.jetpuppi_eta     [i] = jetp4.Eta()
            self.jetpuppi_phi     [i] = jetp4.Phi()
            self.jetpuppi_mass    [i] = jetp4.M()
            self.jetpuppi_idpass  [i] = 0             # DUMMY 

            ### JETID: Jet constituents seem to broken!! For now set all Jet ID to True TO BE FIXED ######

            # compute jet id by looping over jet constituents
            if self.debug : print '   new PUPPI jet: ', item.PT, item.Eta, item.Phi, item.Mass

            p4tot = ROOT.TLorentzVector(0., 0., 0., 0.)

            nconst = 0

            neutralEmEnergy = 0.
            chargedEmEnergy = 0.

            neutralHadEnergy = 0.
            chargedHadEnergy = 0.

            for j in xrange(len(item.Constituents)):
                const = item.Constituents.At(j)
                p4 = ROOT.TLorentzVector(0., 0., 0., 0.)
                if isinstance(const, ROOT.ParticleFlowCandidate):
                    p4 = ROOT.ParticleFlowCandidate(const).P4()
                    nconst +=1
                    p4tot += p4

                    if const.Charge == 0:
                        if const.PID == 22: neutralEmEnergy += const.E
                        if const.PID == 0 : neutralHadEnergy += const.E
                    else:
                        if abs(const.PID) == 11: chargedEmEnergy += const.E
                        elif abs(const.PID) != 13: chargedHadEnergy += const.E
                    
                    if self.debug: print '       PFCandidate: ',const.PID, p4.Pt(), p4.Eta(), p4.Phi(), p4.M()

            EmEnergy  = neutralEmEnergy + chargedEmEnergy
            HadEnergy = neutralHadEnergy + chargedHadEnergy
            
            if EmEnergy > 0.:
                neutralEmEF = neutralEmEnergy / EmEnergy
            else:
                neutralEmEF = 0.;
            
            if HadEnergy > 0.:
                neutralHadEF = neutralHadEnergy / HadEnergy
            else:
                neutralHadEF = 0.;

            if self.debug : print '   -> Nconst: ', nconst
            if self.debug : print '   jet const sum: ', p4tot.Pt(), p4tot.Eta(), p4tot.Phi(), p4tot.M()
            if self.debug : print '   jet          : ', jetp4.Pt(), jetp4.Eta(), jetp4.Phi(), jetp4.M()

            # compute jet Id (0: LOOSE, 1: MEDIUM, 2: TIGHT)

            if nconst > 1 and neutralEmEF < 0.99 and neutralHadEF < 0.99: self.jetpuppi_idpass[i] |= 1 << 0
            if nconst > 1 and neutralEmEF < 0.99 and neutralHadEF < 0.99: self.jetpuppi_idpass[i] |= 1 << 1
            if nconst > 1 and neutralEmEF < 0.90 and neutralHadEF < 0.90: self.jetpuppi_idpass[i] |= 1 << 2

            self.jetpuppi_idpass[i] |= 1 << 0
            self.jetpuppi_idpass[i] |= 1 << 1
            self.jetpuppi_idpass[i] |= 1 << 2

            #### BTagging

            self.jetpuppi_DeepJET [i] = 0.  ## some dummy value
            self.jetpuppi_btag[i] = item.BTag
            
            '''for j in range(3):
                if ( item.BTag & (1 << j) ):
                    self.jetpuppi_btag[i] |= 1 << j
                    print 'jet', i, self.jetpuppi_btag[i] 
            '''
            i += 1
        self.jetpuppi_size[0] = i

    #___________________________________________
    def processCHSJets(self, jets, rho):

        i = 0
        for item in jets:
            jetp4 = item.P4()
            self.jetchs_pt      [i] = jetp4.Pt()
            self.jetchs_eta     [i] = jetp4.Eta()
            self.jetchs_phi     [i] = jetp4.Phi()
            self.jetchs_mass    [i] = jetp4.M()
            self.jetchs_idpass  [i] = 0             # DUMMY 

            ### JETID: Jet constituents seem to broken!! For now set all Jet ID to True TO BE FIXED ######

            # compute jet id by looping over jet constituents
            if self.debug : print '   new CHS jet: ', item.PT, item.Eta, item.Phi, item.Mass

            p4tot = ROOT.TLorentzVector(0., 0., 0., 0.)

            nconst = 0

            neutralEmEnergy = 0.
            chargedEmEnergy = 0.

            neutralHadEnergy = 0.
            chargedHadEnergy = 0.

            for j in xrange(len(item.Constituents)):
                const = item.Constituents.At(j)
                p4 = ROOT.TLorentzVector(0., 0., 0., 0.)
                if isinstance(const, ROOT.ParticleFlowCandidate):
                    p4 = ROOT.ParticleFlowCandidate(const).P4()
                    nconst +=1
                    if self.debug: print '       PFCandidate: ',const.PID, p4.Pt(), p4.Eta(), p4.Phi(), p4.M()
                    p4tot += p4

                    if const.Charge == 0:
                	if const.PID == 22: neutralEmEnergy += const.E
                	if const.PID == 0 : neutralHadEnergy += const.E
                    else:
                	if abs(const.PID) == 11: chargedEmEnergy += const.E
                	elif abs(const.PID) != 13: chargedHadEnergy += const.E

            corr = TLorentzVector()
            for r in rho:
                if item.Eta > r.Edges[0] and item.Eta < r.Edges[1]:
                    corr = item.Area * r.Rho 

	    neutralEmEnergy  -= corr.E()
	    EmEnergy  = neutralEmEnergy + chargedEmEnergy
	    
	    #neutralHadEnergy  -= corr.E()
	    HadEnergy = neutralHadEnergy + chargedHadEnergy

            if EmEnergy > 0.:
                neutralEmEF = neutralEmEnergy / EmEnergy
            else:
                neutralEmEF = 0.;
            
            if HadEnergy > 0.:
                neutralHadEF = neutralHadEnergy / HadEnergy
            else:
                neutralHadEF = 0.;

             
            if self.debug: print '   -> Nconst: ', nconst
            if self.debug: print '   -> Nconst: ', nconst

            sumpTcorr = p4tot - corr
            if self.debug : print '   jet const sum uncorr. : ', p4tot.Pt(), p4tot.Eta(), p4tot.Phi(), p4tot.M()
            if self.debug : print '   jet const sum corr.   : ', sumpTcorr.Pt(), sumpTcorr.Eta(), sumpTcorr.Phi(), sumpTcorr.M()
            if self.debug : print '   jet                   : ', jetp4.Pt(), jetp4.Eta(), jetp4.Phi(), jetp4.M()

            # compute jet Id (0: LOOSE, 1: MEDIUM, 2: TIGHT)

            if nconst > 1 and neutralEmEF < 0.99 and neutralHadEF < 0.99: self.jetchs_idpass[i] |= 1 << 0
            if nconst > 1 and neutralEmEF < 0.99 and neutralHadEF < 0.99: self.jetchs_idpass[i] |= 1 << 1
            if nconst > 1 and neutralEmEF < 0.90 and neutralHadEF < 0.90: self.jetchs_idpass[i] |= 1 << 2

            i += 1
        self.jetchs_size[0] = i


    #___________________________________________
    def processTaus(self, jets):
        ##### For now un-tagged tau's are simply jets for now (TO BE FIXED)

        i = 0
        for item in jets:

            jetp4 = item.P4()
            self.tau_pt          [i]  = jetp4.Pt()
            self.tau_eta         [i]  = jetp4.Eta()
            self.tau_phi         [i]  = jetp4.Phi()
            self.tau_mass        [i]  = jetp4.M()
            self.tau_charge      [i]  = item.Charge
            self.tau_decaymode   [i]  = 0. # dummy for now, has to be implemented in Delphes
            self.tau_chargediso  [i]  = 0. # dummy for now, has to be implemented in Delphes
            self.tau_neutraliso  [i]  = 0. # dummy for now, has to be implemented in Delphes
            self.tau_combinediso [i]  = 0. # dummy for now, has to be implemented in Delphes

            #if ( item.TauTag & (1 << 0) ):
            #    self.tau_isopass    [i] |= 1 << 0

            self.tau_isopass[i] = item.TauTag

            #for j in range(4):
            #    if ( item.TauTag & (1 << j) ):
            #        self.tau_isopass[i] |= 1 << j

            i += 1
        self.tau_size[0] = i

    #___________________________________________
    def processPFMissingET(self, met):
        i = 0
        for item in met:

            self.metpf_pt    [i] = item.MET
            self.metpf_phi   [i] = item.Phi
            i += 1
        self.metpf_size  [0] = i

    #___________________________________________
    def processPuppiMissingET(self, met):
        i = 0
        for item in met:

            self.metpuppi_pt    [i] = item.MET
            self.metpuppi_phi   [i] = item.Phi
            i += 1
        self.metpuppi_size  [0] = i


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
    branchGenMissingET    = treeReader.UseBranch('GenMissingET')   

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

    branchPuppiJet        = treeReader.UseBranch('JetPUPPI')
    branchCHSJet          = treeReader.UseBranch('Jet')

    branchPuppiMissingET  = treeReader.UseBranch('PuppiMissingET')
    branchPFMissingET     = treeReader.UseBranch('MissingET')

    # NEED these branches to access jet constituents
    branchPuppiCandidate  = treeReader.UseBranch('ParticleFlowCandidate')
    branchPFCandidateCHS  = treeReader.UseBranch('ParticleFlowCandidateCHS')

    branchRho             = treeReader.UseBranch('Rho')
    
    treeProducer = TreeProducer(debug)

    if nevents > 0:
        numberOfEntries = nevents

    ################ Start event loop #######################
    for entry in range(0, numberOfEntries):

        # Load selected branches with data from specified event
        treeReader.ReadEntry(entry)

        if (entry+1)%1000 == 0:
            print ' ... processed {} events ...'.format(entry+1)

        treeProducer.processEvent(entry)
        treeProducer.processVertices(branchVertex)
        treeProducer.processGenParticles(branchParticle)
        treeProducer.processGenJets(branchGenJet)
        treeProducer.processGenMissingET(branchGenMissingET)
        treeProducer.processElectrons(branchElectron, branchElectronLoose, branchElectronMedium, branchElectronTight)
        treeProducer.processMuons(branchMuon, branchMuonLoose, branchMuonMedium, branchMuonTight)
        treeProducer.processPhotons(branchPhoton, branchPhotonLoose, branchPhotonMedium, branchPhotonTight)
        treeProducer.processCHSJets(branchCHSJet, branchRho)
        treeProducer.processPuppiJets(branchPuppiJet)
        treeProducer.processTaus(branchPuppiJet)
        treeProducer.processPFMissingET(branchPFMissingET)
        treeProducer.processPuppiMissingET(branchPuppiMissingET)

        ## fill tree 
        treeProducer.fill()

    out_root = ROOT.TFile(outputFile,"RECREATE")
    out_root.mkdir("myana")
    out_root.cd("myana")
    treeProducer.write()
 

#_______________________________________________________________________________________
if __name__ == "__main__":
    main()

