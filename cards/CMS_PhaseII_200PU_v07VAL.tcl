
########################################
#
#  Main authors: Michele Selvaggi (CERN)
#
#  Released on: Nov. 2020
#
#  Version: v07 base 
#
#  Notes: - added DNN tau tagging ( 01/04/2019)
#         - removed BTaggingMTD, and replaced BTagging with latest DeepJet parameterisation (20/07/2019)
#         - adding electrons, muons and photon fakes
#         - removed ele/mu/gamma CHS collections
#         - adding medium WPs for ele/mu/photons
#
#
#######################################
# Order of execution of various modules
#######################################
  

set ExecutionPath {

  PileUpMerger
  ParticlePropagator
  TrackMergerProp

  DenseProp
  DenseMergeTracks
  DenseTrackFilter

  ChargedHadronTrackingEfficiency
  ElectronTrackingEfficiency
  MuonTrackingEfficiency

  ChargedHadronMomentumSmearing
  ElectronEnergySmearing
  MuonMomentumSmearing

  TrackMerger

  ECal
  HCal

  PhotonEnergySmearing
  ElectronFilter

  TrackPileUpSubtractor
  RecoPuFilter

  TowerMerger
  NeutralEFlowMerger

  EFlowMerger
  EFlowMergerCHS
  Rho

  LeptonFilterNoLep
  LeptonFilterLep
  RunPUPPIBase
  RunPUPPIMerger
  RunPUPPI

  EFlowFilterPuppi

  GenParticleFilter
  PhotonFilter

  NeutrinoFilter

  PhotonIsolation
  PhotonScale
  PhotonSmear


  PhotonLooseID
  PhotonMediumID
  PhotonTightID

  ElectronIsolation
  ElectronScale
  ElectronSmear

  ElectronLooseEfficiency
  ElectronMediumEfficiency
  ElectronTightEfficiency

  MuonIsolation
  MuonScale
  MuonSmear

  MuonLooseIdEfficiency
  MuonMediumIdEfficiency
  MuonTightIdEfficiency

  GenMissingET
  GenPileUpMissingET

  GenJetFinder
  GenJetFinderAK8

  PuppiMissingET
  ScalarHT

  FastJetFinderPUPPI
  FastJetFinderPUPPIAK8

  JetScalePUPPI
  JetScalePUPPIAK8
  JetSmearPUPPI
  JetSmearPUPPIAK8

  JetLooseID
  JetTightID

  JetFlavorAssociationPUPPI
  JetFlavorAssociationPUPPIAK8

  BTaggingPUPPILoose
  BTaggingPUPPIMedium
  BTaggingPUPPITight

  BTaggingPUPPIAK8Loose
  BTaggingPUPPIAK8Medium
  BTaggingPUPPIAK8Tight

  TauTaggingPUPPILoose
  TauTaggingPUPPIMedium
  TauTaggingPUPPITight

  JetFakeMakerLoose
  JetFakeMakerMedium
  JetFakeMakerTight

  PhotonFakeMergerLoose
  PhotonFakeMergerMedium
  PhotonFakeMergerTight

  ElectronFakeMergerLoose
  ElectronFakeMergerMedium
  ElectronFakeMergerTight

  MuonFakeMergerLoose
  MuonFakeMergerMedium
  MuonFakeMergerTight

  TreeWriter
}


###############
# PileUp Merger
###############

module PileUpMerger PileUpMerger {
  set InputArray Delphes/stableParticles

  set ParticleOutputArray stableParticles
  set VertexOutputArray vertices

  # pre-generated minbias input file
  set PileUpFile /eos/cms/store/group/upgrade/delphes/PhaseII/MinBias_100k.pileup
  #set PileUpFile MinBias_100k.pileup
  
  # average expected pile up
  set MeanPileUp 200

  # maximum spread in the beam direction in m
  set ZVertexSpread 0.25

  # maximum spread in time in s
  set TVertexSpread 800E-12

  # vertex smearing formula f(z,t) (z,t need to be respectively given in m,s) - {exp(-(t^2/160e-12^2/2))*exp(-(z^2/0.053^2/2))}
  set VertexDistributionFormula {exp(-(t^2/160e-12^2/2))*exp(-(z^2/0.053^2/2))}

}



#####################################
# Track propagation to calorimeters
#####################################

module ParticlePropagator ParticlePropagator {
  set InputArray PileUpMerger/stableParticles

  set OutputArray stableParticles
  set NeutralOutputArray neutralParticles
  set ChargedHadronOutputArray chargedHadrons
  set ElectronOutputArray electrons
  set MuonOutputArray muons

  # radius of the magnetic field coverage, in m
  set Radius 1.29
  # half-length of the magnetic field coverage, in m
  set HalfLength 3.0

  # magnetic field
  set Bz 3.8
}


##############
# Track merger
##############

module Merger TrackMergerProp {
# add InputArray InputArray
  add InputArray ParticlePropagator/chargedHadrons
  add InputArray ParticlePropagator/electrons
  add InputArray ParticlePropagator/muons
  set OutputArray tracks
}


####################################
# Track propagation to pseudo-pixel
####################################

module ParticlePropagator DenseProp {

  set InputArray TrackMergerProp/tracks

  # radius of the first pixel layer
  set Radius 0.3
  set RadiusMax 1.29
  # half-length of the magnetic field coverage, in m
  set HalfLength 0.7
  set HalfLengthMax 3.0

  # magnetic field
  set Bz 3.8
}


####################
# Dense Track merger
###################

module Merger DenseMergeTracks {
# add InputArray InputArray
  add InputArray DenseProp/chargedHadrons
  add InputArray DenseProp/electrons
  add InputArray DenseProp/muons
  set OutputArray tracks
}

######################
#   Dense Track Filter
######################

module DenseTrackFilter DenseTrackFilter {

  set TrackInputArray DenseMergeTracks/tracks

  set TrackOutputArray tracks
  set ChargedHadronOutputArray chargedHadrons
  set ElectronOutputArray electrons
  set MuonOutputArray muons

  set EtaPhiRes 0.003
  set EtaMax 4.0

  set pi [expr {acos(-1)}]

  set nbins_phi [expr {$pi/$EtaPhiRes} ]
  set nbins_phi [expr {int($nbins_phi)} ]

  set PhiBins {}
  for {set i -$nbins_phi} {$i <= $nbins_phi} {incr i} {
    add PhiBins [expr {$i * $pi/$nbins_phi}]
  }

  set nbins_eta [expr {$EtaMax/$EtaPhiRes} ]
  set nbins_eta [expr {int($nbins_eta)} ]

  for {set i -$nbins_eta} {$i <= $nbins_eta} {incr i} {
    set eta [expr {$i * $EtaPhiRes}]
    add EtaPhiBins $eta $PhiBins
  }

}



####################################
# Charged hadron tracking efficiency
####################################

module Efficiency ChargedHadronTrackingEfficiency {
  ## particles after propagation
  set InputArray  DenseTrackFilter/chargedHadrons
  set OutputArray chargedHadrons
  # tracking efficiency formula for charged hadrons
  set EfficiencyFormula {
      (pt <= 0.2) * (0.00) + \
          (abs(eta) <= 1.2) * (pt > 0.2 && pt <= 1.0) * (pt * 0.96) + \
          (abs(eta) <= 1.2) * (pt > 1.0) * (0.97) + \
          (abs(eta) > 1.2 && abs(eta) <= 2.5) * (pt > 0.2 && pt <= 1.0) * (pt*0.85) + \
          (abs(eta) > 1.2 && abs(eta) <= 2.5) * (pt > 1.0) * (0.87) + \
          (abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 0.2 && pt <= 1.0) * (pt*0.8) + \
          (abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 1.0) * (0.82) + \
          (abs(eta) > 4.0) * (0.00)
  }
}


#####################################
# Electron tracking efficiency - ID
####################################

module Efficiency ElectronTrackingEfficiency {
  set InputArray  DenseTrackFilter/electrons
  set OutputArray electrons
  # tracking efficiency formula for electrons
  set EfficiencyFormula {
      (pt <= 0.2) * (0.00) + \
          (abs(eta) <= 1.2) * (pt > 0.2 && pt <= 1.0) * (pt * 0.96) + \
          (abs(eta) <= 1.2) * (pt > 1.0) * (0.97) + \
          (abs(eta) > 1.2 && abs(eta) <= 2.5) * (pt > 0.2 && pt <= 1.0) * (pt*0.85) + \
          (abs(eta) > 1.2 && abs(eta) <= 2.5) * (pt > 1.0 && pt <= 10.0) * (0.82+pt*0.01) + \
          (abs(eta) > 1.2 && abs(eta) <= 2.5) * (pt > 10.0) * (1.0) + \
          (abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 0.2 && pt <= 1.0) * (pt*0.8) + \
          (abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 1.0 && pt <= 10.0) * (0.8+pt*0.01) + \
          (abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 10.0) * (1.0) + \
          (abs(eta) > 4.0) * (0.00)

  }
  # eta 1.2-2.5 had 0.9 and 2.5-4.0 had 0.85 above 10 GeV
}

##########################
# Muon tracking efficiency
##########################

module Efficiency MuonTrackingEfficiency {
  set InputArray DenseTrackFilter/muons
  set OutputArray muons
  # tracking efficiency formula for muons
  set EfficiencyFormula {
      (pt <= 0.2) * (0.00) + \
          (abs(eta) <= 1.2) * (pt > 0.2 && pt <= 1.0) * (pt * 1.00) + \
          (abs(eta) <= 1.2) * (pt > 1.0) * (1.00) + \
          (abs(eta) > 1.2 && abs(eta) <= 2.8) * (pt > 0.2 && pt <= 1.0) * (pt*1.00) + \
          (abs(eta) > 1.2 && abs(eta) <= 2.8) * (pt > 1.0) * (1.00) + \
          (abs(eta) > 2.8 && abs(eta) <= 4.0) * (pt > 0.2 && pt <= 1.0) * (pt*0.95) + \
          (abs(eta) > 2.8 && abs(eta) <= 4.0) * (pt > 1.0) * (0.95) + \
          (abs(eta) > 4.0) * (0.00)

  }
}


########################################
# Momentum resolution for charged tracks
########################################

module MomentumSmearing ChargedHadronMomentumSmearing {
  ## hadrons after having applied the tracking efficiency
  set InputArray  ChargedHadronTrackingEfficiency/chargedHadrons
  set OutputArray chargedHadrons
  # resolution formula for charged hadrons ,

  #
  # Automatically generated tracker resolution formula for layout: OT612IT4025
  #
  #  By Unknown author on: 2017-06-30.17:03:00
  #
  set ResolutionFormula {    (abs(eta) >= 0.0000 && abs(eta) < 0.2000) * (pt >= 0.0000 && pt < 1.0000) * (0.00457888) + \
     (abs(eta) >= 0.0000 && abs(eta) < 0.2000) * (pt >= 1.0000 && pt < 10.0000) * (0.004579 + (pt-1.000000)* 0.000045) + \
     (abs(eta) >= 0.0000 && abs(eta) < 0.2000) * (pt >= 10.0000 && pt < 100.0000) * (0.004983 + (pt-10.000000)* 0.000047) + \
     (abs(eta) >= 0.0000 && abs(eta) < 0.2000) * (pt >= 100.0000) * (0.009244*pt/100.000000) + \
     (abs(eta) >= 0.2000 && abs(eta) < 0.4000) * (pt >= 0.0000 && pt < 1.0000) * (0.00505011) + \
     (abs(eta) >= 0.2000 && abs(eta) < 0.4000) * (pt >= 1.0000 && pt < 10.0000) * (0.005050 + (pt-1.000000)* 0.000033) + \
     (abs(eta) >= 0.2000 && abs(eta) < 0.4000) * (pt >= 10.0000 && pt < 100.0000) * (0.005343 + (pt-10.000000)* 0.000043) + \
     (abs(eta) >= 0.2000 && abs(eta) < 0.4000) * (pt >= 100.0000) * (0.009172*pt/100.000000) + \
     (abs(eta) >= 0.4000 && abs(eta) < 0.6000) * (pt >= 0.0000 && pt < 1.0000) * (0.00510573) + \
     (abs(eta) >= 0.4000 && abs(eta) < 0.6000) * (pt >= 1.0000 && pt < 10.0000) * (0.005106 + (pt-1.000000)* 0.000023) + \
     (abs(eta) >= 0.4000 && abs(eta) < 0.6000) * (pt >= 10.0000 && pt < 100.0000) * (0.005317 + (pt-10.000000)* 0.000042) + \
     (abs(eta) >= 0.4000 && abs(eta) < 0.6000) * (pt >= 100.0000) * (0.009077*pt/100.000000) + \
     (abs(eta) >= 0.6000 && abs(eta) < 0.8000) * (pt >= 0.0000 && pt < 1.0000) * (0.00578020) + \
     (abs(eta) >= 0.6000 && abs(eta) < 0.8000) * (pt >= 1.0000 && pt < 10.0000) * (0.005780 + (pt-1.000000)* -0.000000) + \
     (abs(eta) >= 0.6000 && abs(eta) < 0.8000) * (pt >= 10.0000 && pt < 100.0000) * (0.005779 + (pt-10.000000)* 0.000038) + \
     (abs(eta) >= 0.6000 && abs(eta) < 0.8000) * (pt >= 100.0000) * (0.009177*pt/100.000000) + \
     (abs(eta) >= 0.8000 && abs(eta) < 1.0000) * (pt >= 0.0000 && pt < 1.0000) * (0.00728723) + \
     (abs(eta) >= 0.8000 && abs(eta) < 1.0000) * (pt >= 1.0000 && pt < 10.0000) * (0.007287 + (pt-1.000000)* -0.000031) + \
     (abs(eta) >= 0.8000 && abs(eta) < 1.0000) * (pt >= 10.0000 && pt < 100.0000) * (0.007011 + (pt-10.000000)* 0.000038) + \
     (abs(eta) >= 0.8000 && abs(eta) < 1.0000) * (pt >= 100.0000) * (0.010429*pt/100.000000) + \
     (abs(eta) >= 1.0000 && abs(eta) < 1.2000) * (pt >= 0.0000 && pt < 1.0000) * (0.01045117) + \
     (abs(eta) >= 1.0000 && abs(eta) < 1.2000) * (pt >= 1.0000 && pt < 10.0000) * (0.010451 + (pt-1.000000)* -0.000051) + \
     (abs(eta) >= 1.0000 && abs(eta) < 1.2000) * (pt >= 10.0000 && pt < 100.0000) * (0.009989 + (pt-10.000000)* 0.000043) + \
     (abs(eta) >= 1.0000 && abs(eta) < 1.2000) * (pt >= 100.0000) * (0.013867*pt/100.000000) + \
     (abs(eta) >= 1.2000 && abs(eta) < 1.4000) * (pt >= 0.0000 && pt < 1.0000) * (0.01477199) + \
     (abs(eta) >= 1.2000 && abs(eta) < 1.4000) * (pt >= 1.0000 && pt < 10.0000) * (0.014772 + (pt-1.000000)* -0.000128) + \
     (abs(eta) >= 1.2000 && abs(eta) < 1.4000) * (pt >= 10.0000 && pt < 100.0000) * (0.013616 + (pt-10.000000)* 0.000035) + \
     (abs(eta) >= 1.2000 && abs(eta) < 1.4000) * (pt >= 100.0000) * (0.016800*pt/100.000000) + \
     (abs(eta) >= 1.4000 && abs(eta) < 1.6000) * (pt >= 0.0000 && pt < 1.0000) * (0.01731474) + \
     (abs(eta) >= 1.4000 && abs(eta) < 1.6000) * (pt >= 1.0000 && pt < 10.0000) * (0.017315 + (pt-1.000000)* -0.000208) + \
     (abs(eta) >= 1.4000 && abs(eta) < 1.6000) * (pt >= 10.0000 && pt < 100.0000) * (0.015439 + (pt-10.000000)* 0.000030) + \
     (abs(eta) >= 1.4000 && abs(eta) < 1.6000) * (pt >= 100.0000) * (0.018161*pt/100.000000) + \
     (abs(eta) >= 1.6000 && abs(eta) < 1.8000) * (pt >= 0.0000 && pt < 1.0000) * (0.01942025) + \
     (abs(eta) >= 1.6000 && abs(eta) < 1.8000) * (pt >= 1.0000 && pt < 10.0000) * (0.019420 + (pt-1.000000)* -0.000417) + \
     (abs(eta) >= 1.6000 && abs(eta) < 1.8000) * (pt >= 10.0000 && pt < 100.0000) * (0.015669 + (pt-10.000000)* 0.000026) + \
     (abs(eta) >= 1.6000 && abs(eta) < 1.8000) * (pt >= 100.0000) * (0.018039*pt/100.000000) + \
     (abs(eta) >= 1.8000 && abs(eta) < 2.0000) * (pt >= 0.0000 && pt < 1.0000) * (0.02201432) + \
     (abs(eta) >= 1.8000 && abs(eta) < 2.0000) * (pt >= 1.0000 && pt < 10.0000) * (0.022014 + (pt-1.000000)* -0.000667) + \
     (abs(eta) >= 1.8000 && abs(eta) < 2.0000) * (pt >= 10.0000 && pt < 100.0000) * (0.016012 + (pt-10.000000)* 0.000045) + \
     (abs(eta) >= 1.8000 && abs(eta) < 2.0000) * (pt >= 100.0000) * (0.020098*pt/100.000000) + \
     (abs(eta) >= 2.0000 && abs(eta) < 2.2000) * (pt >= 0.0000 && pt < 1.0000) * (0.02574300) + \
     (abs(eta) >= 2.0000 && abs(eta) < 2.2000) * (pt >= 1.0000 && pt < 10.0000) * (0.025743 + (pt-1.000000)* -0.001118) + \
     (abs(eta) >= 2.0000 && abs(eta) < 2.2000) * (pt >= 10.0000 && pt < 100.0000) * (0.015681 + (pt-10.000000)* 0.000051) + \
     (abs(eta) >= 2.0000 && abs(eta) < 2.2000) * (pt >= 100.0000) * (0.020289*pt/100.000000) + \
     (abs(eta) >= 2.2000 && abs(eta) < 2.4000) * (pt >= 0.0000 && pt < 1.0000) * (0.02885821) + \
     (abs(eta) >= 2.2000 && abs(eta) < 2.4000) * (pt >= 1.0000 && pt < 10.0000) * (0.028858 + (pt-1.000000)* -0.001345) + \
     (abs(eta) >= 2.2000 && abs(eta) < 2.4000) * (pt >= 10.0000 && pt < 100.0000) * (0.016753 + (pt-10.000000)* 0.000053) + \
     (abs(eta) >= 2.2000 && abs(eta) < 2.4000) * (pt >= 100.0000) * (0.021524*pt/100.000000) + \
     (abs(eta) >= 2.4000 && abs(eta) < 2.6000) * (pt >= 0.0000 && pt < 1.0000) * (0.03204812) + \
     (abs(eta) >= 2.4000 && abs(eta) < 2.6000) * (pt >= 1.0000 && pt < 10.0000) * (0.032048 + (pt-1.000000)* -0.001212) + \
     (abs(eta) >= 2.4000 && abs(eta) < 2.6000) * (pt >= 10.0000 && pt < 100.0000) * (0.021138 + (pt-10.000000)* 0.000037) + \
     (abs(eta) >= 2.4000 && abs(eta) < 2.6000) * (pt >= 100.0000) * (0.024477*pt/100.000000) + \
     (abs(eta) >= 2.6000 && abs(eta) < 2.8000) * (pt >= 0.0000 && pt < 1.0000) * (0.03950405) + \
     (abs(eta) >= 2.6000 && abs(eta) < 2.8000) * (pt >= 1.0000 && pt < 10.0000) * (0.039504 + (pt-1.000000)* -0.001386) + \
     (abs(eta) >= 2.6000 && abs(eta) < 2.8000) * (pt >= 10.0000 && pt < 100.0000) * (0.027026 + (pt-10.000000)* 0.000037) + \
     (abs(eta) >= 2.6000 && abs(eta) < 2.8000) * (pt >= 100.0000) * (0.030392*pt/100.000000) + \
     (abs(eta) >= 2.8000 && abs(eta) < 3.0000) * (pt >= 0.0000 && pt < 1.0000) * (0.04084751) + \
     (abs(eta) >= 2.8000 && abs(eta) < 3.0000) * (pt >= 1.0000 && pt < 10.0000) * (0.040848 + (pt-1.000000)* -0.001780) + \
     (abs(eta) >= 2.8000 && abs(eta) < 3.0000) * (pt >= 10.0000 && pt < 100.0000) * (0.024824 + (pt-10.000000)* 0.000029) + \
     (abs(eta) >= 2.8000 && abs(eta) < 3.0000) * (pt >= 100.0000) * (0.027445*pt/100.000000) + \
     (abs(eta) >= 3.0000 && abs(eta) < 3.2000) * (pt >= 0.0000 && pt < 1.0000) * (0.04532425) + \
     (abs(eta) >= 3.0000 && abs(eta) < 3.2000) * (pt >= 1.0000 && pt < 10.0000) * (0.045324 + (pt-1.000000)* -0.002497) + \
     (abs(eta) >= 3.0000 && abs(eta) < 3.2000) * (pt >= 10.0000 && pt < 100.0000) * (0.022851 + (pt-10.000000)* 0.000024) + \
     (abs(eta) >= 3.0000 && abs(eta) < 3.2000) * (pt >= 100.0000) * (0.025053*pt/100.000000) + \
     (abs(eta) >= 3.2000 && abs(eta) < 3.4000) * (pt >= 0.0000 && pt < 1.0000) * (0.06418925) + \
     (abs(eta) >= 3.2000 && abs(eta) < 3.4000) * (pt >= 1.0000 && pt < 10.0000) * (0.064189 + (pt-1.000000)* -0.004055) + \
     (abs(eta) >= 3.2000 && abs(eta) < 3.4000) * (pt >= 10.0000 && pt < 100.0000) * (0.027691 + (pt-10.000000)* 0.000034) + \
     (abs(eta) >= 3.2000 && abs(eta) < 3.4000) * (pt >= 100.0000) * (0.030710*pt/100.000000) + \
     (abs(eta) >= 3.4000 && abs(eta) < 3.6000) * (pt >= 0.0000 && pt < 1.0000) * (0.07682500) + \
     (abs(eta) >= 3.4000 && abs(eta) < 3.6000) * (pt >= 1.0000 && pt < 10.0000) * (0.076825 + (pt-1.000000)* -0.004510) + \
     (abs(eta) >= 3.4000 && abs(eta) < 3.6000) * (pt >= 10.0000 && pt < 100.0000) * (0.036234 + (pt-10.000000)* 0.000049) + \
     (abs(eta) >= 3.4000 && abs(eta) < 3.6000) * (pt >= 100.0000) * (0.040629*pt/100.000000) + \
     (abs(eta) >= 3.6000 && abs(eta) < 3.8000) * (pt >= 0.0000 && pt < 1.0000) * (0.09796358) + \
     (abs(eta) >= 3.6000 && abs(eta) < 3.8000) * (pt >= 1.0000 && pt < 10.0000) * (0.097964 + (pt-1.000000)* -0.005758) + \
     (abs(eta) >= 3.6000 && abs(eta) < 3.8000) * (pt >= 10.0000 && pt < 100.0000) * (0.046145 + (pt-10.000000)* 0.000069) + \
     (abs(eta) >= 3.6000 && abs(eta) < 3.8000) * (pt >= 100.0000) * (0.052345*pt/100.000000) + \
     (abs(eta) >= 3.8000 && abs(eta) < 4.0000) * (pt >= 0.0000 && pt < 1.0000) * (0.13415929) + \
     (abs(eta) >= 3.8000 && abs(eta) < 4.0000) * (pt >= 1.0000 && pt < 10.0000) * (0.134159 + (pt-1.000000)* -0.008283) + \
     (abs(eta) >= 3.8000 && abs(eta) < 4.0000) * (pt >= 10.0000 && pt < 100.0000) * (0.059612 + (pt-10.000000)* 0.000111) + \
     (abs(eta) >= 3.8000 && abs(eta) < 4.0000) * (pt >= 100.0000) * (0.069617*pt/100.000000)
  }


}


#################################
# Energy resolution for electrons
#################################

module EnergySmearing ElectronEnergySmearing {
  set InputArray ElectronTrackingEfficiency/electrons
  set OutputArray electrons

  # set ResolutionFormula {resolution formula as a function of eta and energy}

  # resolution formula for electrons

  # taking something flat in energy for now, ECAL will take over at high energy anyway.
  # inferred from hep-ex/1306.2016 and 1502.02701
  set ResolutionFormula {

                        (abs(eta) <= 1.5)  * (energy*0.028) +
    (abs(eta) > 1.5  && abs(eta) <= 1.75)  * (energy*0.037) +
    (abs(eta) > 1.75  && abs(eta) <= 2.15) * (energy*0.038) +
    (abs(eta) > 2.15  && abs(eta) <= 3.00) * (energy*0.044) +
    (abs(eta) > 3.00  && abs(eta) <= 4.00) * (energy*0.10)}

}

###############################
# Momentum resolution for muons
###############################

module MomentumSmearing MuonMomentumSmearing {
  set InputArray MuonTrackingEfficiency/muons
  set OutputArray muons
  # resolution formula for muons

  # up to |eta| < 2.8 take measurement from tracking + muon chambers
  # for |eta| > 2.8 and pT < 5.0 take measurement from tracking alone taken from
  # http://mersi.web.cern.ch/mersi/layouts/.private/Baseline_tilted_200_Pixel_1_1_1/index.html
  source muonMomentumResolution.tcl
}



##############
# Track merger
##############

module Merger TrackMerger {
# add InputArray InputArray
  add InputArray ChargedHadronMomentumSmearing/chargedHadrons
  add InputArray ElectronEnergySmearing/electrons
  add InputArray MuonMomentumSmearing/muons
  set OutputArray tracks
}


#############
#   ECAL
#############

module SimpleCalorimeter ECal {
  set ParticleInputArray ParticlePropagator/stableParticles
  set TrackInputArray TrackMerger/tracks

  set TowerOutputArray ecalTowers
  set EFlowTrackOutputArray eflowTracks
  set EFlowTowerOutputArray eflowPhotons

  set IsEcal true

  set EnergyMin 0.5
  set EnergySignificanceMin 1.0

  set SmearTowerCenter true

  set pi [expr {acos(-1)}]

  # lists of the edges of each tower in eta and phi
  # each list starts with the lower edge of the first tower
  # the list ends with the higher edged of the last tower

  # assume 0.02 x 0.02 resolution in eta,phi in the barrel |eta| < 1.5

  set PhiBins {}
  for {set i -180} {$i <= 180} {incr i} {
    add PhiBins [expr {$i * $pi/180.0}]
  }

  # 0.02 unit in eta up to eta = 1.5 (barrel)
  for {set i -85} {$i <= 86} {incr i} {
    set eta [expr {$i * 0.0174}]
    add EtaPhiBins $eta $PhiBins
  }

  # assume 0.02 x 0.02 resolution in eta,phi in the endcaps 1.5 < |eta| < 3.0 (HGCAL- ECAL)

  set PhiBins {}
  for {set i -180} {$i <= 180} {incr i} {
    add PhiBins [expr {$i * $pi/180.0}]
  }

  # 0.02 unit in eta up to eta = 3
  for {set i 1} {$i <= 84} {incr i} {
    set eta [expr { -2.958 + $i * 0.0174}]
    add EtaPhiBins $eta $PhiBins
  }

  for {set i 1} {$i <= 84} {incr i} {
    set eta [expr { 1.4964 + $i * 0.0174}]
    add EtaPhiBins $eta $PhiBins
  }

  # take present CMS granularity for HF

  # 0.175 x (0.175 - 0.35) resolution in eta,phi in the HF 3.0 < |eta| < 5.0
  set PhiBins {}
  for {set i -18} {$i <= 18} {incr i} {
    add PhiBins [expr {$i * $pi/18.0}]
  }

  foreach eta {-5 -4.7 -4.525 -4.35 -4.175 -4 -3.825 -3.65 -3.475 -3.3 -3.125 -2.958 3.125 3.3 3.475 3.65 3.825 4 4.175 4.35 4.525 4.7 5} {
    add EtaPhiBins $eta $PhiBins
  }


  add EnergyFraction {0} {0.0}
  # energy fractions for e, gamma and pi0
  add EnergyFraction {11} {1.0}
  add EnergyFraction {22} {1.0}
  add EnergyFraction {111} {1.0}
  # energy fractions for muon, neutrinos and neutralinos
  add EnergyFraction {12} {0.0}
  add EnergyFraction {13} {0.0}
  add EnergyFraction {14} {0.0}
  add EnergyFraction {16} {0.0}
  add EnergyFraction {1000022} {0.0}
  add EnergyFraction {1000023} {0.0}
  add EnergyFraction {1000025} {0.0}
  add EnergyFraction {1000035} {0.0}
  add EnergyFraction {1000045} {0.0}
  # energy fractions for K0short and Lambda
  add EnergyFraction {310} {0.3}
  add EnergyFraction {3122} {0.3}

  # set ResolutionFormula {resolution formula as a function of eta and energy}

  # for the ECAL barrel (|eta| < 1.5), see hep-ex/1306.2016 and 1502.02701
  # for the endcaps (1.5 < |eta| < 3.0), we take HGCAL  see LHCC-P-008, Fig. 3.39, p.117

  set ResolutionFormula {  (abs(eta) <= 1.50)                    * sqrt(energy^2*0.009^2 + energy*0.12^2 + 0.45^2) +
                           (abs(eta) > 1.50 && abs(eta) <= 1.75) * sqrt(energy^2*0.006^2 + energy*0.20^2) + \
                           (abs(eta) > 1.75 && abs(eta) <= 2.15) * sqrt(energy^2*0.007^2 + energy*0.21^2) + \
                           (abs(eta) > 2.15 && abs(eta) <= 3.00) * sqrt(energy^2*0.008^2 + energy*0.24^2) + \
                           (abs(eta) >= 3.0 && abs(eta) <= 5.0)  * sqrt(energy^2*0.08^2 + energy*1.98^2)}

}

#############
#   HCAL
#############

module SimpleCalorimeter HCal {
  set ParticleInputArray ParticlePropagator/stableParticles
  set TrackInputArray ECal/eflowTracks

  set TowerOutputArray hcalTowers
  set EFlowTrackOutputArray eflowTracks
  set EFlowTowerOutputArray eflowNeutralHadrons

  set IsEcal false

  set EnergyMin 1.0
  set EnergySignificanceMin 1.0

  set SmearTowerCenter true

  set pi [expr {acos(-1)}]

  # lists of the edges of each tower in eta and phi
  # each list starts with the lower edge of the first tower
  # the list ends with the higher edged of the last tower

  # assume 0.087 x 0.087 resolution in eta,phi in the barrel |eta| < 1.5

  set PhiBins {}
  for {set i -36} {$i <= 36} {incr i} {
    add PhiBins [expr {$i * $pi/36.0}]
  }
  foreach eta {-1.566 -1.479 -1.392 -1.305 -1.218 -1.131 -1.044 -0.957 -0.87 -0.783 -0.696 -0.609 -0.522 -0.435 -0.348 -0.261 -0.174 -0.087 0 0.087 0.174 0.261 0.348 0.435 0.522 0.609 0.696 0.783 0.87 0.957 1.044 1.131 1.218 1.305 1.392 1.479 1.566 1.65} {
    add EtaPhiBins $eta $PhiBins
  }

  # assume 0.07 x 0.07 resolution in eta,phi in the endcaps 1.5 < |eta| < 3.0 (HGCAL- HCAL)

  set PhiBins {}
  for {set i -45} {$i <= 45} {incr i} {
    add PhiBins [expr {$i * $pi/45.0}]
  }

  # 0.07 unit in eta up to eta = 3
  for {set i 1} {$i <= 21} {incr i} {
    set eta [expr { -2.958 + $i * 0.0696}]
    add EtaPhiBins $eta $PhiBins
  }

  for {set i 1} {$i <= 21} {incr i} {
    set eta [expr { 1.4964 + $i * 0.0696}]
    add EtaPhiBins $eta $PhiBins
  }

  # take present CMS granularity for HF

  # 0.175 x (0.175 - 0.35) resolution in eta,phi in the HF 3.0 < |eta| < 5.0
  set PhiBins {}
  for {set i -18} {$i <= 18} {incr i} {
    add PhiBins [expr {$i * $pi/18.0}]
  }

  foreach eta {-5 -4.7 -4.525 -4.35 -4.175 -4 -3.825 -3.65 -3.475 -3.3 -3.125 -2.958 3.125 3.3 3.475 3.65 3.825 4 4.175 4.35 4.525 4.7 5} {
    add EtaPhiBins $eta $PhiBins
  }


  # default energy fractions {abs(PDG code)} {Fecal Fhcal}
  add EnergyFraction {0} {1.0}
  # energy fractions for e, gamma and pi0
  add EnergyFraction {11} {0.0}
  add EnergyFraction {22} {0.0}
  add EnergyFraction {111} {0.0}
  # energy fractions for muon, neutrinos and neutralinos
  add EnergyFraction {12} {0.0}
  add EnergyFraction {13} {0.0}
  add EnergyFraction {14} {0.0}
  add EnergyFraction {16} {0.0}
  add EnergyFraction {1000022} {0.0}
  add EnergyFraction {1000023} {0.0}
  add EnergyFraction {1000025} {0.0}
  add EnergyFraction {1000035} {0.0}
  add EnergyFraction {1000045} {0.0}
  # energy fractions for K0short and Lambda
  add EnergyFraction {310} {0.7}
  add EnergyFraction {3122} {0.7}

# set ResolutionFormula {resolution formula as a function of eta and energy}
  set ResolutionFormula {                    (abs(eta) <= 1.5) * sqrt(energy^2*0.05^2 + energy*1.00^2) + \
                                                   (abs(eta) > 1.5 && abs(eta) <= 3.0) * sqrt(energy^2*0.05^2 + energy*1.00^2) + \
                                                   (abs(eta) > 3.0 && abs(eta) <= 5.0) * sqrt(energy^2*0.11^2 + energy*2.80^2)}

}

#################################
# Energy resolution for electrons
#################################

module EnergySmearing PhotonEnergySmearing {
  set InputArray ECal/eflowPhotons
  set OutputArray eflowPhotons

  # adding 1% extra photon smearing
  set ResolutionFormula {energy*0.001}

}



#################
# Electron filter
#################

module PdgCodeFilter ElectronFilter {
  set InputArray HCal/eflowTracks
  set OutputArray electrons
  set Invert true
  add PdgCode {11}
  add PdgCode {-11}
}



##########################
# Track pile-up subtractor
##########################

module TrackPileUpSubtractor TrackPileUpSubtractor {
# add InputArray InputArray OutputArray
  add InputArray HCal/eflowTracks eflowTracks
  add InputArray ElectronFilter/electrons electrons
  add InputArray MuonMomentumSmearing/muons muons

  set VertexInputArray PileUpMerger/vertices
  # assume perfect pile-up subtraction for tracks with |z| > fZVertexResolution
  # Z vertex resolution in m (naive guess from tkLayout, assumed flat in pt for now)
  
  set ZVertexResolution {
  
     (pt < 10. ) * ( ( -0.8*log10(pt) + 1. ) * (10^(0.5*abs(eta) + 1.5)) * 1e-06 ) +
     (pt >= 10. ) * ( ( 0.2 ) * (10^(0.5*abs(eta) + 1.5)) * 1e-06 )
     
     }
 }

########################
# Reco PU filter
########################

module RecoPuFilter RecoPuFilter {
  set InputArray HCal/eflowTracks
  set OutputArray eflowTracks
}

###################################################
# Tower Merger (in case not using e-flow algorithm)
###################################################

module Merger TowerMerger {
# add InputArray InputArray
  add InputArray ECal/ecalTowers
  add InputArray HCal/hcalTowers
  set OutputArray towers
}

####################
# Neutral eflow erger
####################

module Merger NeutralEFlowMerger {
# add InputArray InputArray
  add InputArray PhotonEnergySmearing/eflowPhotons
  add InputArray HCal/eflowNeutralHadrons
  set OutputArray eflowTowers
}

#####################
# Energy flow merger
#####################

module Merger EFlowMerger {
# add InputArray InputArray
  add InputArray HCal/eflowTracks
  add InputArray PhotonEnergySmearing/eflowPhotons
  add InputArray HCal/eflowNeutralHadrons
  set OutputArray eflow
}

############################
# Energy flow merger no PU
############################

module Merger EFlowMergerCHS {
# add InputArray InputArray
  add InputArray RecoPuFilter/eflowTracks
  add InputArray PhotonEnergySmearing/eflowPhotons
  add InputArray HCal/eflowNeutralHadrons
  set OutputArray eflow
}

#########################################
### Run the puppi code (to be tuned) ###
#########################################

module PdgCodeFilter LeptonFilterNoLep {
  set InputArray HCal/eflowTracks
  set OutputArray eflowTracksNoLeptons
  set Invert false
  add PdgCode {13}
  add PdgCode {-13}
  add PdgCode {11}
  add PdgCode {-11}
}

module PdgCodeFilter LeptonFilterLep {
  set InputArray HCal/eflowTracks
  set OutputArray eflowTracksLeptons
  set Invert true
  add PdgCode {11}
  add PdgCode {-11}
  add PdgCode {13}
  add PdgCode {-13}
}

module RunPUPPI RunPUPPIBase {
  ## input information
  set TrackInputArray   LeptonFilterNoLep/eflowTracksNoLeptons
  set NeutralInputArray NeutralEFlowMerger/eflowTowers
  set PVInputArray      PileUpMerger/vertices
  set MinPuppiWeight    0.05
  set UseExp            false
  set UseNoLep          false

  ## define puppi algorithm parameters (more than one for the same eta region is possible)
  add EtaMinBin           0.0   1.5   4.0
  add EtaMaxBin           1.5   4.0   10.0
  add PtMinBin            0.0   0.0   0.0
  add ConeSizeBin         0.2   0.2   0.2
  add RMSPtMinBin         0.1   0.5   0.5
  add RMSScaleFactorBin   1.0   1.0   1.0
  add NeutralMinEBin      0.2   0.2   0.5
  add NeutralPtSlope      0.006 0.013 0.067
  add ApplyCHS            true  true  true
  add UseCharged          true  true  false
  add ApplyLowPUCorr      true  true  true
  add MetricId            5     5     5
  add CombId              0     0     0

  ## output name
  set OutputArray         PuppiParticles
  set OutputArrayTracks   puppiTracks
  set OutputArrayNeutrals puppiNeutrals
}

module Merger RunPUPPIMerger {
  add InputArray RunPUPPIBase/PuppiParticles
  add InputArray LeptonFilterLep/eflowTracksLeptons
  set OutputArray PuppiParticles
}

# need this because of leptons that were added back
module RecoPuFilter RunPUPPI {
  set InputArray RunPUPPIMerger/PuppiParticles
  set OutputArray PuppiParticles
}

######################
# EFlowFilterPuppi
######################

module PdgCodeFilter EFlowFilterPuppi {
  set InputArray RunPUPPI/PuppiParticles
  set OutputArray eflow

  add PdgCode {11}
  add PdgCode {-11}
  add PdgCode {13}
  add PdgCode {-13}
}

######################
# EFlowFilterCHS
######################

module PdgCodeFilter EFlowFilterCHS {
  set InputArray EFlowMergerCHS/eflow
  set OutputArray eflow

  add PdgCode {11}
  add PdgCode {-11}
  add PdgCode {13}
  add PdgCode {-13}
}


###################
# Missing ET merger
###################

module Merger MissingET {
# add InputArray InputArray
#  add InputArray RunPUPPI/PuppiParticles
  add InputArray EFlowMerger/eflow
  set MomentumOutputArray momentum
}

module Merger PuppiMissingET {
  #add InputArray InputArray
  add InputArray RunPUPPI/PuppiParticles
  #add InputArray EFlowMerger/eflow
  set MomentumOutputArray momentum
}

#########################
# Ger PileUp Missing ET
#########################

module Merger GenPileUpMissingET {
# add InputArray InputArray
#  add InputArray RunPUPPI/PuppiParticles
  add InputArray ParticlePropagator/stableParticles
  set MomentumOutputArray momentum
}

##################
# Scalar HT merger
##################

module Merger ScalarHT {
# add InputArray InputArray
  add InputArray RunPUPPI/PuppiParticles
  set EnergyOutputArray energy
}

#################
# Neutrino Filter
#################

module PdgCodeFilter NeutrinoFilter {

  set InputArray Delphes/stableParticles
  set OutputArray filteredParticles

  set PTMin 0.0

  add PdgCode {12}
  add PdgCode {14}
  add PdgCode {16}
  add PdgCode {-12}
  add PdgCode {-14}
  add PdgCode {-16}

}

#####################
# MC truth jet finder
#####################

module FastJetFinder GenJetFinder {
  set InputArray NeutrinoFilter/filteredParticles

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.4

  set JetPTMin 15.0
}

module FastJetFinder GenJetFinderAK8 {
  set InputArray NeutrinoFilter/filteredParticles

  set OutputArray jetsAK8

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.8

  set JetPTMin 200.0
}

#########################
# Gen Missing ET merger
########################

module Merger GenMissingET {

# add InputArray InputArray
  add InputArray NeutrinoFilter/filteredParticles
  set MomentumOutputArray momentum
}


#############
# Rho pile-up
#############

module FastJetGridMedianEstimator Rho {

  set InputArray EFlowMergerCHS/eflow
  set RhoOutputArray rho

  # add GridRange rapmin rapmax drap dphi
  # rapmin - the minimum rapidity extent of the grid
  # rapmax - the maximum rapidity extent of the grid
  # drap - the grid spacing in rapidity
  # dphi - the grid spacing in azimuth

  add GridRange -5.0 -4.0 1.0 1.0
  add GridRange -4.0 -1.5 1.0 1.0
  add GridRange -1.5 1.5 1.0 1.0
  add GridRange 1.5 4.0 1.0 1.0
  add GridRange 4.0 5.0 1.0 1.0

}



module FastJetFinder FastJetFinderPUPPI {
#  set InputArray TowerMerger/towers
  set InputArray RunPUPPI/PuppiParticles

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.4

  set JetPTMin 15.0
}


module FastJetFinder FastJetFinderPUPPIAK8 {
#  set InputArray TowerMerger/towers
  set InputArray RunPUPPI/PuppiParticles

  set OutputArray jets

  set JetAlgorithm 6
  set ParameterR 0.8

  set ComputeNsubjettiness 1
  set Beta 1.0
  set AxisMode 4

  set ComputeTrimming 1
  set RTrim 0.2
  set PtFracTrim 0.05

  set ComputePruning 1
  set ZcutPrun 0.1
  set RcutPrun 0.5
  set RPrun 0.8

  set ComputeSoftDrop 1
  set BetaSoftDrop 0.0
  set SymmetryCutSoftDrop 0.1
  set R0SoftDrop 0.8

  set JetPTMin 200.0
}

##################
# Jet Energy Scale
##################


module EnergyScale JetScalePUPPI {
  set InputArray FastJetFinderPUPPI/jets
  set OutputArray jets

 # scale formula for jets
    ### jetpuppi tightID momentum scale
  set ScaleFormula {
   (abs(eta) > 0.0 && abs(eta) <= 1.3) * (pt > 20.0 && pt <= 50.0) * (0.834) +
   (abs(eta) > 0.0 && abs(eta) <= 1.3) * (pt > 50.0 && pt <= 100.0) * (0.825) +
   (abs(eta) > 0.0 && abs(eta) <= 1.3) * (pt > 100.0 && pt <= 200.0) * (0.859) +
   (abs(eta) > 0.0 && abs(eta) <= 1.3) * (pt > 200.0 && pt <= 400.0) * (0.893) +
   (abs(eta) > 0.0 && abs(eta) <= 1.3) * (pt > 400.0 && pt <= 14000.0) * (0.928) +
   (abs(eta) > 1.3 && abs(eta) <= 2.5) * (pt > 20.0 && pt <= 50.0) * (0.816) +
   (abs(eta) > 1.3 && abs(eta) <= 2.5) * (pt > 50.0 && pt <= 100.0) * (0.690) +
   (abs(eta) > 1.3 && abs(eta) <= 2.5) * (pt > 100.0 && pt <= 200.0) * (0.702) +
   (abs(eta) > 1.3 && abs(eta) <= 2.5) * (pt > 200.0 && pt <= 400.0) * (0.719) +
   (abs(eta) > 1.3 && abs(eta) <= 2.5) * (pt > 400.0 && pt <= 14000.0) * (0.483) +
   (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 50.0) * (0.872) +
   (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 100.0) * (0.843) +
   (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 200.0) * (0.589) +
   (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 200.0 && pt <= 400.0) * (0.714) +
   (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 400.0 && pt <= 14000.0) * (0.678) +
   (abs(eta) > 3.0 && abs(eta) <= 5.0) * (pt > 20.0 && pt <= 50.0) * (0.495) +
   (abs(eta) > 3.0 && abs(eta) <= 5.0) * (pt > 50.0 && pt <= 100.0) * (0.460) +
   (abs(eta) > 3.0 && abs(eta) <= 5.0) * (pt > 100.0 && pt <= 200.0) * (0.533) +
   (abs(eta) > 3.0 && abs(eta) <= 5.0) * (pt > 200.0 && pt <= 400.0) * (1.196) +
   (abs(eta) > 3.0 && abs(eta) <= 5.0) * (pt > 400.0 && pt <= 14000.0) * (0.739)  
  }
}

module EnergyScale JetScalePUPPIAK8 {
  set InputArray FastJetFinderPUPPIAK8/jets
  set OutputArray jets

 # scale formula for jets
  ## DUMMY_JETAK8_SCALE
  set ScaleFormula {1.00}
  ## ENDDUMMY_JETAK8_SCALE
}

####################
# Jet Energy Smear #
####################


module MomentumSmearing JetSmearPUPPI {
  set InputArray JetScalePUPPI/jets
  set OutputArray jets

 # scale formula for jets
   ### jetpuppi tightID momentum resolution
  set ResolutionFormula {
   (abs(eta) > 0.0 && abs(eta) <= 1.3) * (pt > 20.0 && pt <= 50.0) * (0.000001) +
   (abs(eta) > 0.0 && abs(eta) <= 1.3) * (pt > 50.0 && pt <= 100.0) * (0.114070) +
   (abs(eta) > 0.0 && abs(eta) <= 1.3) * (pt > 100.0 && pt <= 200.0) * (0.105469) +
   (abs(eta) > 0.0 && abs(eta) <= 1.3) * (pt > 200.0 && pt <= 400.0) * (0.079494) +
   (abs(eta) > 0.0 && abs(eta) <= 1.3) * (pt > 400.0 && pt <= 14000.0) * (0.057811) +
   (abs(eta) > 1.3 && abs(eta) <= 2.5) * (pt > 20.0 && pt <= 50.0) * (0.174301) +
   (abs(eta) > 1.3 && abs(eta) <= 2.5) * (pt > 50.0 && pt <= 100.0) * (0.170512) +
   (abs(eta) > 1.3 && abs(eta) <= 2.5) * (pt > 100.0 && pt <= 200.0) * (0.174692) +
   (abs(eta) > 1.3 && abs(eta) <= 2.5) * (pt > 200.0 && pt <= 400.0) * (0.261556) +
   (abs(eta) > 1.3 && abs(eta) <= 2.5) * (pt > 400.0 && pt <= 14000.0) * (0.520893) +
   (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 50.0) * (0.145647) +
   (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 100.0) * (0.371395) +
   (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 200.0) * (0.112873) +
   (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 200.0 && pt <= 400.0) * (0.099651) +
   (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 400.0 && pt <= 14000.0) * (0.035230) +
   (abs(eta) > 3.0 && abs(eta) <= 5.0) * (pt > 20.0 && pt <= 50.0) * (0.052101) +
   (abs(eta) > 3.0 && abs(eta) <= 5.0) * (pt > 50.0 && pt <= 100.0) * (0.054670) +
   (abs(eta) > 3.0 && abs(eta) <= 5.0) * (pt > 100.0 && pt <= 200.0) * (0.000001) +
   (abs(eta) > 3.0 && abs(eta) <= 5.0) * (pt > 200.0 && pt <= 400.0) * (0.421392) +
   (abs(eta) > 3.0 && abs(eta) <= 5.0) * (pt > 400.0 && pt <= 14000.0) * (0.256315)  
  }
}

module MomentumSmearing JetSmearPUPPIAK8 {
  set InputArray JetScalePUPPIAK8/jets
  set OutputArray jets

 ## DUMMY_JETPUPPIAK8_SMEAR
  set ResolutionFormula {1.00e-10}
 ## ENDDUMMY_JETPUPPIAK8_SMEAR
}



#####################
# Jet Id Loose     #
#####################

module Efficiency JetLooseID {

  ## input particles
  set InputArray JetSmearPUPPI/jets
  ## output particles
  set OutputArray jets
  # set EfficiencyFormula {efficiency formula as a function of eta and pt}
  # efficiency formula for jets


    ### jetpuppi loose ID 

  set EfficiencyFormula {

   (abs(eta) > 0.0 && abs(eta) <= 0.26) * (pt > 20.0 && pt <= 26.0) * (0.646444691389) +
   (abs(eta) > 0.0 && abs(eta) <= 0.26) * (pt > 26.0 && pt <= 32.0) * (0.744558823529) +
   (abs(eta) > 0.0 && abs(eta) <= 0.26) * (pt > 32.0 && pt <= 38.0) * (0.855654956249) +
   (abs(eta) > 0.0 && abs(eta) <= 0.26) * (pt > 38.0 && pt <= 44.0) * (0.939890710383) +
   (abs(eta) > 0.0 && abs(eta) <= 0.26) * (pt > 44.0 && pt <= 50.0) * (0.89912471853) +
   (abs(eta) > 0.0 && abs(eta) <= 0.26) * (pt > 50.0 && pt <= 60.0) * (0.95388160347) +
   (abs(eta) > 0.0 && abs(eta) <= 0.26) * (pt > 60.0 && pt <= 70.0) * (0.921771413418) +
   (abs(eta) > 0.0 && abs(eta) <= 0.26) * (pt > 70.0 && pt <= 80.0) * (0.969866071429) +
   (abs(eta) > 0.0 && abs(eta) <= 0.26) * (pt > 80.0 && pt <= 90.0) * (0.973088007094) +
   (abs(eta) > 0.0 && abs(eta) <= 0.26) * (pt > 90.0 && pt <= 100.0) * (0.99461557726) +
   (abs(eta) > 0.0 && abs(eta) <= 0.26) * (pt > 100.0 && pt <= 150.0) * (0.985838101588) +
   (abs(eta) > 0.0 && abs(eta) <= 0.26) * (pt > 150.0 && pt <= 200.0) * (0.992843002893) +
   (abs(eta) > 0.0 && abs(eta) <= 0.26) * (pt > 200.0 && pt <= 300.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.26) * (pt > 300.0 && pt <= 400.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.26) * (pt > 400.0) * (0.997006618796) +
   (abs(eta) > 0.26 && abs(eta) <= 0.52) * (pt > 20.0 && pt <= 26.0) * (0.605123595506) +
   (abs(eta) > 0.26 && abs(eta) <= 0.52) * (pt > 26.0 && pt <= 32.0) * (0.715038186483) +
   (abs(eta) > 0.26 && abs(eta) <= 0.52) * (pt > 32.0 && pt <= 38.0) * (0.853598881901) +
   (abs(eta) > 0.26 && abs(eta) <= 0.52) * (pt > 38.0 && pt <= 44.0) * (0.934017329044) +
   (abs(eta) > 0.26 && abs(eta) <= 0.52) * (pt > 44.0 && pt <= 50.0) * (0.898433008809) +
   (abs(eta) > 0.26 && abs(eta) <= 0.52) * (pt > 50.0 && pt <= 60.0) * (0.965772432932) +
   (abs(eta) > 0.26 && abs(eta) <= 0.52) * (pt > 60.0 && pt <= 70.0) * (0.903114186851) +
   (abs(eta) > 0.26 && abs(eta) <= 0.52) * (pt > 70.0 && pt <= 80.0) * (0.960365853659) +
   (abs(eta) > 0.26 && abs(eta) <= 0.52) * (pt > 80.0 && pt <= 90.0) * (0.953780752533) +
   (abs(eta) > 0.26 && abs(eta) <= 0.52) * (pt > 90.0 && pt <= 100.0) * (0.990068580137) +
   (abs(eta) > 0.26 && abs(eta) <= 0.52) * (pt > 100.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 0.26 && abs(eta) <= 0.52) * (pt > 150.0 && pt <= 200.0) * (0.994303162044) +
   (abs(eta) > 0.26 && abs(eta) <= 0.52) * (pt > 200.0 && pt <= 300.0) * (0.999123243408) +
   (abs(eta) > 0.26 && abs(eta) <= 0.52) * (pt > 300.0 && pt <= 400.0) * (1.0) +
   (abs(eta) > 0.26 && abs(eta) <= 0.52) * (pt > 400.0) * (0.999102720867) +
   (abs(eta) > 0.52 && abs(eta) <= 0.78) * (pt > 20.0 && pt <= 26.0) * (0.510210558628) +
   (abs(eta) > 0.52 && abs(eta) <= 0.78) * (pt > 26.0 && pt <= 32.0) * (0.668975535168) +
   (abs(eta) > 0.52 && abs(eta) <= 0.78) * (pt > 32.0 && pt <= 38.0) * (0.847222222222) +
   (abs(eta) > 0.52 && abs(eta) <= 0.78) * (pt > 38.0 && pt <= 44.0) * (0.844137098467) +
   (abs(eta) > 0.52 && abs(eta) <= 0.78) * (pt > 44.0 && pt <= 50.0) * (0.973127380449) +
   (abs(eta) > 0.52 && abs(eta) <= 0.78) * (pt > 50.0 && pt <= 60.0) * (0.949882903981) +
   (abs(eta) > 0.52 && abs(eta) <= 0.78) * (pt > 60.0 && pt <= 70.0) * (0.981816413847) +
   (abs(eta) > 0.52 && abs(eta) <= 0.78) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 0.52 && abs(eta) <= 0.78) * (pt > 80.0 && pt <= 90.0) * (0.995914198161) +
   (abs(eta) > 0.52 && abs(eta) <= 0.78) * (pt > 90.0 && pt <= 100.0) * (0.986799209021) +
   (abs(eta) > 0.52 && abs(eta) <= 0.78) * (pt > 100.0 && pt <= 150.0) * (0.976331882068) +
   (abs(eta) > 0.52 && abs(eta) <= 0.78) * (pt > 150.0 && pt <= 200.0) * (0.983865963388) +
   (abs(eta) > 0.52 && abs(eta) <= 0.78) * (pt > 200.0 && pt <= 300.0) * (0.99041189093) +
   (abs(eta) > 0.52 && abs(eta) <= 0.78) * (pt > 300.0 && pt <= 400.0) * (0.987772780277) +
   (abs(eta) > 0.52 && abs(eta) <= 0.78) * (pt > 400.0) * (0.999882623002) +
   (abs(eta) > 0.78 && abs(eta) <= 1.04) * (pt > 20.0 && pt <= 26.0) * (0.567091087169) +
   (abs(eta) > 0.78 && abs(eta) <= 1.04) * (pt > 26.0 && pt <= 32.0) * (0.905369127517) +
   (abs(eta) > 0.78 && abs(eta) <= 1.04) * (pt > 32.0 && pt <= 38.0) * (0.938603107558) +
   (abs(eta) > 0.78 && abs(eta) <= 1.04) * (pt > 38.0 && pt <= 44.0) * (0.930284640948) +
   (abs(eta) > 0.78 && abs(eta) <= 1.04) * (pt > 44.0 && pt <= 50.0) * (0.9300155521) +
   (abs(eta) > 0.78 && abs(eta) <= 1.04) * (pt > 50.0 && pt <= 60.0) * (0.944572649573) +
   (abs(eta) > 0.78 && abs(eta) <= 1.04) * (pt > 60.0 && pt <= 70.0) * (0.944375259444) +
   (abs(eta) > 0.78 && abs(eta) <= 1.04) * (pt > 70.0 && pt <= 80.0) * (0.956829348182) +
   (abs(eta) > 0.78 && abs(eta) <= 1.04) * (pt > 80.0 && pt <= 90.0) * (0.986268892621) +
   (abs(eta) > 0.78 && abs(eta) <= 1.04) * (pt > 90.0 && pt <= 100.0) * (0.980725120237) +
   (abs(eta) > 0.78 && abs(eta) <= 1.04) * (pt > 100.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 0.78 && abs(eta) <= 1.04) * (pt > 150.0 && pt <= 200.0) * (0.995559453059) +
   (abs(eta) > 0.78 && abs(eta) <= 1.04) * (pt > 200.0 && pt <= 300.0) * (1.0) +
   (abs(eta) > 0.78 && abs(eta) <= 1.04) * (pt > 300.0 && pt <= 400.0) * (0.99445767845) +
   (abs(eta) > 0.78 && abs(eta) <= 1.04) * (pt > 400.0) * (0.99952728788) +
   (abs(eta) > 1.04 && abs(eta) <= 1.3) * (pt > 20.0 && pt <= 26.0) * (0.714522022516) +
   (abs(eta) > 1.04 && abs(eta) <= 1.3) * (pt > 26.0 && pt <= 32.0) * (0.760850559158) +
   (abs(eta) > 1.04 && abs(eta) <= 1.3) * (pt > 32.0 && pt <= 38.0) * (0.852779889864) +
   (abs(eta) > 1.04 && abs(eta) <= 1.3) * (pt > 38.0 && pt <= 44.0) * (0.908245367978) +
   (abs(eta) > 1.04 && abs(eta) <= 1.3) * (pt > 44.0 && pt <= 50.0) * (0.86599967197) +
   (abs(eta) > 1.04 && abs(eta) <= 1.3) * (pt > 50.0 && pt <= 60.0) * (0.969510390149) +
   (abs(eta) > 1.04 && abs(eta) <= 1.3) * (pt > 60.0 && pt <= 70.0) * (0.99462575728) +
   (abs(eta) > 1.04 && abs(eta) <= 1.3) * (pt > 70.0 && pt <= 80.0) * (0.989883573889) +
   (abs(eta) > 1.04 && abs(eta) <= 1.3) * (pt > 80.0 && pt <= 90.0) * (0.982286056254) +
   (abs(eta) > 1.04 && abs(eta) <= 1.3) * (pt > 90.0 && pt <= 100.0) * (0.972712167689) +
   (abs(eta) > 1.04 && abs(eta) <= 1.3) * (pt > 100.0 && pt <= 150.0) * (0.999074478674) +
   (abs(eta) > 1.04 && abs(eta) <= 1.3) * (pt > 150.0 && pt <= 200.0) * (1.0) +
   (abs(eta) > 1.04 && abs(eta) <= 1.3) * (pt > 200.0 && pt <= 300.0) * (1.0) +
   (abs(eta) > 1.04 && abs(eta) <= 1.3) * (pt > 300.0 && pt <= 400.0) * (1.0) +
   (abs(eta) > 1.04 && abs(eta) <= 1.3) * (pt > 400.0) * (1.0) +
   (abs(eta) > 1.3 && abs(eta) <= 1.54) * (pt > 20.0 && pt <= 26.0) * (0.655359565807) +
   (abs(eta) > 1.3 && abs(eta) <= 1.54) * (pt > 26.0 && pt <= 32.0) * (0.856044341235) +
   (abs(eta) > 1.3 && abs(eta) <= 1.54) * (pt > 32.0 && pt <= 38.0) * (0.834927466922) +
   (abs(eta) > 1.3 && abs(eta) <= 1.54) * (pt > 38.0 && pt <= 44.0) * (0.912185341947) +
   (abs(eta) > 1.3 && abs(eta) <= 1.54) * (pt > 44.0 && pt <= 50.0) * (0.884293785311) +
   (abs(eta) > 1.3 && abs(eta) <= 1.54) * (pt > 50.0 && pt <= 60.0) * (0.843465188387) +
   (abs(eta) > 1.3 && abs(eta) <= 1.54) * (pt > 60.0 && pt <= 70.0) * (0.964149504195) +
   (abs(eta) > 1.3 && abs(eta) <= 1.54) * (pt > 70.0 && pt <= 80.0) * (0.887861271676) +
   (abs(eta) > 1.3 && abs(eta) <= 1.54) * (pt > 80.0 && pt <= 90.0) * (0.912464728254) +
   (abs(eta) > 1.3 && abs(eta) <= 1.54) * (pt > 90.0 && pt <= 100.0) * (0.875769285817) +
   (abs(eta) > 1.3 && abs(eta) <= 1.54) * (pt > 100.0 && pt <= 150.0) * (0.907308580627) +
   (abs(eta) > 1.3 && abs(eta) <= 1.54) * (pt > 150.0 && pt <= 200.0) * (0.969284332689) +
   (abs(eta) > 1.3 && abs(eta) <= 1.54) * (pt > 200.0 && pt <= 300.0) * (0.905125228894) +
   (abs(eta) > 1.3 && abs(eta) <= 1.54) * (pt > 300.0 && pt <= 400.0) * (0.935385633896) +
   (abs(eta) > 1.3 && abs(eta) <= 1.54) * (pt > 400.0) * (0.844901394901) +
   (abs(eta) > 1.54 && abs(eta) <= 1.78) * (pt > 20.0 && pt <= 26.0) * (0.627330827068) +
   (abs(eta) > 1.54 && abs(eta) <= 1.78) * (pt > 26.0 && pt <= 32.0) * (0.651911468813) +
   (abs(eta) > 1.54 && abs(eta) <= 1.78) * (pt > 32.0 && pt <= 38.0) * (0.690726817043) +
   (abs(eta) > 1.54 && abs(eta) <= 1.78) * (pt > 38.0 && pt <= 44.0) * (0.752654648722) +
   (abs(eta) > 1.54 && abs(eta) <= 1.78) * (pt > 44.0 && pt <= 50.0) * (0.733203505355) +
   (abs(eta) > 1.54 && abs(eta) <= 1.78) * (pt > 50.0 && pt <= 60.0) * (0.824566943342) +
   (abs(eta) > 1.54 && abs(eta) <= 1.78) * (pt > 60.0 && pt <= 70.0) * (0.829189589505) +
   (abs(eta) > 1.54 && abs(eta) <= 1.78) * (pt > 70.0 && pt <= 80.0) * (0.71083172147) +
   (abs(eta) > 1.54 && abs(eta) <= 1.78) * (pt > 80.0 && pt <= 90.0) * (0.794794532512) +
   (abs(eta) > 1.54 && abs(eta) <= 1.78) * (pt > 90.0 && pt <= 100.0) * (0.787025485653) +
   (abs(eta) > 1.54 && abs(eta) <= 1.78) * (pt > 100.0 && pt <= 150.0) * (0.775512180009) +
   (abs(eta) > 1.54 && abs(eta) <= 1.78) * (pt > 150.0 && pt <= 200.0) * (0.786866177463) +
   (abs(eta) > 1.54 && abs(eta) <= 1.78) * (pt > 200.0 && pt <= 300.0) * (0.73989316329) +
   (abs(eta) > 1.54 && abs(eta) <= 1.78) * (pt > 300.0 && pt <= 400.0) * (0.623988387857) +
   (abs(eta) > 1.54 && abs(eta) <= 1.78) * (pt > 400.0) * (0.596799473371) +
   (abs(eta) > 1.78 && abs(eta) <= 2.02) * (pt > 20.0 && pt <= 26.0) * (0.417905288628) +
   (abs(eta) > 1.78 && abs(eta) <= 2.02) * (pt > 26.0 && pt <= 32.0) * (0.445069953364) +
   (abs(eta) > 1.78 && abs(eta) <= 2.02) * (pt > 32.0 && pt <= 38.0) * (0.784336645237) +
   (abs(eta) > 1.78 && abs(eta) <= 2.02) * (pt > 38.0 && pt <= 44.0) * (0.628650954309) +
   (abs(eta) > 1.78 && abs(eta) <= 2.02) * (pt > 44.0 && pt <= 50.0) * (0.521505376344) +
   (abs(eta) > 1.78 && abs(eta) <= 2.02) * (pt > 50.0 && pt <= 60.0) * (0.726041666667) +
   (abs(eta) > 1.78 && abs(eta) <= 2.02) * (pt > 60.0 && pt <= 70.0) * (0.821764911057) +
   (abs(eta) > 1.78 && abs(eta) <= 2.02) * (pt > 70.0 && pt <= 80.0) * (0.687034793145) +
   (abs(eta) > 1.78 && abs(eta) <= 2.02) * (pt > 80.0 && pt <= 90.0) * (0.670454545455) +
   (abs(eta) > 1.78 && abs(eta) <= 2.02) * (pt > 90.0 && pt <= 100.0) * (0.636713735558) +
   (abs(eta) > 1.78 && abs(eta) <= 2.02) * (pt > 100.0 && pt <= 150.0) * (0.723881492143) +
   (abs(eta) > 1.78 && abs(eta) <= 2.02) * (pt > 150.0 && pt <= 200.0) * (0.775281276238) +
   (abs(eta) > 1.78 && abs(eta) <= 2.02) * (pt > 200.0 && pt <= 300.0) * (0.655152135561) +
   (abs(eta) > 1.78 && abs(eta) <= 2.02) * (pt > 300.0 && pt <= 400.0) * (0.608727528667) +
   (abs(eta) > 1.78 && abs(eta) <= 2.02) * (pt > 400.0) * (0.578990434807) +
   (abs(eta) > 2.02 && abs(eta) <= 2.26) * (pt > 20.0 && pt <= 26.0) * (0.366119149064) +
   (abs(eta) > 2.02 && abs(eta) <= 2.26) * (pt > 26.0 && pt <= 32.0) * (0.447434113151) +
   (abs(eta) > 2.02 && abs(eta) <= 2.26) * (pt > 32.0 && pt <= 38.0) * (0.566811091854) +
   (abs(eta) > 2.02 && abs(eta) <= 2.26) * (pt > 38.0 && pt <= 44.0) * (0.631193643425) +
   (abs(eta) > 2.02 && abs(eta) <= 2.26) * (pt > 44.0 && pt <= 50.0) * (0.653942428035) +
   (abs(eta) > 2.02 && abs(eta) <= 2.26) * (pt > 50.0 && pt <= 60.0) * (0.758013329102) +
   (abs(eta) > 2.02 && abs(eta) <= 2.26) * (pt > 60.0 && pt <= 70.0) * (0.710869206663) +
   (abs(eta) > 2.02 && abs(eta) <= 2.26) * (pt > 70.0 && pt <= 80.0) * (0.808382363623) +
   (abs(eta) > 2.02 && abs(eta) <= 2.26) * (pt > 80.0 && pt <= 90.0) * (0.829787234043) +
   (abs(eta) > 2.02 && abs(eta) <= 2.26) * (pt > 90.0 && pt <= 100.0) * (0.683744465528) +
   (abs(eta) > 2.02 && abs(eta) <= 2.26) * (pt > 100.0 && pt <= 150.0) * (0.680985221675) +
   (abs(eta) > 2.02 && abs(eta) <= 2.26) * (pt > 150.0 && pt <= 200.0) * (0.614848793748) +
   (abs(eta) > 2.02 && abs(eta) <= 2.26) * (pt > 200.0 && pt <= 300.0) * (0.637782214714) +
   (abs(eta) > 2.02 && abs(eta) <= 2.26) * (pt > 300.0 && pt <= 400.0) * (0.584644194757) +
   (abs(eta) > 2.02 && abs(eta) <= 2.26) * (pt > 400.0) * (0.544117647059) +
   (abs(eta) > 2.26 && abs(eta) <= 2.5) * (pt > 20.0 && pt <= 26.0) * (0.445076586433) +
   (abs(eta) > 2.26 && abs(eta) <= 2.5) * (pt > 26.0 && pt <= 32.0) * (0.478021978022) +
   (abs(eta) > 2.26 && abs(eta) <= 2.5) * (pt > 32.0 && pt <= 38.0) * (0.696459016393) +
   (abs(eta) > 2.26 && abs(eta) <= 2.5) * (pt > 38.0 && pt <= 44.0) * (0.746429381093) +
   (abs(eta) > 2.26 && abs(eta) <= 2.5) * (pt > 44.0 && pt <= 50.0) * (0.716761939886) +
   (abs(eta) > 2.26 && abs(eta) <= 2.5) * (pt > 50.0 && pt <= 60.0) * (0.743937563876) +
   (abs(eta) > 2.26 && abs(eta) <= 2.5) * (pt > 60.0 && pt <= 70.0) * (0.668055555556) +
   (abs(eta) > 2.26 && abs(eta) <= 2.5) * (pt > 70.0 && pt <= 80.0) * (0.735176619008) +
   (abs(eta) > 2.26 && abs(eta) <= 2.5) * (pt > 80.0 && pt <= 90.0) * (0.67365967366) +
   (abs(eta) > 2.26 && abs(eta) <= 2.5) * (pt > 90.0 && pt <= 100.0) * (0.760406424123) +
   (abs(eta) > 2.26 && abs(eta) <= 2.5) * (pt > 100.0 && pt <= 150.0) * (0.749638242894) +
   (abs(eta) > 2.26 && abs(eta) <= 2.5) * (pt > 150.0 && pt <= 200.0) * (0.622641509434) +
   (abs(eta) > 2.26 && abs(eta) <= 2.5) * (pt > 200.0 && pt <= 300.0) * (0.670033670034) +
   (abs(eta) > 2.26 && abs(eta) <= 2.5) * (pt > 300.0 && pt <= 400.0) * (0.686488169365) +
   (abs(eta) > 2.26 && abs(eta) <= 2.5) * (pt > 400.0) * (0.583333333333) +
   (abs(eta) > 2.5 && abs(eta) <= 2.6) * (pt > 20.0 && pt <= 26.0) * (0.883458646617) +
   (abs(eta) > 2.5 && abs(eta) <= 2.6) * (pt > 26.0 && pt <= 32.0) * (0.578534031414) +
   (abs(eta) > 2.5 && abs(eta) <= 2.6) * (pt > 32.0 && pt <= 38.0) * (0.60960960961) +
   (abs(eta) > 2.5 && abs(eta) <= 2.6) * (pt > 38.0 && pt <= 44.0) * (0.87910501624) +
   (abs(eta) > 2.5 && abs(eta) <= 2.6) * (pt > 44.0 && pt <= 50.0) * (0.797720797721) +
   (abs(eta) > 2.5 && abs(eta) <= 2.6) * (pt > 50.0 && pt <= 60.0) * (0.719476744186) +
   (abs(eta) > 2.5 && abs(eta) <= 2.6) * (pt > 60.0 && pt <= 70.0) * (0.799515445185) +
   (abs(eta) > 2.5 && abs(eta) <= 2.6) * (pt > 70.0 && pt <= 80.0) * (0.611650485437) +
   (abs(eta) > 2.5 && abs(eta) <= 2.6) * (pt > 80.0 && pt <= 90.0) * (0.823529411765) +
   (abs(eta) > 2.5 && abs(eta) <= 2.6) * (pt > 90.0 && pt <= 100.0) * (0.336879432624) +
   (abs(eta) > 2.5 && abs(eta) <= 2.6) * (pt > 100.0 && pt <= 150.0) * (0.766780025058) +
   (abs(eta) > 2.5 && abs(eta) <= 2.6) * (pt > 150.0 && pt <= 200.0) * (0.653846153846) +
   (abs(eta) > 2.5 && abs(eta) <= 2.6) * (pt > 200.0 && pt <= 300.0) * (0.673611111111) +
   (abs(eta) > 2.5 && abs(eta) <= 2.6) * (pt > 300.0 && pt <= 400.0) * (0.5) +
   (abs(eta) > 2.5 && abs(eta) <= 2.6) * (pt > 400.0) * (0.333333333333) +
   (abs(eta) > 2.6 && abs(eta) <= 2.7) * (pt > 20.0 && pt <= 26.0) * (0.586092715232) +
   (abs(eta) > 2.6 && abs(eta) <= 2.7) * (pt > 26.0 && pt <= 32.0) * (0.748412040873) +
   (abs(eta) > 2.6 && abs(eta) <= 2.7) * (pt > 32.0 && pt <= 38.0) * (0.812154696133) +
   (abs(eta) > 2.6 && abs(eta) <= 2.7) * (pt > 38.0 && pt <= 44.0) * (0.644016227181) +
   (abs(eta) > 2.6 && abs(eta) <= 2.7) * (pt > 44.0 && pt <= 50.0) * (0.971715328467) +
   (abs(eta) > 2.6 && abs(eta) <= 2.7) * (pt > 50.0 && pt <= 60.0) * (0.617045454545) +
   (abs(eta) > 2.6 && abs(eta) <= 2.7) * (pt > 60.0 && pt <= 70.0) * (0.821829163072) +
   (abs(eta) > 2.6 && abs(eta) <= 2.7) * (pt > 70.0 && pt <= 80.0) * (0.681818181818) +
   (abs(eta) > 2.6 && abs(eta) <= 2.7) * (pt > 80.0 && pt <= 90.0) * (0.815384615385) +
   (abs(eta) > 2.6 && abs(eta) <= 2.7) * (pt > 90.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 2.6 && abs(eta) <= 2.7) * (pt > 100.0 && pt <= 150.0) * (0.720190779014) +
   (abs(eta) > 2.6 && abs(eta) <= 2.7) * (pt > 150.0 && pt <= 200.0) * (0.592079969243) +
   (abs(eta) > 2.6 && abs(eta) <= 2.7) * (pt > 200.0 && pt <= 300.0) * (0.5) +
   (abs(eta) > 2.6 && abs(eta) <= 2.7) * (pt > 300.0 && pt <= 400.0) * (0.338235294118) +
   (abs(eta) > 2.6 && abs(eta) <= 2.7) * (pt > 400.0) * (1.0) +
   (abs(eta) > 2.7 && abs(eta) <= 2.8) * (pt > 20.0 && pt <= 26.0) * (0.605580843064) +
   (abs(eta) > 2.7 && abs(eta) <= 2.8) * (pt > 26.0 && pt <= 32.0) * (0.660406403941) +
   (abs(eta) > 2.7 && abs(eta) <= 2.8) * (pt > 32.0 && pt <= 38.0) * (0.862090290662) +
   (abs(eta) > 2.7 && abs(eta) <= 2.8) * (pt > 38.0 && pt <= 44.0) * (0.928034371643) +
   (abs(eta) > 2.7 && abs(eta) <= 2.8) * (pt > 44.0 && pt <= 50.0) * (0.825958702065) +
   (abs(eta) > 2.7 && abs(eta) <= 2.8) * (pt > 50.0 && pt <= 60.0) * (0.798214285714) +
   (abs(eta) > 2.7 && abs(eta) <= 2.8) * (pt > 60.0 && pt <= 70.0) * (0.801346801347) +
   (abs(eta) > 2.7 && abs(eta) <= 2.8) * (pt > 70.0 && pt <= 80.0) * (0.782726045884) +
   (abs(eta) > 2.7 && abs(eta) <= 2.8) * (pt > 80.0 && pt <= 90.0) * (0.72188449848) +
   (abs(eta) > 2.7 && abs(eta) <= 2.8) * (pt > 90.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 2.7 && abs(eta) <= 2.8) * (pt > 100.0 && pt <= 150.0) * (0.757411114288) +
   (abs(eta) > 2.7 && abs(eta) <= 2.8) * (pt > 150.0 && pt <= 200.0) * (0.555555555556) +
   (abs(eta) > 2.7 && abs(eta) <= 2.8) * (pt > 200.0 && pt <= 300.0) * (0.696985446985) +
   (abs(eta) > 2.7 && abs(eta) <= 2.8) * (pt > 300.0 && pt <= 400.0) * (0.5) +
   (abs(eta) > 2.7 && abs(eta) <= 2.8) * (pt > 400.0) * (0.0) +
   (abs(eta) > 2.8 && abs(eta) <= 2.9) * (pt > 20.0 && pt <= 26.0) * (0.198859166011) +
   (abs(eta) > 2.8 && abs(eta) <= 2.9) * (pt > 26.0 && pt <= 32.0) * (0.647429519071) +
   (abs(eta) > 2.8 && abs(eta) <= 2.9) * (pt > 32.0 && pt <= 38.0) * (0.847368421053) +
   (abs(eta) > 2.8 && abs(eta) <= 2.9) * (pt > 38.0 && pt <= 44.0) * (0.740236686391) +
   (abs(eta) > 2.8 && abs(eta) <= 2.9) * (pt > 44.0 && pt <= 50.0) * (0.738362760835) +
   (abs(eta) > 2.8 && abs(eta) <= 2.9) * (pt > 50.0 && pt <= 60.0) * (0.857843137255) +
   (abs(eta) > 2.8 && abs(eta) <= 2.9) * (pt > 60.0 && pt <= 70.0) * (0.727131782946) +
   (abs(eta) > 2.8 && abs(eta) <= 2.9) * (pt > 70.0 && pt <= 80.0) * (0.633802816901) +
   (abs(eta) > 2.8 && abs(eta) <= 2.9) * (pt > 80.0 && pt <= 90.0) * (0.714285714286) +
   (abs(eta) > 2.8 && abs(eta) <= 2.9) * (pt > 90.0 && pt <= 100.0) * (0.810256410256) +
   (abs(eta) > 2.8 && abs(eta) <= 2.9) * (pt > 100.0 && pt <= 150.0) * (0.540887040887) +
   (abs(eta) > 2.8 && abs(eta) <= 2.9) * (pt > 150.0 && pt <= 200.0) * (0.666666666667) +
   (abs(eta) > 2.8 && abs(eta) <= 2.9) * (pt > 200.0 && pt <= 300.0) * (0.428571428571) +
   (abs(eta) > 2.8 && abs(eta) <= 2.9) * (pt > 300.0 && pt <= 400.0) * (0.5) +
   (abs(eta) > 2.8 && abs(eta) <= 2.9) * (pt > 400.0) * (0.0) +
   (abs(eta) > 2.9 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 26.0) * (0.398034398034) +
   (abs(eta) > 2.9 && abs(eta) <= 3.0) * (pt > 26.0 && pt <= 32.0) * (0.515294117647) +
   (abs(eta) > 2.9 && abs(eta) <= 3.0) * (pt > 32.0 && pt <= 38.0) * (0.493665540541) +
   (abs(eta) > 2.9 && abs(eta) <= 3.0) * (pt > 38.0 && pt <= 44.0) * (0.643570170288) +
   (abs(eta) > 2.9 && abs(eta) <= 3.0) * (pt > 44.0 && pt <= 50.0) * (0.727272727273) +
   (abs(eta) > 2.9 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 60.0) * (0.514492753623) +
   (abs(eta) > 2.9 && abs(eta) <= 3.0) * (pt > 60.0 && pt <= 70.0) * (0.433613445378) +
   (abs(eta) > 2.9 && abs(eta) <= 3.0) * (pt > 70.0 && pt <= 80.0) * (0.477453580902) +
   (abs(eta) > 2.9 && abs(eta) <= 3.0) * (pt > 80.0 && pt <= 90.0) * (0.926470588235) +
   (abs(eta) > 2.9 && abs(eta) <= 3.0) * (pt > 90.0 && pt <= 100.0) * (0.445054945055) +
   (abs(eta) > 2.9 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 150.0) * (0.737662337662) +
   (abs(eta) > 2.9 && abs(eta) <= 3.0) * (pt > 150.0 && pt <= 200.0) * (0.691588785047) +
   (abs(eta) > 2.9 && abs(eta) <= 3.0) * (pt > 200.0 && pt <= 300.0) * (0.333333333333) +
   (abs(eta) > 2.9 && abs(eta) <= 3.0) * (pt > 300.0 && pt <= 400.0) * (0.666666666667) +
   (abs(eta) > 2.9 && abs(eta) <= 3.0) * (pt > 400.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 20.0 && pt <= 26.0) * (0.0995847464078) +
   (abs(eta) > 3.0) * (pt > 26.0 && pt <= 32.0) * (0.142522011121) +
   (abs(eta) > 3.0) * (pt > 32.0 && pt <= 38.0) * (0.204736303123) +
   (abs(eta) > 3.0) * (pt > 38.0 && pt <= 44.0) * (0.355622590917) +
   (abs(eta) > 3.0) * (pt > 44.0 && pt <= 50.0) * (0.306850603035) +
   (abs(eta) > 3.0) * (pt > 50.0 && pt <= 60.0) * (0.386136916801) +
   (abs(eta) > 3.0) * (pt > 60.0 && pt <= 70.0) * (0.394656840166) +
   (abs(eta) > 3.0) * (pt > 70.0 && pt <= 80.0) * (0.374466950959) +
   (abs(eta) > 3.0) * (pt > 80.0 && pt <= 90.0) * (0.318655216517) +
   (abs(eta) > 3.0) * (pt > 90.0 && pt <= 100.0) * (0.286337843736) +
   (abs(eta) > 3.0) * (pt > 100.0 && pt <= 150.0) * (0.253841890462) +
   (abs(eta) > 3.0) * (pt > 150.0 && pt <= 200.0) * (0.277247956403) +
   (abs(eta) > 3.0) * (pt > 200.0 && pt <= 300.0) * (0.347826086957) +
   (abs(eta) > 3.0) * (pt > 300.0 && pt <= 400.0) * (0.5) +
   (abs(eta) > 3.0) * (pt > 400.0) * (1.0)  
  }

}



#####################
# Jet Id Tight      #
#####################

module Efficiency JetTightID {

  ## input particles
  set InputArray JetSmearPUPPI/jets
  ## output particles
  set OutputArray jets
  # set EfficiencyFormula {efficiency formula as a function of eta and pt}
  # efficiency formula for jets


    ### jetpuppi tight ID 

  set EfficiencyFormula {

   (abs(eta) > 0.0 && abs(eta) <= 0.26) * (pt > 20.0 && pt <= 26.0) * (0.646444691389) +
   (abs(eta) > 0.0 && abs(eta) <= 0.26) * (pt > 26.0 && pt <= 32.0) * (0.744558823529) +
   (abs(eta) > 0.0 && abs(eta) <= 0.26) * (pt > 32.0 && pt <= 38.0) * (0.855654956249) +
   (abs(eta) > 0.0 && abs(eta) <= 0.26) * (pt > 38.0 && pt <= 44.0) * (0.939890710383) +
   (abs(eta) > 0.0 && abs(eta) <= 0.26) * (pt > 44.0 && pt <= 50.0) * (0.89912471853) +
   (abs(eta) > 0.0 && abs(eta) <= 0.26) * (pt > 50.0 && pt <= 60.0) * (0.95388160347) +
   (abs(eta) > 0.0 && abs(eta) <= 0.26) * (pt > 60.0 && pt <= 70.0) * (0.921771413418) +
   (abs(eta) > 0.0 && abs(eta) <= 0.26) * (pt > 70.0 && pt <= 80.0) * (0.969866071429) +
   (abs(eta) > 0.0 && abs(eta) <= 0.26) * (pt > 80.0 && pt <= 90.0) * (0.973088007094) +
   (abs(eta) > 0.0 && abs(eta) <= 0.26) * (pt > 90.0 && pt <= 100.0) * (0.99461557726) +
   (abs(eta) > 0.0 && abs(eta) <= 0.26) * (pt > 100.0 && pt <= 150.0) * (0.985838101588) +
   (abs(eta) > 0.0 && abs(eta) <= 0.26) * (pt > 150.0 && pt <= 200.0) * (0.992843002893) +
   (abs(eta) > 0.0 && abs(eta) <= 0.26) * (pt > 200.0 && pt <= 300.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.26) * (pt > 300.0 && pt <= 400.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.26) * (pt > 400.0) * (0.991932793764) +
   (abs(eta) > 0.26 && abs(eta) <= 0.52) * (pt > 20.0 && pt <= 26.0) * (0.605123595506) +
   (abs(eta) > 0.26 && abs(eta) <= 0.52) * (pt > 26.0 && pt <= 32.0) * (0.715038186483) +
   (abs(eta) > 0.26 && abs(eta) <= 0.52) * (pt > 32.0 && pt <= 38.0) * (0.84343699045) +
   (abs(eta) > 0.26 && abs(eta) <= 0.52) * (pt > 38.0 && pt <= 44.0) * (0.934017329044) +
   (abs(eta) > 0.26 && abs(eta) <= 0.52) * (pt > 44.0 && pt <= 50.0) * (0.898433008809) +
   (abs(eta) > 0.26 && abs(eta) <= 0.52) * (pt > 50.0 && pt <= 60.0) * (0.965772432932) +
   (abs(eta) > 0.26 && abs(eta) <= 0.52) * (pt > 60.0 && pt <= 70.0) * (0.903114186851) +
   (abs(eta) > 0.26 && abs(eta) <= 0.52) * (pt > 70.0 && pt <= 80.0) * (0.960365853659) +
   (abs(eta) > 0.26 && abs(eta) <= 0.52) * (pt > 80.0 && pt <= 90.0) * (0.953780752533) +
   (abs(eta) > 0.26 && abs(eta) <= 0.52) * (pt > 90.0 && pt <= 100.0) * (0.990068580137) +
   (abs(eta) > 0.26 && abs(eta) <= 0.52) * (pt > 100.0 && pt <= 150.0) * (0.995895398875) +
   (abs(eta) > 0.26 && abs(eta) <= 0.52) * (pt > 150.0 && pt <= 200.0) * (0.994303162044) +
   (abs(eta) > 0.26 && abs(eta) <= 0.52) * (pt > 200.0 && pt <= 300.0) * (0.993946439038) +
   (abs(eta) > 0.26 && abs(eta) <= 0.52) * (pt > 300.0 && pt <= 400.0) * (1.0) +
   (abs(eta) > 0.26 && abs(eta) <= 0.52) * (pt > 400.0) * (0.996507648865) +
   (abs(eta) > 0.52 && abs(eta) <= 0.78) * (pt > 20.0 && pt <= 26.0) * (0.510210558628) +
   (abs(eta) > 0.52 && abs(eta) <= 0.78) * (pt > 26.0 && pt <= 32.0) * (0.668975535168) +
   (abs(eta) > 0.52 && abs(eta) <= 0.78) * (pt > 32.0 && pt <= 38.0) * (0.847222222222) +
   (abs(eta) > 0.52 && abs(eta) <= 0.78) * (pt > 38.0 && pt <= 44.0) * (0.844137098467) +
   (abs(eta) > 0.52 && abs(eta) <= 0.78) * (pt > 44.0 && pt <= 50.0) * (0.973127380449) +
   (abs(eta) > 0.52 && abs(eta) <= 0.78) * (pt > 50.0 && pt <= 60.0) * (0.949882903981) +
   (abs(eta) > 0.52 && abs(eta) <= 0.78) * (pt > 60.0 && pt <= 70.0) * (0.981816413847) +
   (abs(eta) > 0.52 && abs(eta) <= 0.78) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 0.52 && abs(eta) <= 0.78) * (pt > 80.0 && pt <= 90.0) * (0.995914198161) +
   (abs(eta) > 0.52 && abs(eta) <= 0.78) * (pt > 90.0 && pt <= 100.0) * (0.986799209021) +
   (abs(eta) > 0.52 && abs(eta) <= 0.78) * (pt > 100.0 && pt <= 150.0) * (0.976331882068) +
   (abs(eta) > 0.52 && abs(eta) <= 0.78) * (pt > 150.0 && pt <= 200.0) * (0.983865963388) +
   (abs(eta) > 0.52 && abs(eta) <= 0.78) * (pt > 200.0 && pt <= 300.0) * (0.99041189093) +
   (abs(eta) > 0.52 && abs(eta) <= 0.78) * (pt > 300.0 && pt <= 400.0) * (0.987772780277) +
   (abs(eta) > 0.52 && abs(eta) <= 0.78) * (pt > 400.0) * (0.999882623002) +
   (abs(eta) > 0.78 && abs(eta) <= 1.04) * (pt > 20.0 && pt <= 26.0) * (0.567091087169) +
   (abs(eta) > 0.78 && abs(eta) <= 1.04) * (pt > 26.0 && pt <= 32.0) * (0.905369127517) +
   (abs(eta) > 0.78 && abs(eta) <= 1.04) * (pt > 32.0 && pt <= 38.0) * (0.938603107558) +
   (abs(eta) > 0.78 && abs(eta) <= 1.04) * (pt > 38.0 && pt <= 44.0) * (0.918508885999) +
   (abs(eta) > 0.78 && abs(eta) <= 1.04) * (pt > 44.0 && pt <= 50.0) * (0.9300155521) +
   (abs(eta) > 0.78 && abs(eta) <= 1.04) * (pt > 50.0 && pt <= 60.0) * (0.944572649573) +
   (abs(eta) > 0.78 && abs(eta) <= 1.04) * (pt > 60.0 && pt <= 70.0) * (0.944375259444) +
   (abs(eta) > 0.78 && abs(eta) <= 1.04) * (pt > 70.0 && pt <= 80.0) * (0.956829348182) +
   (abs(eta) > 0.78 && abs(eta) <= 1.04) * (pt > 80.0 && pt <= 90.0) * (0.986268892621) +
   (abs(eta) > 0.78 && abs(eta) <= 1.04) * (pt > 90.0 && pt <= 100.0) * (0.980725120237) +
   (abs(eta) > 0.78 && abs(eta) <= 1.04) * (pt > 100.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 0.78 && abs(eta) <= 1.04) * (pt > 150.0 && pt <= 200.0) * (0.995559453059) +
   (abs(eta) > 0.78 && abs(eta) <= 1.04) * (pt > 200.0 && pt <= 300.0) * (1.0) +
   (abs(eta) > 0.78 && abs(eta) <= 1.04) * (pt > 300.0 && pt <= 400.0) * (0.99445767845) +
   (abs(eta) > 0.78 && abs(eta) <= 1.04) * (pt > 400.0) * (0.99952728788) +
   (abs(eta) > 1.04 && abs(eta) <= 1.3) * (pt > 20.0 && pt <= 26.0) * (0.714522022516) +
   (abs(eta) > 1.04 && abs(eta) <= 1.3) * (pt > 26.0 && pt <= 32.0) * (0.760850559158) +
   (abs(eta) > 1.04 && abs(eta) <= 1.3) * (pt > 32.0 && pt <= 38.0) * (0.852779889864) +
   (abs(eta) > 1.04 && abs(eta) <= 1.3) * (pt > 38.0 && pt <= 44.0) * (0.908245367978) +
   (abs(eta) > 1.04 && abs(eta) <= 1.3) * (pt > 44.0 && pt <= 50.0) * (0.86599967197) +
   (abs(eta) > 1.04 && abs(eta) <= 1.3) * (pt > 50.0 && pt <= 60.0) * (0.969510390149) +
   (abs(eta) > 1.04 && abs(eta) <= 1.3) * (pt > 60.0 && pt <= 70.0) * (0.99462575728) +
   (abs(eta) > 1.04 && abs(eta) <= 1.3) * (pt > 70.0 && pt <= 80.0) * (0.989883573889) +
   (abs(eta) > 1.04 && abs(eta) <= 1.3) * (pt > 80.0 && pt <= 90.0) * (0.982286056254) +
   (abs(eta) > 1.04 && abs(eta) <= 1.3) * (pt > 90.0 && pt <= 100.0) * (0.972712167689) +
   (abs(eta) > 1.04 && abs(eta) <= 1.3) * (pt > 100.0 && pt <= 150.0) * (0.999074478674) +
   (abs(eta) > 1.04 && abs(eta) <= 1.3) * (pt > 150.0 && pt <= 200.0) * (1.0) +
   (abs(eta) > 1.04 && abs(eta) <= 1.3) * (pt > 200.0 && pt <= 300.0) * (0.991326048506) +
   (abs(eta) > 1.04 && abs(eta) <= 1.3) * (pt > 300.0 && pt <= 400.0) * (1.0) +
   (abs(eta) > 1.04 && abs(eta) <= 1.3) * (pt > 400.0) * (1.0) +
   (abs(eta) > 1.3 && abs(eta) <= 1.54) * (pt > 20.0 && pt <= 26.0) * (0.655359565807) +
   (abs(eta) > 1.3 && abs(eta) <= 1.54) * (pt > 26.0 && pt <= 32.0) * (0.856044341235) +
   (abs(eta) > 1.3 && abs(eta) <= 1.54) * (pt > 32.0 && pt <= 38.0) * (0.818228917583) +
   (abs(eta) > 1.3 && abs(eta) <= 1.54) * (pt > 38.0 && pt <= 44.0) * (0.912185341947) +
   (abs(eta) > 1.3 && abs(eta) <= 1.54) * (pt > 44.0 && pt <= 50.0) * (0.884293785311) +
   (abs(eta) > 1.3 && abs(eta) <= 1.54) * (pt > 50.0 && pt <= 60.0) * (0.843465188387) +
   (abs(eta) > 1.3 && abs(eta) <= 1.54) * (pt > 60.0 && pt <= 70.0) * (0.964149504195) +
   (abs(eta) > 1.3 && abs(eta) <= 1.54) * (pt > 70.0 && pt <= 80.0) * (0.887861271676) +
   (abs(eta) > 1.3 && abs(eta) <= 1.54) * (pt > 80.0 && pt <= 90.0) * (0.912464728254) +
   (abs(eta) > 1.3 && abs(eta) <= 1.54) * (pt > 90.0 && pt <= 100.0) * (0.845570344926) +
   (abs(eta) > 1.3 && abs(eta) <= 1.54) * (pt > 100.0 && pt <= 150.0) * (0.900685890258) +
   (abs(eta) > 1.3 && abs(eta) <= 1.54) * (pt > 150.0 && pt <= 200.0) * (0.969284332689) +
   (abs(eta) > 1.3 && abs(eta) <= 1.54) * (pt > 200.0 && pt <= 300.0) * (0.905125228894) +
   (abs(eta) > 1.3 && abs(eta) <= 1.54) * (pt > 300.0 && pt <= 400.0) * (0.935385633896) +
   (abs(eta) > 1.3 && abs(eta) <= 1.54) * (pt > 400.0) * (0.844901394901) +
   (abs(eta) > 1.54 && abs(eta) <= 1.78) * (pt > 20.0 && pt <= 26.0) * (0.627330827068) +
   (abs(eta) > 1.54 && abs(eta) <= 1.78) * (pt > 26.0 && pt <= 32.0) * (0.651911468813) +
   (abs(eta) > 1.54 && abs(eta) <= 1.78) * (pt > 32.0 && pt <= 38.0) * (0.690726817043) +
   (abs(eta) > 1.54 && abs(eta) <= 1.78) * (pt > 38.0 && pt <= 44.0) * (0.752654648722) +
   (abs(eta) > 1.54 && abs(eta) <= 1.78) * (pt > 44.0 && pt <= 50.0) * (0.733203505355) +
   (abs(eta) > 1.54 && abs(eta) <= 1.78) * (pt > 50.0 && pt <= 60.0) * (0.824566943342) +
   (abs(eta) > 1.54 && abs(eta) <= 1.78) * (pt > 60.0 && pt <= 70.0) * (0.829189589505) +
   (abs(eta) > 1.54 && abs(eta) <= 1.78) * (pt > 70.0 && pt <= 80.0) * (0.71083172147) +
   (abs(eta) > 1.54 && abs(eta) <= 1.78) * (pt > 80.0 && pt <= 90.0) * (0.794794532512) +
   (abs(eta) > 1.54 && abs(eta) <= 1.78) * (pt > 90.0 && pt <= 100.0) * (0.787025485653) +
   (abs(eta) > 1.54 && abs(eta) <= 1.78) * (pt > 100.0 && pt <= 150.0) * (0.775512180009) +
   (abs(eta) > 1.54 && abs(eta) <= 1.78) * (pt > 150.0 && pt <= 200.0) * (0.786866177463) +
   (abs(eta) > 1.54 && abs(eta) <= 1.78) * (pt > 200.0 && pt <= 300.0) * (0.73989316329) +
   (abs(eta) > 1.54 && abs(eta) <= 1.78) * (pt > 300.0 && pt <= 400.0) * (0.623988387857) +
   (abs(eta) > 1.54 && abs(eta) <= 1.78) * (pt > 400.0) * (0.596799473371) +
   (abs(eta) > 1.78 && abs(eta) <= 2.02) * (pt > 20.0 && pt <= 26.0) * (0.417905288628) +
   (abs(eta) > 1.78 && abs(eta) <= 2.02) * (pt > 26.0 && pt <= 32.0) * (0.445069953364) +
   (abs(eta) > 1.78 && abs(eta) <= 2.02) * (pt > 32.0 && pt <= 38.0) * (0.784336645237) +
   (abs(eta) > 1.78 && abs(eta) <= 2.02) * (pt > 38.0 && pt <= 44.0) * (0.628650954309) +
   (abs(eta) > 1.78 && abs(eta) <= 2.02) * (pt > 44.0 && pt <= 50.0) * (0.521505376344) +
   (abs(eta) > 1.78 && abs(eta) <= 2.02) * (pt > 50.0 && pt <= 60.0) * (0.726041666667) +
   (abs(eta) > 1.78 && abs(eta) <= 2.02) * (pt > 60.0 && pt <= 70.0) * (0.821764911057) +
   (abs(eta) > 1.78 && abs(eta) <= 2.02) * (pt > 70.0 && pt <= 80.0) * (0.687034793145) +
   (abs(eta) > 1.78 && abs(eta) <= 2.02) * (pt > 80.0 && pt <= 90.0) * (0.670454545455) +
   (abs(eta) > 1.78 && abs(eta) <= 2.02) * (pt > 90.0 && pt <= 100.0) * (0.636713735558) +
   (abs(eta) > 1.78 && abs(eta) <= 2.02) * (pt > 100.0 && pt <= 150.0) * (0.723881492143) +
   (abs(eta) > 1.78 && abs(eta) <= 2.02) * (pt > 150.0 && pt <= 200.0) * (0.775281276238) +
   (abs(eta) > 1.78 && abs(eta) <= 2.02) * (pt > 200.0 && pt <= 300.0) * (0.644915383443) +
   (abs(eta) > 1.78 && abs(eta) <= 2.02) * (pt > 300.0 && pt <= 400.0) * (0.608727528667) +
   (abs(eta) > 1.78 && abs(eta) <= 2.02) * (pt > 400.0) * (0.566928134082) +
   (abs(eta) > 2.02 && abs(eta) <= 2.26) * (pt > 20.0 && pt <= 26.0) * (0.366119149064) +
   (abs(eta) > 2.02 && abs(eta) <= 2.26) * (pt > 26.0 && pt <= 32.0) * (0.447434113151) +
   (abs(eta) > 2.02 && abs(eta) <= 2.26) * (pt > 32.0 && pt <= 38.0) * (0.566811091854) +
   (abs(eta) > 2.02 && abs(eta) <= 2.26) * (pt > 38.0 && pt <= 44.0) * (0.631193643425) +
   (abs(eta) > 2.02 && abs(eta) <= 2.26) * (pt > 44.0 && pt <= 50.0) * (0.653942428035) +
   (abs(eta) > 2.02 && abs(eta) <= 2.26) * (pt > 50.0 && pt <= 60.0) * (0.758013329102) +
   (abs(eta) > 2.02 && abs(eta) <= 2.26) * (pt > 60.0 && pt <= 70.0) * (0.710869206663) +
   (abs(eta) > 2.02 && abs(eta) <= 2.26) * (pt > 70.0 && pt <= 80.0) * (0.808382363623) +
   (abs(eta) > 2.02 && abs(eta) <= 2.26) * (pt > 80.0 && pt <= 90.0) * (0.829787234043) +
   (abs(eta) > 2.02 && abs(eta) <= 2.26) * (pt > 90.0 && pt <= 100.0) * (0.683744465528) +
   (abs(eta) > 2.02 && abs(eta) <= 2.26) * (pt > 100.0 && pt <= 150.0) * (0.680985221675) +
   (abs(eta) > 2.02 && abs(eta) <= 2.26) * (pt > 150.0 && pt <= 200.0) * (0.614848793748) +
   (abs(eta) > 2.02 && abs(eta) <= 2.26) * (pt > 200.0 && pt <= 300.0) * (0.637782214714) +
   (abs(eta) > 2.02 && abs(eta) <= 2.26) * (pt > 300.0 && pt <= 400.0) * (0.584644194757) +
   (abs(eta) > 2.02 && abs(eta) <= 2.26) * (pt > 400.0) * (0.544117647059) +
   (abs(eta) > 2.26 && abs(eta) <= 2.5) * (pt > 20.0 && pt <= 26.0) * (0.445076586433) +
   (abs(eta) > 2.26 && abs(eta) <= 2.5) * (pt > 26.0 && pt <= 32.0) * (0.478021978022) +
   (abs(eta) > 2.26 && abs(eta) <= 2.5) * (pt > 32.0 && pt <= 38.0) * (0.696459016393) +
   (abs(eta) > 2.26 && abs(eta) <= 2.5) * (pt > 38.0 && pt <= 44.0) * (0.746429381093) +
   (abs(eta) > 2.26 && abs(eta) <= 2.5) * (pt > 44.0 && pt <= 50.0) * (0.716761939886) +
   (abs(eta) > 2.26 && abs(eta) <= 2.5) * (pt > 50.0 && pt <= 60.0) * (0.743937563876) +
   (abs(eta) > 2.26 && abs(eta) <= 2.5) * (pt > 60.0 && pt <= 70.0) * (0.668055555556) +
   (abs(eta) > 2.26 && abs(eta) <= 2.5) * (pt > 70.0 && pt <= 80.0) * (0.735176619008) +
   (abs(eta) > 2.26 && abs(eta) <= 2.5) * (pt > 80.0 && pt <= 90.0) * (0.67365967366) +
   (abs(eta) > 2.26 && abs(eta) <= 2.5) * (pt > 90.0 && pt <= 100.0) * (0.760406424123) +
   (abs(eta) > 2.26 && abs(eta) <= 2.5) * (pt > 100.0 && pt <= 150.0) * (0.749638242894) +
   (abs(eta) > 2.26 && abs(eta) <= 2.5) * (pt > 150.0 && pt <= 200.0) * (0.622641509434) +
   (abs(eta) > 2.26 && abs(eta) <= 2.5) * (pt > 200.0 && pt <= 300.0) * (0.670033670034) +
   (abs(eta) > 2.26 && abs(eta) <= 2.5) * (pt > 300.0 && pt <= 400.0) * (0.686488169365) +
   (abs(eta) > 2.26 && abs(eta) <= 2.5) * (pt > 400.0) * (0.583333333333) +
   (abs(eta) > 2.5 && abs(eta) <= 2.6) * (pt > 20.0 && pt <= 26.0) * (0.883458646617) +
   (abs(eta) > 2.5 && abs(eta) <= 2.6) * (pt > 26.0 && pt <= 32.0) * (0.578534031414) +
   (abs(eta) > 2.5 && abs(eta) <= 2.6) * (pt > 32.0 && pt <= 38.0) * (0.60960960961) +
   (abs(eta) > 2.5 && abs(eta) <= 2.6) * (pt > 38.0 && pt <= 44.0) * (0.87910501624) +
   (abs(eta) > 2.5 && abs(eta) <= 2.6) * (pt > 44.0 && pt <= 50.0) * (0.797720797721) +
   (abs(eta) > 2.5 && abs(eta) <= 2.6) * (pt > 50.0 && pt <= 60.0) * (0.719476744186) +
   (abs(eta) > 2.5 && abs(eta) <= 2.6) * (pt > 60.0 && pt <= 70.0) * (0.799515445185) +
   (abs(eta) > 2.5 && abs(eta) <= 2.6) * (pt > 70.0 && pt <= 80.0) * (0.611650485437) +
   (abs(eta) > 2.5 && abs(eta) <= 2.6) * (pt > 80.0 && pt <= 90.0) * (0.823529411765) +
   (abs(eta) > 2.5 && abs(eta) <= 2.6) * (pt > 90.0 && pt <= 100.0) * (0.336879432624) +
   (abs(eta) > 2.5 && abs(eta) <= 2.6) * (pt > 100.0 && pt <= 150.0) * (0.766780025058) +
   (abs(eta) > 2.5 && abs(eta) <= 2.6) * (pt > 150.0 && pt <= 200.0) * (0.653846153846) +
   (abs(eta) > 2.5 && abs(eta) <= 2.6) * (pt > 200.0 && pt <= 300.0) * (0.673611111111) +
   (abs(eta) > 2.5 && abs(eta) <= 2.6) * (pt > 300.0 && pt <= 400.0) * (0.5) +
   (abs(eta) > 2.5 && abs(eta) <= 2.6) * (pt > 400.0) * (0.333333333333) +
   (abs(eta) > 2.6 && abs(eta) <= 2.7) * (pt > 20.0 && pt <= 26.0) * (0.586092715232) +
   (abs(eta) > 2.6 && abs(eta) <= 2.7) * (pt > 26.0 && pt <= 32.0) * (0.748412040873) +
   (abs(eta) > 2.6 && abs(eta) <= 2.7) * (pt > 32.0 && pt <= 38.0) * (0.812154696133) +
   (abs(eta) > 2.6 && abs(eta) <= 2.7) * (pt > 38.0 && pt <= 44.0) * (0.644016227181) +
   (abs(eta) > 2.6 && abs(eta) <= 2.7) * (pt > 44.0 && pt <= 50.0) * (0.971715328467) +
   (abs(eta) > 2.6 && abs(eta) <= 2.7) * (pt > 50.0 && pt <= 60.0) * (0.617045454545) +
   (abs(eta) > 2.6 && abs(eta) <= 2.7) * (pt > 60.0 && pt <= 70.0) * (0.821829163072) +
   (abs(eta) > 2.6 && abs(eta) <= 2.7) * (pt > 70.0 && pt <= 80.0) * (0.681818181818) +
   (abs(eta) > 2.6 && abs(eta) <= 2.7) * (pt > 80.0 && pt <= 90.0) * (0.815384615385) +
   (abs(eta) > 2.6 && abs(eta) <= 2.7) * (pt > 90.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 2.6 && abs(eta) <= 2.7) * (pt > 100.0 && pt <= 150.0) * (0.720190779014) +
   (abs(eta) > 2.6 && abs(eta) <= 2.7) * (pt > 150.0 && pt <= 200.0) * (0.592079969243) +
   (abs(eta) > 2.6 && abs(eta) <= 2.7) * (pt > 200.0 && pt <= 300.0) * (0.5) +
   (abs(eta) > 2.6 && abs(eta) <= 2.7) * (pt > 300.0 && pt <= 400.0) * (0.338235294118) +
   (abs(eta) > 2.6 && abs(eta) <= 2.7) * (pt > 400.0) * (1.0) +
   (abs(eta) > 2.7 && abs(eta) <= 2.8) * (pt > 20.0 && pt <= 26.0) * (0.605580843064) +
   (abs(eta) > 2.7 && abs(eta) <= 2.8) * (pt > 26.0 && pt <= 32.0) * (0.660406403941) +
   (abs(eta) > 2.7 && abs(eta) <= 2.8) * (pt > 32.0 && pt <= 38.0) * (0.862090290662) +
   (abs(eta) > 2.7 && abs(eta) <= 2.8) * (pt > 38.0 && pt <= 44.0) * (0.928034371643) +
   (abs(eta) > 2.7 && abs(eta) <= 2.8) * (pt > 44.0 && pt <= 50.0) * (0.825958702065) +
   (abs(eta) > 2.7 && abs(eta) <= 2.8) * (pt > 50.0 && pt <= 60.0) * (0.798214285714) +
   (abs(eta) > 2.7 && abs(eta) <= 2.8) * (pt > 60.0 && pt <= 70.0) * (0.801346801347) +
   (abs(eta) > 2.7 && abs(eta) <= 2.8) * (pt > 70.0 && pt <= 80.0) * (0.782726045884) +
   (abs(eta) > 2.7 && abs(eta) <= 2.8) * (pt > 80.0 && pt <= 90.0) * (0.72188449848) +
   (abs(eta) > 2.7 && abs(eta) <= 2.8) * (pt > 90.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 2.7 && abs(eta) <= 2.8) * (pt > 100.0 && pt <= 150.0) * (0.757411114288) +
   (abs(eta) > 2.7 && abs(eta) <= 2.8) * (pt > 150.0 && pt <= 200.0) * (0.555555555556) +
   (abs(eta) > 2.7 && abs(eta) <= 2.8) * (pt > 200.0 && pt <= 300.0) * (0.696985446985) +
   (abs(eta) > 2.7 && abs(eta) <= 2.8) * (pt > 300.0 && pt <= 400.0) * (0.5) +
   (abs(eta) > 2.7 && abs(eta) <= 2.8) * (pt > 400.0) * (0.0) +
   (abs(eta) > 2.8 && abs(eta) <= 2.9) * (pt > 20.0 && pt <= 26.0) * (0.198859166011) +
   (abs(eta) > 2.8 && abs(eta) <= 2.9) * (pt > 26.0 && pt <= 32.0) * (0.647429519071) +
   (abs(eta) > 2.8 && abs(eta) <= 2.9) * (pt > 32.0 && pt <= 38.0) * (0.847368421053) +
   (abs(eta) > 2.8 && abs(eta) <= 2.9) * (pt > 38.0 && pt <= 44.0) * (0.740236686391) +
   (abs(eta) > 2.8 && abs(eta) <= 2.9) * (pt > 44.0 && pt <= 50.0) * (0.738362760835) +
   (abs(eta) > 2.8 && abs(eta) <= 2.9) * (pt > 50.0 && pt <= 60.0) * (0.857843137255) +
   (abs(eta) > 2.8 && abs(eta) <= 2.9) * (pt > 60.0 && pt <= 70.0) * (0.727131782946) +
   (abs(eta) > 2.8 && abs(eta) <= 2.9) * (pt > 70.0 && pt <= 80.0) * (0.633802816901) +
   (abs(eta) > 2.8 && abs(eta) <= 2.9) * (pt > 80.0 && pt <= 90.0) * (0.714285714286) +
   (abs(eta) > 2.8 && abs(eta) <= 2.9) * (pt > 90.0 && pt <= 100.0) * (0.810256410256) +
   (abs(eta) > 2.8 && abs(eta) <= 2.9) * (pt > 100.0 && pt <= 150.0) * (0.540887040887) +
   (abs(eta) > 2.8 && abs(eta) <= 2.9) * (pt > 150.0 && pt <= 200.0) * (0.666666666667) +
   (abs(eta) > 2.8 && abs(eta) <= 2.9) * (pt > 200.0 && pt <= 300.0) * (0.428571428571) +
   (abs(eta) > 2.8 && abs(eta) <= 2.9) * (pt > 300.0 && pt <= 400.0) * (0.5) +
   (abs(eta) > 2.8 && abs(eta) <= 2.9) * (pt > 400.0) * (0.0) +
   (abs(eta) > 2.9 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 26.0) * (0.398034398034) +
   (abs(eta) > 2.9 && abs(eta) <= 3.0) * (pt > 26.0 && pt <= 32.0) * (0.515294117647) +
   (abs(eta) > 2.9 && abs(eta) <= 3.0) * (pt > 32.0 && pt <= 38.0) * (0.493665540541) +
   (abs(eta) > 2.9 && abs(eta) <= 3.0) * (pt > 38.0 && pt <= 44.0) * (0.643570170288) +
   (abs(eta) > 2.9 && abs(eta) <= 3.0) * (pt > 44.0 && pt <= 50.0) * (0.727272727273) +
   (abs(eta) > 2.9 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 60.0) * (0.514492753623) +
   (abs(eta) > 2.9 && abs(eta) <= 3.0) * (pt > 60.0 && pt <= 70.0) * (0.433613445378) +
   (abs(eta) > 2.9 && abs(eta) <= 3.0) * (pt > 70.0 && pt <= 80.0) * (0.477453580902) +
   (abs(eta) > 2.9 && abs(eta) <= 3.0) * (pt > 80.0 && pt <= 90.0) * (0.926470588235) +
   (abs(eta) > 2.9 && abs(eta) <= 3.0) * (pt > 90.0 && pt <= 100.0) * (0.445054945055) +
   (abs(eta) > 2.9 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 150.0) * (0.737662337662) +
   (abs(eta) > 2.9 && abs(eta) <= 3.0) * (pt > 150.0 && pt <= 200.0) * (0.691588785047) +
   (abs(eta) > 2.9 && abs(eta) <= 3.0) * (pt > 200.0 && pt <= 300.0) * (0.333333333333) +
   (abs(eta) > 2.9 && abs(eta) <= 3.0) * (pt > 300.0 && pt <= 400.0) * (0.666666666667) +
   (abs(eta) > 2.9 && abs(eta) <= 3.0) * (pt > 400.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 20.0 && pt <= 26.0) * (0.0863067802201) +
   (abs(eta) > 3.0) * (pt > 26.0 && pt <= 32.0) * (0.142522011121) +
   (abs(eta) > 3.0) * (pt > 32.0 && pt <= 38.0) * (0.176817716334) +
   (abs(eta) > 3.0) * (pt > 38.0 && pt <= 44.0) * (0.345163102949) +
   (abs(eta) > 3.0) * (pt > 44.0 && pt <= 50.0) * (0.266826611334) +
   (abs(eta) > 3.0) * (pt > 50.0 && pt <= 60.0) * (0.351033560728) +
   (abs(eta) > 3.0) * (pt > 60.0 && pt <= 70.0) * (0.378212805159) +
   (abs(eta) > 3.0) * (pt > 70.0 && pt <= 80.0) * (0.355743603412) +
   (abs(eta) > 3.0) * (pt > 80.0 && pt <= 90.0) * (0.318655216517) +
   (abs(eta) > 3.0) * (pt > 90.0 && pt <= 100.0) * (0.286337843736) +
   (abs(eta) > 3.0) * (pt > 100.0 && pt <= 150.0) * (0.253841890462) +
   (abs(eta) > 3.0) * (pt > 150.0 && pt <= 200.0) * (0.277247956403) +
   (abs(eta) > 3.0) * (pt > 200.0 && pt <= 300.0) * (0.347826086957) +
   (abs(eta) > 3.0) * (pt > 300.0 && pt <= 400.0) * (0.5) +
   (abs(eta) > 3.0) * (pt > 400.0) * (1.0)  
  }

}



#################
# Photon filter
#################

module PdgCodeFilter PhotonFilter {
  set InputArray PhotonEnergySmearing/eflowPhotons
  set OutputArray photons
  set Invert true
  set PTMin 5.0
  add PdgCode {22}
}



####################
# Photon isolation #
####################

module Isolation PhotonIsolation {

  # particle for which calculate the isolation
  set CandidateInputArray PhotonFilter/photons

  # isolation collection
  set IsolationInputArray EFlowFilterPuppi/eflow

  # output array
  set OutputArray photons

  # veto isolation cand. based on proximity to input cand.
  set DeltaRMin 0.01
  set UseMiniCone true

  # isolation cone
  set DeltaRMax 0.3

  # minimum pT
  set PTMin     0.0

  # iso ratio to cut
  set PTRatioMax 9999.

}

##################
# Photon scale #
##################


module EnergyScale PhotonScale {
  set InputArray PhotonIsolation/photons
  set OutputArray photons

    ### photon looseIDISO momentum scale
  set ScaleFormula {
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 20.0) * (1.004) +
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 50.0) * (0.995) +
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 100.0) * (0.997) +
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 150.0) * (0.998) +
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 150.0 && pt <= 14000.0) * (0.997) +
   (abs(eta) > 1.5 && abs(eta) <= 3.0) * (pt > 10.0 && pt <= 20.0) * (1.229) +
   (abs(eta) > 1.5 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 50.0) * (1.034) +
   (abs(eta) > 1.5 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 100.0) * (0.987) +
   (abs(eta) > 1.5 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 150.0) * (0.991) +
   (abs(eta) > 1.5 && abs(eta) <= 3.0) * (pt > 150.0 && pt <= 14000.0) * (0.921)  
  }

}


##################
# Photon smear #
##################

module MomentumSmearing PhotonSmear {

  set InputArray PhotonScale/photons
  set OutputArray photons

    ### photon looseIDISO momentum resolution
  set ResolutionFormula {
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 20.0) * (0.029120) +
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 50.0) * (0.000001) +
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 100.0) * (0.000001) +
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 150.0) * (0.000001) +
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 150.0 && pt <= 14000.0) * (0.000001) +
   (abs(eta) > 1.5 && abs(eta) <= 3.0) * (pt > 10.0 && pt <= 20.0) * (0.076837) +
   (abs(eta) > 1.5 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 50.0) * (0.058847) +
   (abs(eta) > 1.5 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 100.0) * (0.041388) +
   (abs(eta) > 1.5 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 150.0) * (0.010490) +
   (abs(eta) > 1.5 && abs(eta) <= 3.0) * (pt > 150.0 && pt <= 14000.0) * (0.106903)  
  }

}




#####################
# Photon Id Loose   #
#####################

module Efficiency PhotonLooseID {

  ## input particles
  set InputArray PhotonSmear/photons
  ## output particles
  set OutputArray photons
  # set EfficiencyFormula {efficiency formula as a function of eta and pt}
  # efficiency formula for photons


    ### photon loose ID 

  set EfficiencyFormula {

   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 10.0 && pt <= 12.0) * (0.614285714286) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 12.0 && pt <= 14.0) * (0.794117647059) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 14.0 && pt <= 16.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 16.0 && pt <= 18.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 18.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 20.0 && pt <= 26.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 26.0 && pt <= 32.0) * (0.993372660403) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 32.0 && pt <= 38.0) * (0.961348363321) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 38.0 && pt <= 44.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 44.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 50.0 && pt <= 60.0) * (0.995573149449) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 60.0 && pt <= 70.0) * (0.985067081106) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 70.0 && pt <= 80.0) * (0.997375328084) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 90.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 100.0 && pt <= 125.0) * (0.955219080428) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 125.0 && pt <= 150.0) * (0.975811764706) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 150.0) * (0.976725333137) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 10.0 && pt <= 12.0) * (0.907407407407) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 12.0 && pt <= 14.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 14.0 && pt <= 16.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 16.0 && pt <= 18.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 18.0 && pt <= 20.0) * (0.915254237288) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 20.0 && pt <= 26.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 26.0 && pt <= 32.0) * (0.977951635846) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 32.0 && pt <= 38.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 38.0 && pt <= 44.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 44.0 && pt <= 50.0) * (0.946387238036) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 50.0 && pt <= 60.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 60.0 && pt <= 70.0) * (0.967253376995) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 80.0 && pt <= 90.0) * (0.999750918488) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 90.0 && pt <= 100.0) * (0.990056569362) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 100.0 && pt <= 125.0) * (0.971517703266) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 125.0 && pt <= 150.0) * (0.9598766122) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 150.0) * (0.978509445901) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 10.0 && pt <= 12.0) * (0.658536585366) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 12.0 && pt <= 14.0) * (0.877611940299) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 14.0 && pt <= 16.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 16.0 && pt <= 18.0) * (0.963662790698) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 18.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 20.0 && pt <= 26.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 26.0 && pt <= 32.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 32.0 && pt <= 38.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 38.0 && pt <= 44.0) * (0.973509933775) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 44.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 50.0 && pt <= 60.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 60.0 && pt <= 70.0) * (0.940157700541) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 70.0 && pt <= 80.0) * (0.984961743683) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 80.0 && pt <= 90.0) * (0.925310692753) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 90.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 100.0 && pt <= 125.0) * (0.982647096828) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 125.0 && pt <= 150.0) * (0.980462936224) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 150.0) * (0.959225579349) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 10.0 && pt <= 12.0) * (0.7) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 12.0 && pt <= 14.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 14.0 && pt <= 16.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 16.0 && pt <= 18.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 18.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 20.0 && pt <= 26.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 26.0 && pt <= 32.0) * (0.986008610086) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 32.0 && pt <= 38.0) * (0.955434782609) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 38.0 && pt <= 44.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 44.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 50.0 && pt <= 60.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 60.0 && pt <= 70.0) * (0.929412314304) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 70.0 && pt <= 80.0) * (0.983133833773) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 80.0 && pt <= 90.0) * (0.99203187251) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 90.0 && pt <= 100.0) * (0.986634600221) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 100.0 && pt <= 125.0) * (0.94835832317) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 125.0 && pt <= 150.0) * (0.965795508898) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 150.0) * (0.98146880442) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 12.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 12.0 && pt <= 14.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 14.0 && pt <= 16.0) * (0.659932659933) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 16.0 && pt <= 18.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 18.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 26.0) * (0.916) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 26.0 && pt <= 32.0) * (0.968850698174) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 32.0 && pt <= 38.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 38.0 && pt <= 44.0) * (0.999689296256) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 44.0 && pt <= 50.0) * (0.970893643801) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 60.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 60.0 && pt <= 70.0) * (0.981879194631) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 70.0 && pt <= 80.0) * (0.961545738341) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 80.0 && pt <= 90.0) * (0.982362474974) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 90.0 && pt <= 100.0) * (0.987160955701) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 125.0) * (0.928305871957) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 125.0 && pt <= 150.0) * (0.973671497585) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 150.0) * (0.942050998803) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 10.0 && pt <= 12.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 12.0 && pt <= 14.0) * (0.141949152542) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 14.0 && pt <= 16.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 16.0 && pt <= 18.0) * (0.285294117647) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 18.0 && pt <= 20.0) * (0.226373626374) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 20.0 && pt <= 26.0) * (0.324552341598) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 26.0 && pt <= 32.0) * (0.381155555556) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 32.0 && pt <= 38.0) * (0.383098591549) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 38.0 && pt <= 44.0) * (0.440639269406) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 44.0 && pt <= 50.0) * (0.341989047871) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60.0) * (0.22720917226) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70.0) * (0.156269691241) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80.0) * (0.21052975374) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90.0) * (0.10952420841) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100.0) * (0.0634457864779) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 100.0 && pt <= 125.0) * (0.0968098111124) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 125.0 && pt <= 150.0) * (0.124135741499) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 150.0) * (0.081693152298) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 10.0 && pt <= 12.0) * (0.0) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 12.0 && pt <= 14.0) * (0.0) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 14.0 && pt <= 16.0) * (0.449315068493) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 16.0 && pt <= 18.0) * (0.572368421053) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 18.0 && pt <= 20.0) * (0.0) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 20.0 && pt <= 26.0) * (0.144885679904) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 26.0 && pt <= 32.0) * (0.030784030784) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 32.0 && pt <= 38.0) * (0.100647993311) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 38.0 && pt <= 44.0) * (0.0) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 44.0 && pt <= 50.0) * (0.0) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 50.0 && pt <= 60.0) * (0.0198549789302) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 60.0 && pt <= 70.0) * (0.0386243386243) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 70.0 && pt <= 80.0) * (0.025843373494) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 80.0 && pt <= 90.0) * (0.0) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 100.0 && pt <= 125.0) * (0.0143965714523) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 125.0 && pt <= 150.0) * (0.0175768147055) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 150.0) * (0.0) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 10.0 && pt <= 12.0) * (0.0) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 12.0 && pt <= 14.0) * (0.0) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 14.0 && pt <= 16.0) * (0.0) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 16.0 && pt <= 18.0) * (0.0) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 18.0 && pt <= 20.0) * (0.157384987893) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 26.0) * (0.0) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 26.0 && pt <= 32.0) * (0.124803767661) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 32.0 && pt <= 38.0) * (0.0) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 38.0 && pt <= 44.0) * (0.0) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 44.0 && pt <= 50.0) * (0.0435267857143) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60.0) * (0.0276087617889) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70.0) * (0.0385078219013) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80.0) * (0.0) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90.0) * (0.0523076923077) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100.0) * (0.0408888888889) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 125.0) * (0.0455657492355) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 150.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 10.0 && pt <= 12.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 12.0 && pt <= 14.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 14.0 && pt <= 16.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 16.0 && pt <= 18.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 18.0 && pt <= 20.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 20.0 && pt <= 26.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 26.0 && pt <= 32.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 32.0 && pt <= 38.0) * (0.11377245509) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 38.0 && pt <= 44.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 44.0 && pt <= 50.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 50.0 && pt <= 60.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 60.0 && pt <= 70.0) * (0.0645380434783) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 70.0 && pt <= 80.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 80.0 && pt <= 90.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 100.0 && pt <= 125.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 150.0) * (0.0467607105538) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 10.0 && pt <= 12.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 12.0 && pt <= 14.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 14.0 && pt <= 16.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 16.0 && pt <= 18.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 18.0 && pt <= 20.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 26.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 26.0 && pt <= 32.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 32.0 && pt <= 38.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 38.0 && pt <= 44.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 44.0 && pt <= 50.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 60.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 60.0 && pt <= 70.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 70.0 && pt <= 80.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 80.0 && pt <= 90.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 125.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 150.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 10.0 && pt <= 12.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 12.0 && pt <= 14.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 14.0 && pt <= 16.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 16.0 && pt <= 18.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 18.0 && pt <= 20.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 20.0 && pt <= 26.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 26.0 && pt <= 32.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 32.0 && pt <= 38.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 38.0 && pt <= 44.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 44.0 && pt <= 50.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 50.0 && pt <= 60.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 60.0 && pt <= 70.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 70.0 && pt <= 80.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 80.0 && pt <= 90.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 100.0 && pt <= 125.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 150.0) * (0.0)  
  }

}


#####################
# Photon Id Medium   #
#####################

module Efficiency PhotonMediumID {

  ## input particles
  set InputArray PhotonSmear/photons
  ## output particles
  set OutputArray photons
  # set EfficiencyFormula {efficiency formula as a function of eta and pt}
  # efficiency formula for photons

    ### photon medium ID 

  set EfficiencyFormula {

   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 10.0 && pt <= 12.0) * (0.651515151515) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 12.0 && pt <= 14.0) * (0.692307692308) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 14.0 && pt <= 16.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 16.0 && pt <= 18.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 18.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 20.0 && pt <= 26.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 26.0 && pt <= 32.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 32.0 && pt <= 38.0) * (0.973660714286) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 38.0 && pt <= 44.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 44.0 && pt <= 50.0) * (0.983055975794) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 50.0 && pt <= 60.0) * (0.95568350909) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 60.0 && pt <= 70.0) * (0.968289856174) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 70.0 && pt <= 80.0) * (0.99432012698) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 80.0 && pt <= 90.0) * (0.998484376622) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 90.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 100.0 && pt <= 125.0) * (0.945705841965) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 125.0 && pt <= 150.0) * (0.934771784232) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 150.0) * (0.967616071532) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 10.0 && pt <= 12.0) * (0.494949494949) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 12.0 && pt <= 14.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 14.0 && pt <= 16.0) * (0.949119373777) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 16.0 && pt <= 18.0) * (0.839215686275) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 18.0 && pt <= 20.0) * (0.947368421053) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 20.0 && pt <= 26.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 26.0 && pt <= 32.0) * (0.999563826694) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 32.0 && pt <= 38.0) * (0.997037037037) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 38.0 && pt <= 44.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 44.0 && pt <= 50.0) * (0.96726342711) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 50.0 && pt <= 60.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 60.0 && pt <= 70.0) * (0.97169158361) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 70.0 && pt <= 80.0) * (0.996466140697) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 80.0 && pt <= 90.0) * (0.967800805504) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 90.0 && pt <= 100.0) * (0.987975077882) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 100.0 && pt <= 125.0) * (0.966803997195) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 125.0 && pt <= 150.0) * (0.925094920625) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 150.0) * (0.967476578134) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 10.0 && pt <= 12.0) * (0.3375) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 12.0 && pt <= 14.0) * (0.904615384615) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 14.0 && pt <= 16.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 16.0 && pt <= 18.0) * (0.963662790698) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 18.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 20.0 && pt <= 26.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 26.0 && pt <= 32.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 32.0 && pt <= 38.0) * (0.933356377463) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 38.0 && pt <= 44.0) * (0.953468697124) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 44.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 50.0 && pt <= 60.0) * (0.969512195122) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 60.0 && pt <= 70.0) * (0.942164440564) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 70.0 && pt <= 80.0) * (0.953584791182) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 80.0 && pt <= 90.0) * (0.906192065258) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 90.0 && pt <= 100.0) * (0.969960963534) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 100.0 && pt <= 125.0) * (0.971049364741) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 125.0 && pt <= 150.0) * (0.960470561998) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 150.0) * (0.95) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 10.0 && pt <= 12.0) * (0.742424242424) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 12.0 && pt <= 14.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 14.0 && pt <= 16.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 16.0 && pt <= 18.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 18.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 20.0 && pt <= 26.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 26.0 && pt <= 32.0) * (0.98708994709) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 32.0 && pt <= 38.0) * (0.944290976059) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 38.0 && pt <= 44.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 44.0 && pt <= 50.0) * (0.990315006848) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 50.0 && pt <= 60.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 60.0 && pt <= 70.0) * (0.94183459311) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 70.0 && pt <= 80.0) * (0.970556489186) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 90.0 && pt <= 100.0) * (0.95756020529) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 100.0 && pt <= 125.0) * (0.935853875997) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 125.0 && pt <= 150.0) * (0.935929803718) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 150.0) * (0.964528326199) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 12.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 12.0 && pt <= 14.0) * (0.948979591837) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 14.0 && pt <= 16.0) * (0.702508960573) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 16.0 && pt <= 18.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 18.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 26.0) * (0.889896373057) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 26.0 && pt <= 32.0) * (0.999538276849) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 32.0 && pt <= 38.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 38.0 && pt <= 44.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 44.0 && pt <= 50.0) * (0.986097794823) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 60.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 60.0 && pt <= 70.0) * (0.967828418231) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 70.0 && pt <= 80.0) * (0.935245628246) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 80.0 && pt <= 90.0) * (0.976509097911) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 90.0 && pt <= 100.0) * (0.993581514763) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 125.0) * (0.914949154386) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 125.0 && pt <= 150.0) * (0.962491853645) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 150.0) * (0.942275965634) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 10.0 && pt <= 12.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 12.0 && pt <= 14.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 14.0 && pt <= 16.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 16.0 && pt <= 18.0) * (0.146084337349) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 18.0 && pt <= 20.0) * (0.236781609195) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 20.0 && pt <= 26.0) * (0.264909969258) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 26.0 && pt <= 32.0) * (0.340781108084) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 32.0 && pt <= 38.0) * (0.250712250712) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 38.0 && pt <= 44.0) * (0.269553072626) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 44.0 && pt <= 50.0) * (0.281499911143) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60.0) * (0.175957207207) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70.0) * (0.110155080078) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80.0) * (0.154474759457) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90.0) * (0.0881398778176) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 100.0 && pt <= 125.0) * (0.0781069795188) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 125.0 && pt <= 150.0) * (0.0833155763904) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 150.0) * (0.0445279357576) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 10.0 && pt <= 12.0) * (0.0) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 12.0 && pt <= 14.0) * (0.0) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 14.0 && pt <= 16.0) * (0.230985915493) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 16.0 && pt <= 18.0) * (0.453125) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 18.0 && pt <= 20.0) * (0.0) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 20.0 && pt <= 26.0) * (0.0372986369269) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 26.0 && pt <= 32.0) * (0.03120429059) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 32.0 && pt <= 38.0) * (0.0682397959184) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 38.0 && pt <= 44.0) * (0.0) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 44.0 && pt <= 50.0) * (0.0) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 50.0 && pt <= 60.0) * (0.0) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 60.0 && pt <= 70.0) * (0.0) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 70.0 && pt <= 80.0) * (0.0260948905109) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 80.0 && pt <= 90.0) * (0.0) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 100.0 && pt <= 125.0) * (0.0144821379989) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 150.0) * (0.0) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 10.0 && pt <= 12.0) * (0.0) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 12.0 && pt <= 14.0) * (0.0) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 14.0 && pt <= 16.0) * (0.0) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 16.0 && pt <= 18.0) * (0.0) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 18.0 && pt <= 20.0) * (0.160098522167) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 26.0) * (0.0) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 26.0 && pt <= 32.0) * (0.0) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 32.0 && pt <= 38.0) * (0.0) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 38.0 && pt <= 44.0) * (0.0) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 44.0 && pt <= 50.0) * (0.0441176470588) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60.0) * (0.0276087617889) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70.0) * (0.0389294403893) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80.0) * (0.0) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90.0) * (0.0) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100.0) * (0.0412556053812) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 125.0) * (0.045987654321) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 150.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 10.0 && pt <= 12.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 12.0 && pt <= 14.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 14.0 && pt <= 16.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 16.0 && pt <= 18.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 18.0 && pt <= 20.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 20.0 && pt <= 26.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 26.0 && pt <= 32.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 32.0 && pt <= 38.0) * (0.115151515152) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 38.0 && pt <= 44.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 44.0 && pt <= 50.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 50.0 && pt <= 60.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 60.0 && pt <= 70.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 70.0 && pt <= 80.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 80.0 && pt <= 90.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 100.0 && pt <= 125.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 150.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 10.0 && pt <= 12.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 12.0 && pt <= 14.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 14.0 && pt <= 16.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 16.0 && pt <= 18.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 18.0 && pt <= 20.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 26.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 26.0 && pt <= 32.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 32.0 && pt <= 38.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 38.0 && pt <= 44.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 44.0 && pt <= 50.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 60.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 60.0 && pt <= 70.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 70.0 && pt <= 80.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 80.0 && pt <= 90.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 125.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 150.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 10.0 && pt <= 12.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 12.0 && pt <= 14.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 14.0 && pt <= 16.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 16.0 && pt <= 18.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 18.0 && pt <= 20.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 20.0 && pt <= 26.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 26.0 && pt <= 32.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 32.0 && pt <= 38.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 38.0 && pt <= 44.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 44.0 && pt <= 50.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 50.0 && pt <= 60.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 60.0 && pt <= 70.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 70.0 && pt <= 80.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 80.0 && pt <= 90.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 100.0 && pt <= 125.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 150.0) * (0.0)  
  }

}


#####################
# Photon Id Tight   #
#####################

module Efficiency PhotonTightID {

  ## input particles
  set InputArray PhotonSmear/photons
  ## output particles
  set OutputArray photons
  # set EfficiencyFormula {efficiency formula as a function of eta and pt}
  # efficiency formula for photons

    ### photon tight ID 

  set EfficiencyFormula {

   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 10.0 && pt <= 12.0) * (0.530864197531) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 12.0 && pt <= 14.0) * (0.44262295082) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 14.0 && pt <= 16.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 16.0 && pt <= 18.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 18.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 20.0 && pt <= 26.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 26.0 && pt <= 32.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 32.0 && pt <= 38.0) * (0.968880056114) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 38.0 && pt <= 44.0) * (0.988898500577) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 44.0 && pt <= 50.0) * (0.911041009464) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 50.0 && pt <= 60.0) * (0.925660382977) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 60.0 && pt <= 70.0) * (0.903766676156) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 70.0 && pt <= 80.0) * (0.915827010985) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 80.0 && pt <= 90.0) * (0.961668521758) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 90.0 && pt <= 100.0) * (0.92035749337) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 100.0 && pt <= 125.0) * (0.840747374222) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 125.0 && pt <= 150.0) * (0.910769230769) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 150.0) * (0.909116751364) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 10.0 && pt <= 12.0) * (0.544444444444) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 12.0 && pt <= 14.0) * (0.959064327485) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 14.0 && pt <= 16.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 16.0 && pt <= 18.0) * (0.571581196581) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 18.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 20.0 && pt <= 26.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 26.0 && pt <= 32.0) * (0.980516194332) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 32.0 && pt <= 38.0) * (0.929638262972) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 38.0 && pt <= 44.0) * (0.976511056511) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 44.0 && pt <= 50.0) * (0.946129053476) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 50.0 && pt <= 60.0) * (0.943818838664) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 60.0 && pt <= 70.0) * (0.913428571429) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 70.0 && pt <= 80.0) * (0.915750915751) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 80.0 && pt <= 90.0) * (0.903782135263) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 90.0 && pt <= 100.0) * (0.897323271346) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 100.0 && pt <= 125.0) * (0.905244468118) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 125.0 && pt <= 150.0) * (0.90610537089) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 150.0) * (0.907336956522) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 10.0 && pt <= 12.0) * (0.364864864865) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 12.0 && pt <= 14.0) * (0.736842105263) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 14.0 && pt <= 16.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 16.0 && pt <= 18.0) * (0.865740740741) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 18.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 20.0 && pt <= 26.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 26.0 && pt <= 32.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 32.0 && pt <= 38.0) * (0.884003673095) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 38.0 && pt <= 44.0) * (0.995036101083) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 44.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 50.0 && pt <= 60.0) * (0.938883034773) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 60.0 && pt <= 70.0) * (0.912609114143) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 70.0 && pt <= 80.0) * (0.916411837877) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 80.0 && pt <= 90.0) * (0.893196503231) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 90.0 && pt <= 100.0) * (0.938399539436) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 100.0 && pt <= 125.0) * (0.951210022566) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 125.0 && pt <= 150.0) * (0.927251446193) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 150.0) * (0.916888494528) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 10.0 && pt <= 12.0) * (0.790322580645) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 12.0 && pt <= 14.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 14.0 && pt <= 16.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 16.0 && pt <= 18.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 18.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 20.0 && pt <= 26.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 26.0 && pt <= 32.0) * (0.945782103825) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 32.0 && pt <= 38.0) * (0.917416829746) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 38.0 && pt <= 44.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 44.0 && pt <= 50.0) * (0.985931174089) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 50.0 && pt <= 60.0) * (0.986556603774) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 60.0 && pt <= 70.0) * (0.92332325209) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 70.0 && pt <= 80.0) * (0.922776148583) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 80.0 && pt <= 90.0) * (0.931880108992) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 90.0 && pt <= 100.0) * (0.897219543555) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 100.0 && pt <= 125.0) * (0.875093118398) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 125.0 && pt <= 150.0) * (0.8924801947) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 150.0) * (0.912808635727) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 12.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 12.0 && pt <= 14.0) * (0.989361702128) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 14.0 && pt <= 16.0) * (0.734082397004) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 16.0 && pt <= 18.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 18.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 26.0) * (0.861142061281) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 26.0 && pt <= 32.0) * (0.989939637827) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 32.0 && pt <= 38.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 38.0 && pt <= 44.0) * (0.976705946835) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 44.0 && pt <= 50.0) * (0.928905206943) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 60.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 60.0 && pt <= 70.0) * (0.901650618982) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 70.0 && pt <= 80.0) * (0.874903993856) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 80.0 && pt <= 90.0) * (0.909427347098) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 90.0 && pt <= 100.0) * (0.894823549814) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 125.0) * (0.871111403624) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 125.0 && pt <= 150.0) * (0.937555816686) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 150.0) * (0.864389776105) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 10.0 && pt <= 12.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 12.0 && pt <= 14.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 14.0 && pt <= 16.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 16.0 && pt <= 18.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 18.0 && pt <= 20.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 20.0 && pt <= 26.0) * (0.173094582185) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 26.0 && pt <= 32.0) * (0.127619047619) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 32.0 && pt <= 38.0) * (0.0718562874251) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 38.0 && pt <= 44.0) * (0.125406107862) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 44.0 && pt <= 50.0) * (0.194726166329) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60.0) * (0.108318890815) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70.0) * (0.0481366459627) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80.0) * (0.039176360476) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90.0) * (0.0675768800069) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 100.0 && pt <= 125.0) * (0.0691528424695) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 125.0 && pt <= 150.0) * (0.0556187766714) +
   (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 150.0) * (0.0149498897991) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 10.0 && pt <= 12.0) * (0.0) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 12.0 && pt <= 14.0) * (0.0) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 14.0 && pt <= 16.0) * (0.0) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 16.0 && pt <= 18.0) * (0.162313432836) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 18.0 && pt <= 20.0) * (0.0) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 20.0 && pt <= 26.0) * (0.038441890166) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 26.0 && pt <= 32.0) * (0.0) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 32.0 && pt <= 38.0) * (0.0) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 38.0 && pt <= 44.0) * (0.0) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 44.0 && pt <= 50.0) * (0.0) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 50.0 && pt <= 60.0) * (0.0) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 60.0 && pt <= 70.0) * (0.0) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 70.0 && pt <= 80.0) * (0.0) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 80.0 && pt <= 90.0) * (0.0) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 100.0 && pt <= 125.0) * (0.0) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 150.0) * (0.0) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 10.0 && pt <= 12.0) * (0.0) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 12.0 && pt <= 14.0) * (0.0) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 14.0 && pt <= 16.0) * (0.0) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 16.0 && pt <= 18.0) * (0.0) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 18.0 && pt <= 20.0) * (0.0) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 26.0) * (0.0) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 26.0 && pt <= 32.0) * (0.0) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 32.0 && pt <= 38.0) * (0.0) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 38.0 && pt <= 44.0) * (0.0) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 44.0 && pt <= 50.0) * (0.0) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60.0) * (0.0) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70.0) * (0.0) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80.0) * (0.0) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90.0) * (0.0) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 125.0) * (0.0) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 150.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 10.0 && pt <= 12.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 12.0 && pt <= 14.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 14.0 && pt <= 16.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 16.0 && pt <= 18.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 18.0 && pt <= 20.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 20.0 && pt <= 26.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 26.0 && pt <= 32.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 32.0 && pt <= 38.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 38.0 && pt <= 44.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 44.0 && pt <= 50.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 50.0 && pt <= 60.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 60.0 && pt <= 70.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 70.0 && pt <= 80.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 80.0 && pt <= 90.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 100.0 && pt <= 125.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 150.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 10.0 && pt <= 12.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 12.0 && pt <= 14.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 14.0 && pt <= 16.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 16.0 && pt <= 18.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 18.0 && pt <= 20.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 26.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 26.0 && pt <= 32.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 32.0 && pt <= 38.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 38.0 && pt <= 44.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 44.0 && pt <= 50.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 60.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 60.0 && pt <= 70.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 70.0 && pt <= 80.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 80.0 && pt <= 90.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 125.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 150.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 10.0 && pt <= 12.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 12.0 && pt <= 14.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 14.0 && pt <= 16.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 16.0 && pt <= 18.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 18.0 && pt <= 20.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 20.0 && pt <= 26.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 26.0 && pt <= 32.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 32.0 && pt <= 38.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 38.0 && pt <= 44.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 44.0 && pt <= 50.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 50.0 && pt <= 60.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 60.0 && pt <= 70.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 70.0 && pt <= 80.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 80.0 && pt <= 90.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 100.0 && pt <= 125.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 150.0) * (0.0)  
  }

}


######################
# Electron isolation #
######################

module Isolation ElectronIsolation {

  set CandidateInputArray ElectronFilter/electrons

  # isolation collection
  set IsolationInputArray EFlowFilterPuppi/eflow

  set OutputArray electrons

  set DeltaRMax 0.3
  set PTMin 0.0
  set PTRatioMax 9999.

}


##################
# Electron scale #
##################


module EnergyScale ElectronScale {
  set InputArray ElectronIsolation/electrons
  set OutputArray electrons

    ### electron looseIDISO momentum scale
  set ScaleFormula {
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 20.0) * (1.002) +
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 50.0) * (0.999) +
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 100.0) * (0.997) +
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 150.0) * (0.774) +
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 150.0 && pt <= 14000.0) * (0.985) +
   (abs(eta) > 1.5 && abs(eta) <= 3.0) * (pt > 10.0 && pt <= 20.0) * (1.298) +
   (abs(eta) > 1.5 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 50.0) * (1.074) +
   (abs(eta) > 1.5 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 100.0) * (1.034) +
   (abs(eta) > 1.5 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 150.0) * (0.995) +
   (abs(eta) > 1.5 && abs(eta) <= 3.0) * (pt > 150.0 && pt <= 14000.0) * (1.027)  
  }

}

##################
# Electron smear #
##################

module MomentumSmearing ElectronSmear {

  set InputArray ElectronScale/electrons
  set OutputArray electrons

    ### electron looseIDISO momentum resolution
  set ResolutionFormula {
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 20.0) * (0.019891) +
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 50.0) * (0.022481) +
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 100.0) * (0.015265) +
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 150.0) * (0.013492) +
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 150.0 && pt <= 14000.0) * (0.050900) +
   (abs(eta) > 1.5 && abs(eta) <= 3.0) * (pt > 10.0 && pt <= 20.0) * (0.161362) +
   (abs(eta) > 1.5 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 50.0) * (0.068773) +
   (abs(eta) > 1.5 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 100.0) * (0.045916) +
   (abs(eta) > 1.5 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 150.0) * (0.005769) +
   (abs(eta) > 1.5 && abs(eta) <= 3.0) * (pt > 150.0 && pt <= 14000.0) * (0.246918)  
  }

}



#######################
# Electron loose ID efficiency #
#######################

module Efficiency ElectronLooseEfficiency {

  set InputArray ElectronSmear/electrons
  set OutputArray electrons

    ### electron loose ID 

  set EfficiencyFormula {

   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 10.0 && pt <= 12.0) * (0.931286549708) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 12.0 && pt <= 14.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 14.0 && pt <= 16.0) * (0.897777777778) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 16.0 && pt <= 18.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 18.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 20.0 && pt <= 26.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 26.0 && pt <= 32.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 32.0 && pt <= 38.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 38.0 && pt <= 44.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 44.0 && pt <= 50.0) * (0.987153284672) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 50.0 && pt <= 60.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 60.0 && pt <= 70.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 90.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 100.0 && pt <= 125.0) * (0.75) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 150.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 10.0 && pt <= 12.0) * (0.876232201533) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 12.0 && pt <= 14.0) * (0.871212121212) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 14.0 && pt <= 16.0) * (0.982720588235) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 16.0 && pt <= 18.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 18.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 20.0 && pt <= 26.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 26.0 && pt <= 32.0) * (0.987460815047) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 32.0 && pt <= 38.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 38.0 && pt <= 44.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 44.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 50.0 && pt <= 60.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 60.0 && pt <= 70.0) * (0.850694444444) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 90.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 150.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 10.0 && pt <= 12.0) * (0.601694915254) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 12.0 && pt <= 14.0) * (0.852409638554) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 14.0 && pt <= 16.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 16.0 && pt <= 18.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 18.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 20.0 && pt <= 26.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 26.0 && pt <= 32.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 32.0 && pt <= 38.0) * (0.983394886364) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 38.0 && pt <= 44.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 44.0 && pt <= 50.0) * (0.964411190631) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 50.0 && pt <= 60.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 60.0 && pt <= 70.0) * (0.887931034483) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 90.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 150.0) * (0.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 10.0 && pt <= 12.0) * (0.670377241806) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 12.0 && pt <= 14.0) * (0.697619047619) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 14.0 && pt <= 16.0) * (0.968582375479) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 16.0 && pt <= 18.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 18.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 20.0 && pt <= 26.0) * (0.900071495353) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 26.0 && pt <= 32.0) * (0.983352627923) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 32.0 && pt <= 38.0) * (0.997842810691) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 38.0 && pt <= 44.0) * (0.974627891001) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 44.0 && pt <= 50.0) * (0.9326321516) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 50.0 && pt <= 60.0) * (0.97647702407) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 60.0 && pt <= 70.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 70.0 && pt <= 80.0) * (0.823300970874) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 90.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 150.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 12.0) * (0.638766519824) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 12.0 && pt <= 14.0) * (0.47952047952) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 14.0 && pt <= 16.0) * (0.532423208191) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 16.0 && pt <= 18.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 18.0 && pt <= 20.0) * (0.77652733119) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 26.0) * (0.72772006561) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 26.0 && pt <= 32.0) * (0.891864997098) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 32.0 && pt <= 38.0) * (0.80878477306) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 38.0 && pt <= 44.0) * (0.86904875946) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 44.0 && pt <= 50.0) * (0.92488039476) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 60.0) * (0.690144230769) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 60.0 && pt <= 70.0) * (0.83820662768) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 90.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 125.0) * (0.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 150.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 10.0 && pt <= 12.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 12.0 && pt <= 14.0) * (0.532051282051) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 14.0 && pt <= 16.0) * (0.210204081633) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 16.0 && pt <= 18.0) * (0.345794392523) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 18.0 && pt <= 20.0) * (0.772613065327) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 20.0 && pt <= 26.0) * (0.599828620394) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 26.0 && pt <= 32.0) * (0.707735247209) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 32.0 && pt <= 38.0) * (0.771581359817) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 38.0 && pt <= 44.0) * (0.842105263158) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 44.0 && pt <= 50.0) * (0.682926829268) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 50.0 && pt <= 60.0) * (0.788066723695) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 60.0 && pt <= 70.0) * (0.671604938272) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 70.0 && pt <= 80.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 100.0 && pt <= 125.0) * (0.5) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 150.0) * (0.0) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 10.0 && pt <= 12.0) * (0.0) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 12.0 && pt <= 14.0) * (0.589339339339) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 14.0 && pt <= 16.0) * (0.636666666667) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 16.0 && pt <= 18.0) * (0.875209380235) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 18.0 && pt <= 20.0) * (0.764778325123) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 20.0 && pt <= 26.0) * (0.771216318308) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 26.0 && pt <= 32.0) * (0.815113350126) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 32.0 && pt <= 38.0) * (0.801593650012) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 38.0 && pt <= 44.0) * (0.768194504079) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 44.0 && pt <= 50.0) * (0.814719945355) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 50.0 && pt <= 60.0) * (0.877339572193) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 60.0 && pt <= 70.0) * (1.0) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 100.0 && pt <= 125.0) * (0.5) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 150.0) * (0.0) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 10.0 && pt <= 12.0) * (0.321961620469) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 12.0 && pt <= 14.0) * (0.429139072848) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 14.0 && pt <= 16.0) * (0.512345679012) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 16.0 && pt <= 18.0) * (0.616042780749) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 18.0 && pt <= 20.0) * (0.551288283909) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 20.0 && pt <= 26.0) * (0.812597200622) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 26.0 && pt <= 32.0) * (0.842105263158) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 32.0 && pt <= 38.0) * (0.760898458267) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 38.0 && pt <= 44.0) * (0.716301572369) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 44.0 && pt <= 50.0) * (0.876216968011) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 50.0 && pt <= 60.0) * (0.875) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 60.0 && pt <= 70.0) * (1.0) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 70.0 && pt <= 80.0) * (0.333333333333) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 90.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 150.0) * (0.0) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 10.0 && pt <= 12.0) * (1.0) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 12.0 && pt <= 14.0) * (0.669642857143) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 14.0 && pt <= 16.0) * (1.0) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 16.0 && pt <= 18.0) * (0.528301886792) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 18.0 && pt <= 20.0) * (0.698924731183) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 20.0 && pt <= 26.0) * (0.827408256881) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 26.0 && pt <= 32.0) * (0.841000945911) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 32.0 && pt <= 38.0) * (0.815321180556) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 38.0 && pt <= 44.0) * (0.878514056225) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 44.0 && pt <= 50.0) * (0.842105263158) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 50.0 && pt <= 60.0) * (1.0) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 60.0 && pt <= 70.0) * (0.857142857143) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 100.0 && pt <= 125.0) * (0.0) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 150.0) * (0.0) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 10.0 && pt <= 12.0) * (0.710924369748) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 12.0 && pt <= 14.0) * (1.0) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 14.0 && pt <= 16.0) * (1.0) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 16.0 && pt <= 18.0) * (0.679824561404) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 18.0 && pt <= 20.0) * (0.876734693878) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 20.0 && pt <= 26.0) * (0.977759009009) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 26.0 && pt <= 32.0) * (0.8312995076) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 32.0 && pt <= 38.0) * (0.895330112721) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 38.0 && pt <= 44.0) * (0.896632471008) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 44.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 50.0 && pt <= 60.0) * (0.891747052519) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 60.0 && pt <= 70.0) * (1.0) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 100.0 && pt <= 125.0) * (0.0) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 150.0) * (0.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 10.0 && pt <= 12.0) * (0.910714285714) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 12.0 && pt <= 14.0) * (1.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 14.0 && pt <= 16.0) * (0.798) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 16.0 && pt <= 18.0) * (0.675159235669) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 18.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 20.0 && pt <= 26.0) * (0.765597920277) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 26.0 && pt <= 32.0) * (0.931067669173) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 32.0 && pt <= 38.0) * (1.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 38.0 && pt <= 44.0) * (0.973769430052) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 44.0 && pt <= 50.0) * (0.959821428571) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 50.0 && pt <= 60.0) * (1.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 60.0 && pt <= 70.0) * (1.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 150.0) * (1.0) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 10.0 && pt <= 12.0) * (0.397849462366) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 12.0 && pt <= 14.0) * (0.87476635514) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 14.0 && pt <= 16.0) * (0.887272727273) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 16.0 && pt <= 18.0) * (1.0) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 18.0 && pt <= 20.0) * (0.859173126615) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 26.0) * (0.914977870953) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 26.0 && pt <= 32.0) * (0.980363984674) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 32.0 && pt <= 38.0) * (1.0) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 38.0 && pt <= 44.0) * (1.0) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 44.0 && pt <= 50.0) * (0.91966721222) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 60.0) * (0.803791469194) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 60.0 && pt <= 70.0) * (1.0) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 80.0 && pt <= 90.0) * (0.0) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 125.0) * (0.0) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 150.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 10.0 && pt <= 12.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 12.0 && pt <= 14.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 14.0 && pt <= 16.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 16.0 && pt <= 18.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 18.0 && pt <= 20.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 20.0 && pt <= 26.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 26.0 && pt <= 32.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 32.0 && pt <= 38.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 38.0 && pt <= 44.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 44.0 && pt <= 50.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 50.0 && pt <= 60.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 60.0 && pt <= 70.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 70.0 && pt <= 80.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 80.0 && pt <= 90.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 100.0 && pt <= 125.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 150.0) * (0.0)  
  }


}

#######################
# Electron medium ID efficiency #
#######################

##FIXME!!! sourcing LooseId tcl file because medium does not exists (yet ...)
module Efficiency ElectronMediumEfficiency {

  set InputArray ElectronSmear/electrons
  set OutputArray electrons

    ### electron medium ID 

  set EfficiencyFormula {

   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 10.0 && pt <= 12.0) * (0.708722741433) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 12.0 && pt <= 14.0) * (0.856275303644) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 14.0 && pt <= 16.0) * (0.471595330739) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 16.0 && pt <= 18.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 18.0 && pt <= 20.0) * (0.954861111111) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 20.0 && pt <= 26.0) * (0.802285385059) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 26.0 && pt <= 32.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 32.0 && pt <= 38.0) * (0.910882594158) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 38.0 && pt <= 44.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 44.0 && pt <= 50.0) * (0.992587155963) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 50.0 && pt <= 60.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 60.0 && pt <= 70.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 90.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 100.0 && pt <= 125.0) * (0.75) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 150.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 10.0 && pt <= 12.0) * (0.561167227834) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 12.0 && pt <= 14.0) * (0.530769230769) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 14.0 && pt <= 16.0) * (0.896603773585) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 16.0 && pt <= 18.0) * (0.862359550562) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 18.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 20.0 && pt <= 26.0) * (0.730981887512) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 26.0 && pt <= 32.0) * (0.916932907348) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 32.0 && pt <= 38.0) * (0.950977189831) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 38.0 && pt <= 44.0) * (0.975859339665) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 44.0 && pt <= 50.0) * (0.978531672409) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 50.0 && pt <= 60.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 60.0 && pt <= 70.0) * (0.846286701209) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 90.0 && pt <= 100.0) * (0.666666666667) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 150.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 10.0 && pt <= 12.0) * (0.620087336245) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 12.0 && pt <= 14.0) * (0.730888429752) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 14.0 && pt <= 16.0) * (0.774269005848) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 16.0 && pt <= 18.0) * (0.92027972028) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 18.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 20.0 && pt <= 26.0) * (0.8910933082) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 26.0 && pt <= 32.0) * (0.969176503597) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 32.0 && pt <= 38.0) * (0.887604017217) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 38.0 && pt <= 44.0) * (0.995904365177) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 44.0 && pt <= 50.0) * (0.864595295121) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 50.0 && pt <= 60.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 60.0 && pt <= 70.0) * (0.892326732673) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 90.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 150.0) * (0.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 10.0 && pt <= 12.0) * (0.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 12.0 && pt <= 14.0) * (0.711740890688) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 14.0 && pt <= 16.0) * (0.613830613831) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 16.0 && pt <= 18.0) * (0.807569296375) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 18.0 && pt <= 20.0) * (0.912071535022) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 20.0 && pt <= 26.0) * (0.728348909657) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 26.0 && pt <= 32.0) * (0.820519480519) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 32.0 && pt <= 38.0) * (0.989468599034) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 38.0 && pt <= 44.0) * (0.947283143337) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 44.0 && pt <= 50.0) * (0.887996531977) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 50.0 && pt <= 60.0) * (0.853752759382) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 60.0 && pt <= 70.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 70.0 && pt <= 80.0) * (0.617475728155) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 80.0 && pt <= 90.0) * (0.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 90.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 150.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 12.0) * (0.263636363636) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 12.0 && pt <= 14.0) * (0.29417773238) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 14.0 && pt <= 16.0) * (0.543554006969) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 16.0 && pt <= 18.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 18.0 && pt <= 20.0) * (0.703825136612) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 26.0) * (0.676208513148) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 26.0 && pt <= 32.0) * (0.840200241865) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 32.0 && pt <= 38.0) * (0.737079318013) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 38.0 && pt <= 44.0) * (0.841153846154) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 44.0 && pt <= 50.0) * (0.906075392511) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 60.0) * (0.63226744186) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 60.0 && pt <= 70.0) * (0.502923976608) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 90.0 && pt <= 100.0) * (0.5) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 125.0) * (0.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 150.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 10.0 && pt <= 12.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 12.0 && pt <= 14.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 14.0 && pt <= 16.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 16.0 && pt <= 18.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 18.0 && pt <= 20.0) * (0.792525773196) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 20.0 && pt <= 26.0) * (0.289729103288) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 26.0 && pt <= 32.0) * (0.573042776433) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 32.0 && pt <= 38.0) * (0.660798191875) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 38.0 && pt <= 44.0) * (0.58358306799) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 44.0 && pt <= 50.0) * (0.614706755999) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 50.0 && pt <= 60.0) * (0.718575718576) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 60.0 && pt <= 70.0) * (0.335802469136) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 70.0 && pt <= 80.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 100.0 && pt <= 125.0) * (0.5) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 150.0) * (0.0) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 10.0 && pt <= 12.0) * (0.0) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 12.0 && pt <= 14.0) * (0.358447488584) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 14.0 && pt <= 16.0) * (0.43908045977) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 16.0 && pt <= 18.0) * (0.897766323024) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 18.0 && pt <= 20.0) * (0.646875) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 20.0 && pt <= 26.0) * (0.661566156616) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 26.0 && pt <= 32.0) * (0.661251596424) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 32.0 && pt <= 38.0) * (0.69225251076) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 38.0 && pt <= 44.0) * (0.594717225733) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 44.0 && pt <= 50.0) * (0.818072702332) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 50.0 && pt <= 60.0) * (0.884433962264) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 60.0 && pt <= 70.0) * (1.0) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 80.0 && pt <= 90.0) * (0.0) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 100.0 && pt <= 125.0) * (0.5) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 150.0) * (0.0) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 10.0 && pt <= 12.0) * (0.169853768279) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 12.0 && pt <= 14.0) * (0.45) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 14.0 && pt <= 16.0) * (0.350210970464) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 16.0 && pt <= 18.0) * (0.533333333333) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 18.0 && pt <= 20.0) * (0.466897233202) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 20.0 && pt <= 26.0) * (0.709512578616) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 26.0 && pt <= 32.0) * (0.850632911392) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 32.0 && pt <= 38.0) * (0.704662519232) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 38.0 && pt <= 44.0) * (0.599172624668) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 44.0 && pt <= 50.0) * (0.757363253857) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 50.0 && pt <= 60.0) * (0.880537974684) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 60.0 && pt <= 70.0) * (1.0) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 70.0 && pt <= 80.0) * (0.333333333333) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 90.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 150.0) * (0.0) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 10.0 && pt <= 12.0) * (1.0) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 12.0 && pt <= 14.0) * (0.552147239264) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 14.0 && pt <= 16.0) * (0.917933130699) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 16.0 && pt <= 18.0) * (0.535031847134) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 18.0 && pt <= 20.0) * (0.698924731183) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 20.0 && pt <= 26.0) * (0.776827371695) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 26.0 && pt <= 32.0) * (0.782627673505) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 32.0 && pt <= 38.0) * (0.824868266979) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 38.0 && pt <= 44.0) * (0.786757301108) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 44.0 && pt <= 50.0) * (0.846770666278) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 50.0 && pt <= 60.0) * (1.0) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 60.0 && pt <= 70.0) * (0.714285714286) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 100.0 && pt <= 125.0) * (0.0) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 150.0) * (0.0) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 10.0 && pt <= 12.0) * (0.482051282051) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 12.0 && pt <= 14.0) * (0.651282051282) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 14.0 && pt <= 16.0) * (1.0) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 16.0 && pt <= 18.0) * (0.35632183908) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 18.0 && pt <= 20.0) * (0.892026578073) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 20.0 && pt <= 26.0) * (0.909638554217) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 26.0 && pt <= 32.0) * (0.802794022092) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 32.0 && pt <= 38.0) * (0.86745353412) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 38.0 && pt <= 44.0) * (0.902367859948) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 44.0 && pt <= 50.0) * (0.971212121212) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 50.0 && pt <= 60.0) * (0.780278670954) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 60.0 && pt <= 70.0) * (1.0) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 100.0 && pt <= 125.0) * (0.0) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 150.0) * (0.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 10.0 && pt <= 12.0) * (0.316489361702) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 12.0 && pt <= 14.0) * (0.540540540541) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 14.0 && pt <= 16.0) * (0.54958677686) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 16.0 && pt <= 18.0) * (0.716216216216) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 18.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 20.0 && pt <= 26.0) * (0.784635879218) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 26.0 && pt <= 32.0) * (0.95549382716) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 32.0 && pt <= 38.0) * (1.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 38.0 && pt <= 44.0) * (0.983965968586) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 44.0 && pt <= 50.0) * (0.958333333333) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 50.0 && pt <= 60.0) * (1.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 60.0 && pt <= 70.0) * (1.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 150.0) * (1.0) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 10.0 && pt <= 12.0) * (0.203296703297) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 12.0 && pt <= 14.0) * (0.668571428571) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 14.0 && pt <= 16.0) * (0.228037383178) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 16.0 && pt <= 18.0) * (0.861328125) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 18.0 && pt <= 20.0) * (0.908469945355) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 26.0) * (0.86409087273) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 26.0 && pt <= 32.0) * (1.0) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 32.0 && pt <= 38.0) * (1.0) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 38.0 && pt <= 44.0) * (1.0) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 44.0 && pt <= 50.0) * (0.933416389812) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 60.0) * (0.803791469194) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 60.0 && pt <= 70.0) * (1.0) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 80.0 && pt <= 90.0) * (0.0) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 125.0) * (0.0) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 150.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 10.0 && pt <= 12.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 12.0 && pt <= 14.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 14.0 && pt <= 16.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 16.0 && pt <= 18.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 18.0 && pt <= 20.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 20.0 && pt <= 26.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 26.0 && pt <= 32.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 32.0 && pt <= 38.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 38.0 && pt <= 44.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 44.0 && pt <= 50.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 50.0 && pt <= 60.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 60.0 && pt <= 70.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 70.0 && pt <= 80.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 80.0 && pt <= 90.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 100.0 && pt <= 125.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 150.0) * (0.0)  
  }
}

#######################
# Electron tight ID efficiency #
#######################

module Efficiency ElectronTightEfficiency {

  set InputArray ElectronSmear/electrons
  set OutputArray electrons

    ### electron tight ID 

  set EfficiencyFormula {

   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 10.0 && pt <= 12.0) * (0.446078431373) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 12.0 && pt <= 14.0) * (0.30652173913) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 14.0 && pt <= 16.0) * (0.26008583691) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 16.0 && pt <= 18.0) * (0.315134099617) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 18.0 && pt <= 20.0) * (0.340768277571) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 20.0 && pt <= 26.0) * (0.311331300813) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 26.0 && pt <= 32.0) * (0.520345475571) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 32.0 && pt <= 38.0) * (0.642006269592) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 38.0 && pt <= 44.0) * (0.693425641026) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 44.0 && pt <= 50.0) * (0.727096774194) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 50.0 && pt <= 60.0) * (0.6) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 60.0 && pt <= 70.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 80.0 && pt <= 90.0) * (0.377777777778) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 90.0 && pt <= 100.0) * (0.548387096774) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 100.0 && pt <= 125.0) * (0.75) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 150.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 10.0 && pt <= 12.0) * (0.117554858934) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 12.0 && pt <= 14.0) * (0.0927419354839) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 14.0 && pt <= 16.0) * (0.365163934426) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 16.0 && pt <= 18.0) * (0.303359683794) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 18.0 && pt <= 20.0) * (0.339948783611) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 20.0 && pt <= 26.0) * (0.40947476828) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 26.0 && pt <= 32.0) * (0.35205364627) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 32.0 && pt <= 38.0) * (0.592125779626) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 38.0 && pt <= 44.0) * (0.699958434772) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 44.0 && pt <= 50.0) * (0.707641196013) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 50.0 && pt <= 60.0) * (0.627314814815) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 60.0 && pt <= 70.0) * (0.687719298246) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 70.0 && pt <= 80.0) * (0.856470588235) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 90.0 && pt <= 100.0) * (0.686868686869) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 100.0 && pt <= 125.0) * (0.526315789474) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 150.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 10.0 && pt <= 12.0) * (0.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 12.0 && pt <= 14.0) * (0.155837004405) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 14.0 && pt <= 16.0) * (0.404151404151) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 16.0 && pt <= 18.0) * (0.247368421053) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 18.0 && pt <= 20.0) * (0.403823178017) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 20.0 && pt <= 26.0) * (0.404334677419) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 26.0 && pt <= 32.0) * (0.516794306168) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 32.0 && pt <= 38.0) * (0.555011303693) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 38.0 && pt <= 44.0) * (0.728677537779) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 44.0 && pt <= 50.0) * (0.626733371386) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 50.0 && pt <= 60.0) * (0.85235277543) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 60.0 && pt <= 70.0) * (0.520202020202) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 70.0 && pt <= 80.0) * (0.875776397516) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 100.0 && pt <= 125.0) * (0.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 150.0) * (0.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 10.0 && pt <= 12.0) * (0.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 12.0 && pt <= 14.0) * (0.500854700855) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 14.0 && pt <= 16.0) * (0.39012345679) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 16.0 && pt <= 18.0) * (0.168426903835) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 18.0 && pt <= 20.0) * (0.433811802233) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 20.0 && pt <= 26.0) * (0.207361419069) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 26.0 && pt <= 32.0) * (0.382371874226) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 32.0 && pt <= 38.0) * (0.619023440806) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 38.0 && pt <= 44.0) * (0.473377914332) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 44.0 && pt <= 50.0) * (0.629523045714) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 50.0 && pt <= 60.0) * (0.407534246575) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 60.0 && pt <= 70.0) * (0.659298780488) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 70.0 && pt <= 80.0) * (0.209900990099) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 80.0 && pt <= 90.0) * (0.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 90.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 150.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 12.0) * (0.13679245283) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 12.0 && pt <= 14.0) * (0.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 14.0 && pt <= 16.0) * (0.293233082707) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 16.0 && pt <= 18.0) * (0.418269230769) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 18.0 && pt <= 20.0) * (0.379269729093) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 26.0) * (0.249901652242) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 26.0 && pt <= 32.0) * (0.561725130736) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 32.0 && pt <= 38.0) * (0.452786885246) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 38.0 && pt <= 44.0) * (0.712087912088) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 44.0 && pt <= 50.0) * (0.687251153593) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 60.0) * (0.450443786982) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 60.0 && pt <= 70.0) * (0.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 125.0) * (0.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 150.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 10.0 && pt <= 12.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 12.0 && pt <= 14.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 14.0 && pt <= 16.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 16.0 && pt <= 18.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 18.0 && pt <= 20.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 20.0 && pt <= 26.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 26.0 && pt <= 32.0) * (0.121523320496) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 32.0 && pt <= 38.0) * (0.156182346756) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 38.0 && pt <= 44.0) * (0.137273432199) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 44.0 && pt <= 50.0) * (0.126185636856) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 50.0 && pt <= 60.0) * (0.293602103418) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 60.0 && pt <= 70.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 70.0 && pt <= 80.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 80.0 && pt <= 90.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 100.0 && pt <= 125.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 150.0) * (0.0) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 10.0 && pt <= 12.0) * (0.0) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 12.0 && pt <= 14.0) * (0.123719464145) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 14.0 && pt <= 16.0) * (0.234355828221) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 16.0 && pt <= 18.0) * (0.0) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 18.0 && pt <= 20.0) * (0.134765625) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 20.0 && pt <= 26.0) * (0.130249867092) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 26.0 && pt <= 32.0) * (0.259225634179) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 32.0 && pt <= 38.0) * (0.31121001032) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 38.0 && pt <= 44.0) * (0.244625895684) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 44.0 && pt <= 50.0) * (0.128142458101) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 50.0 && pt <= 60.0) * (0.128777472527) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 60.0 && pt <= 70.0) * (0.507936507937) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 70.0 && pt <= 80.0) * (0.688888888889) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 80.0 && pt <= 90.0) * (0.0) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 100.0 && pt <= 125.0) * (0.0) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 150.0) * (0.0) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 10.0 && pt <= 12.0) * (0.192602040816) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 12.0 && pt <= 14.0) * (0.117391304348) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 14.0 && pt <= 16.0) * (0.0) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 16.0 && pt <= 18.0) * (0.219428571429) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 18.0 && pt <= 20.0) * (0.0970724191063) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 20.0 && pt <= 26.0) * (0.313014827018) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 26.0 && pt <= 32.0) * (0.259946949602) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 32.0 && pt <= 38.0) * (0.125742132537) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 38.0 && pt <= 44.0) * (0.296853244222) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 44.0 && pt <= 50.0) * (0.214899713467) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 50.0 && pt <= 60.0) * (0.511254019293) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 60.0 && pt <= 70.0) * (0.205504587156) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 70.0 && pt <= 80.0) * (0.0) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 80.0 && pt <= 90.0) * (0.0) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 100.0 && pt <= 125.0) * (0.0) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 150.0) * (0.0) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 10.0 && pt <= 12.0) * (0.564356435644) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 12.0 && pt <= 14.0) * (0.288461538462) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 14.0 && pt <= 16.0) * (0.319576719577) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 16.0 && pt <= 18.0) * (0.0) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 18.0 && pt <= 20.0) * (0.0) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 20.0 && pt <= 26.0) * (0.548599670511) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 26.0 && pt <= 32.0) * (0.408912607188) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 32.0 && pt <= 38.0) * (0.562385531136) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 38.0 && pt <= 44.0) * (0.457635983264) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 44.0 && pt <= 50.0) * (0.597039473684) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 50.0 && pt <= 60.0) * (0.722371967655) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 60.0 && pt <= 70.0) * (0.577380952381) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 70.0 && pt <= 80.0) * (0.5) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 100.0 && pt <= 125.0) * (0.0) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 150.0) * (0.0) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 10.0 && pt <= 12.0) * (0.0) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 12.0 && pt <= 14.0) * (0.226785714286) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 14.0 && pt <= 16.0) * (0.486857142857) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 16.0 && pt <= 18.0) * (0.0) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 18.0 && pt <= 20.0) * (0.470639789658) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 20.0 && pt <= 26.0) * (0.273550724638) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 26.0 && pt <= 32.0) * (0.478319783198) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 32.0 && pt <= 38.0) * (0.66638465877) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 38.0 && pt <= 44.0) * (0.76431209603) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 44.0 && pt <= 50.0) * (0.723709677419) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 50.0 && pt <= 60.0) * (0.560949298813) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 60.0 && pt <= 70.0) * (0.833333333333) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 70.0 && pt <= 80.0) * (0.5) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 100.0 && pt <= 125.0) * (0.0) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 150.0) * (0.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 10.0 && pt <= 12.0) * (0.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 12.0 && pt <= 14.0) * (0.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 14.0 && pt <= 16.0) * (0.279411764706) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 16.0 && pt <= 18.0) * (0.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 18.0 && pt <= 20.0) * (0.799798792757) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 20.0 && pt <= 26.0) * (0.633917589176) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 26.0 && pt <= 32.0) * (0.566278317152) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 32.0 && pt <= 38.0) * (0.733150730412) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 38.0 && pt <= 44.0) * (0.719390507012) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 44.0 && pt <= 50.0) * (0.68145800317) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 50.0 && pt <= 60.0) * (0.916239316239) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 60.0 && pt <= 70.0) * (1.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 70.0 && pt <= 80.0) * (0.333333333333) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 150.0) * (1.0) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 10.0 && pt <= 12.0) * (0.0) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 12.0 && pt <= 14.0) * (0.238775510204) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 14.0 && pt <= 16.0) * (0.0) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 16.0 && pt <= 18.0) * (0.294) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 18.0 && pt <= 20.0) * (0.568376068376) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 26.0) * (0.372817008352) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 26.0 && pt <= 32.0) * (0.796555354994) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 32.0 && pt <= 38.0) * (0.802384500745) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 38.0 && pt <= 44.0) * (0.709594303382) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 44.0 && pt <= 50.0) * (0.812279151943) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 60.0) * (0.83137254902) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 60.0 && pt <= 70.0) * (1.0) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 80.0 && pt <= 90.0) * (0.0) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 125.0) * (0.0) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 150.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 10.0 && pt <= 12.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 12.0 && pt <= 14.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 14.0 && pt <= 16.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 16.0 && pt <= 18.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 18.0 && pt <= 20.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 20.0 && pt <= 26.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 26.0 && pt <= 32.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 32.0 && pt <= 38.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 38.0 && pt <= 44.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 44.0 && pt <= 50.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 50.0 && pt <= 60.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 60.0 && pt <= 70.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 70.0 && pt <= 80.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 80.0 && pt <= 90.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 100.0 && pt <= 125.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 3.0) * (pt > 150.0) * (0.0)  
  }
}


##################
# Muon isolation #
##################

module Isolation MuonIsolation {
  set CandidateInputArray MuonMomentumSmearing/muons

  # isolation collection
  set IsolationInputArray EFlowFilterPuppi/eflow

  set OutputArray muons

  set DeltaRMax 0.3
  set PTMin 0.0
  set PTRatioMax 9999.

}


##############
# Muon scale #
##############


module EnergyScale MuonScale {
  set InputArray MuonIsolation/muons
  set OutputArray muons

    ### muon looseIDISO momentum scale
  set ScaleFormula {
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 20.0) * (1.106) +
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 50.0) * (1.103) +
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 100.0) * (1.000) +
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 150.0) * (0.985) +
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 150.0 && pt <= 14000.0) * (0.985) +
   (abs(eta) > 1.5 && abs(eta) <= 2.8) * (pt > 10.0 && pt <= 20.0) * (0.996) +
   (abs(eta) > 1.5 && abs(eta) <= 2.8) * (pt > 20.0 && pt <= 50.0) * (1.002) +
   (abs(eta) > 1.5 && abs(eta) <= 2.8) * (pt > 50.0 && pt <= 100.0) * (0.996) +
   (abs(eta) > 1.5 && abs(eta) <= 2.8) * (pt > 100.0 && pt <= 150.0) * (0.992) +
   (abs(eta) > 1.5 && abs(eta) <= 2.8) * (pt > 150.0 && pt <= 14000.0) * (0.132)  
  }

}


##############
# Muon smear #
##############

module MomentumSmearing MuonSmear {

  set InputArray MuonScale/muons
  set OutputArray muons

    ### muon looseIDISO momentum resolution
  set ResolutionFormula {
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 20.0) * (0.000001) +
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 50.0) * (0.000001) +
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 100.0) * (0.000019) +
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 150.0) * (0.050900) +
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 150.0 && pt <= 14000.0) * (0.050900) +
   (abs(eta) > 1.5 && abs(eta) <= 2.8) * (pt > 10.0 && pt <= 20.0) * (0.000001) +
   (abs(eta) > 1.5 && abs(eta) <= 2.8) * (pt > 20.0 && pt <= 50.0) * (0.002806) +
   (abs(eta) > 1.5 && abs(eta) <= 2.8) * (pt > 50.0 && pt <= 100.0) * (0.000001) +
   (abs(eta) > 1.5 && abs(eta) <= 2.8) * (pt > 100.0 && pt <= 150.0) * (0.023766) +
   (abs(eta) > 1.5 && abs(eta) <= 2.8) * (pt > 150.0 && pt <= 14000.0) * (0.197689)  
  }

}


##################
# Muon Loose Id  #
##################

module Efficiency MuonLooseIdEfficiency {
    set InputArray MuonSmear/muons
    set OutputArray muons
    
      ### muon loose ID 

  set EfficiencyFormula {

   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 10.0 && pt <= 12.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 12.0 && pt <= 14.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 14.0 && pt <= 16.0) * (0.97829506134) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 16.0 && pt <= 18.0) * (0.984410937899) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 18.0 && pt <= 20.0) * (0.880836236934) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 20.0 && pt <= 26.0) * (0.941622609813) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 26.0 && pt <= 32.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 32.0 && pt <= 38.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 38.0 && pt <= 44.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 44.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 50.0 && pt <= 60.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 60.0 && pt <= 70.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 70.0 && pt <= 80.0) * (0.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 90.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 150.0) * (0.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 10.0 && pt <= 12.0) * (0.808015513898) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 12.0 && pt <= 14.0) * (0.951359084406) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 14.0 && pt <= 16.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 16.0 && pt <= 18.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 18.0 && pt <= 20.0) * (0.831893687708) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 20.0 && pt <= 26.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 26.0 && pt <= 32.0) * (0.967479350385) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 32.0 && pt <= 38.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 38.0 && pt <= 44.0) * (0.990420560748) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 44.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 50.0 && pt <= 60.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 60.0 && pt <= 70.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 100.0 && pt <= 125.0) * (0.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 150.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 10.0 && pt <= 12.0) * (0.824647455549) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 12.0 && pt <= 14.0) * (0.929171668667) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 14.0 && pt <= 16.0) * (0.970889708897) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 16.0 && pt <= 18.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 18.0 && pt <= 20.0) * (0.96130952381) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 20.0 && pt <= 26.0) * (0.961125890329) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 26.0 && pt <= 32.0) * (0.995452324576) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 32.0 && pt <= 38.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 38.0 && pt <= 44.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 44.0 && pt <= 50.0) * (0.980914538822) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 50.0 && pt <= 60.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 60.0 && pt <= 70.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 90.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 150.0) * (0.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 10.0 && pt <= 12.0) * (0.789049919485) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 12.0 && pt <= 14.0) * (0.846581196581) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 14.0 && pt <= 16.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 16.0 && pt <= 18.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 18.0 && pt <= 20.0) * (0.899297423888) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 20.0 && pt <= 26.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 26.0 && pt <= 32.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 32.0 && pt <= 38.0) * (0.994935770751) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 38.0 && pt <= 44.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 44.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 50.0 && pt <= 60.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 60.0 && pt <= 70.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 150.0) * (0.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 12.0) * (0.751461988304) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 12.0 && pt <= 14.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 14.0 && pt <= 16.0) * (0.825745682889) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 16.0 && pt <= 18.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 18.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 26.0) * (0.982509505703) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 26.0 && pt <= 32.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 32.0 && pt <= 38.0) * (0.992863576526) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 38.0 && pt <= 44.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 44.0 && pt <= 50.0) * (0.975572188395) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 60.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 60.0 && pt <= 70.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 90.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 150.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 10.0 && pt <= 12.0) * (0.497554157932) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 12.0 && pt <= 14.0) * (0.730837789661) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 14.0 && pt <= 16.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 16.0 && pt <= 18.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 18.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 20.0 && pt <= 26.0) * (0.944187836798) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 26.0 && pt <= 32.0) * (0.963597328244) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 32.0 && pt <= 38.0) * (0.982492209141) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 38.0 && pt <= 44.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 44.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 50.0 && pt <= 60.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 60.0 && pt <= 70.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 80.0 && pt <= 90.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 150.0) * (1.0) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 10.0 && pt <= 12.0) * (0.925925925926) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 12.0 && pt <= 14.0) * (1.0) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 14.0 && pt <= 16.0) * (1.0) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 16.0 && pt <= 18.0) * (0.865231259968) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 18.0 && pt <= 20.0) * (0.848348348348) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 20.0 && pt <= 26.0) * (0.975296442688) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 26.0 && pt <= 32.0) * (0.983988355167) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 32.0 && pt <= 38.0) * (0.983784006256) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 38.0 && pt <= 44.0) * (0.973464052288) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 44.0 && pt <= 50.0) * (0.976013513514) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 50.0 && pt <= 60.0) * (1.0) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 60.0 && pt <= 70.0) * (1.0) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 90.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 150.0) * (0.0) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 10.0 && pt <= 12.0) * (0.448484848485) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 12.0 && pt <= 14.0) * (1.0) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 14.0 && pt <= 16.0) * (0.971704623879) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 16.0 && pt <= 18.0) * (1.0) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 18.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 20.0 && pt <= 26.0) * (1.0) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 26.0 && pt <= 32.0) * (1.0) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 32.0 && pt <= 38.0) * (1.0) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 38.0 && pt <= 44.0) * (1.0) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 44.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 50.0 && pt <= 60.0) * (1.0) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 60.0 && pt <= 70.0) * (1.0) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 150.0) * (0.0) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 10.0 && pt <= 12.0) * (1.0) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 12.0 && pt <= 14.0) * (1.0) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 14.0 && pt <= 16.0) * (1.0) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 16.0 && pt <= 18.0) * (1.0) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 18.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 20.0 && pt <= 26.0) * (0.984347136913) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 26.0 && pt <= 32.0) * (0.928483897848) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 32.0 && pt <= 38.0) * (0.985714285714) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 38.0 && pt <= 44.0) * (0.980057991597) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 44.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 50.0 && pt <= 60.0) * (0.891542288557) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 60.0 && pt <= 70.0) * (0.8) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 150.0) * (0.0) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 10.0 && pt <= 12.0) * (0.946788990826) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 12.0 && pt <= 14.0) * (1.0) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 14.0 && pt <= 16.0) * (0.820987654321) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 16.0 && pt <= 18.0) * (1.0) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 18.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 20.0 && pt <= 26.0) * (1.0) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 26.0 && pt <= 32.0) * (1.0) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 32.0 && pt <= 38.0) * (1.0) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 38.0 && pt <= 44.0) * (1.0) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 44.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 50.0 && pt <= 60.0) * (1.0) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 60.0 && pt <= 70.0) * (1.0) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 90.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 150.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 10.0 && pt <= 12.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 12.0 && pt <= 14.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 14.0 && pt <= 16.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 16.0 && pt <= 18.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 18.0 && pt <= 20.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 20.0 && pt <= 26.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 26.0 && pt <= 32.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 32.0 && pt <= 38.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 38.0 && pt <= 44.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 44.0 && pt <= 50.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 50.0 && pt <= 60.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 60.0 && pt <= 70.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 70.0 && pt <= 80.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 80.0 && pt <= 90.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 100.0 && pt <= 125.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 150.0) * (0.0)  
  }
}

##################
# Muon Medium Id  #
##################

##FIXME!!! sourcing LooseId tcl file because medium does not exists (yet ...)
module Efficiency MuonMediumIdEfficiency {
    set InputArray MuonSmear/muons
    set OutputArray muons

      ### muon medium ID 

  set EfficiencyFormula {

   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 10.0 && pt <= 12.0) * (0.864285714286) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 12.0 && pt <= 14.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 14.0 && pt <= 16.0) * (0.912023460411) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 16.0 && pt <= 18.0) * (0.927014964558) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 18.0 && pt <= 20.0) * (0.89014084507) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 20.0 && pt <= 26.0) * (0.865514312834) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 26.0 && pt <= 32.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 32.0 && pt <= 38.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 38.0 && pt <= 44.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 44.0 && pt <= 50.0) * (0.986882868631) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 50.0 && pt <= 60.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 60.0 && pt <= 70.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 70.0 && pt <= 80.0) * (0.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 90.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 150.0) * (0.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 10.0 && pt <= 12.0) * (0.834445927904) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 12.0 && pt <= 14.0) * (0.972222222222) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 14.0 && pt <= 16.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 16.0 && pt <= 18.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 18.0 && pt <= 20.0) * (0.840268456376) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 20.0 && pt <= 26.0) * (0.987310246679) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 26.0 && pt <= 32.0) * (0.979359373107) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 32.0 && pt <= 38.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 38.0 && pt <= 44.0) * (0.997998822837) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 44.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 50.0 && pt <= 60.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 60.0 && pt <= 70.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 100.0 && pt <= 125.0) * (0.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 150.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 10.0 && pt <= 12.0) * (0.677155443675) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 12.0 && pt <= 14.0) * (0.784194528875) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 14.0 && pt <= 16.0) * (0.989139515455) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 16.0 && pt <= 18.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 18.0 && pt <= 20.0) * (0.977172717272) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 20.0 && pt <= 26.0) * (0.952028072623) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 26.0 && pt <= 32.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 32.0 && pt <= 38.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 38.0 && pt <= 44.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 44.0 && pt <= 50.0) * (0.987432243399) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 50.0 && pt <= 60.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 60.0 && pt <= 70.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 90.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 150.0) * (0.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 10.0 && pt <= 12.0) * (0.680555555556) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 12.0 && pt <= 14.0) * (0.736979166667) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 14.0 && pt <= 16.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 16.0 && pt <= 18.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 18.0 && pt <= 20.0) * (0.920421860019) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 20.0 && pt <= 26.0) * (0.994632535095) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 26.0 && pt <= 32.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 32.0 && pt <= 38.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 38.0 && pt <= 44.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 44.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 50.0 && pt <= 60.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 60.0 && pt <= 70.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 150.0) * (0.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 12.0) * (0.634567901235) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 12.0 && pt <= 14.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 14.0 && pt <= 16.0) * (0.835982199619) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 16.0 && pt <= 18.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 18.0 && pt <= 20.0) * (0.913165266106) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 26.0) * (0.966762728146) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 26.0 && pt <= 32.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 32.0 && pt <= 38.0) * (0.997799688057) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 38.0 && pt <= 44.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 44.0 && pt <= 50.0) * (0.979935212493) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 60.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 60.0 && pt <= 70.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 90.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 150.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 10.0 && pt <= 12.0) * (0.50070323488) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 12.0 && pt <= 14.0) * (0.742753623188) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 14.0 && pt <= 16.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 16.0 && pt <= 18.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 18.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 20.0 && pt <= 26.0) * (0.948569218871) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 26.0 && pt <= 32.0) * (0.973818707811) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 32.0 && pt <= 38.0) * (0.962269083304) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 38.0 && pt <= 44.0) * (0.988433586399) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 44.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 50.0 && pt <= 60.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 60.0 && pt <= 70.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 80.0 && pt <= 90.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 150.0) * (1.0) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 10.0 && pt <= 12.0) * (0.938189845475) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 12.0 && pt <= 14.0) * (1.0) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 14.0 && pt <= 16.0) * (1.0) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 16.0 && pt <= 18.0) * (0.873590982287) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 18.0 && pt <= 20.0) * (0.8599695586) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 20.0 && pt <= 26.0) * (0.9836120031) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 26.0 && pt <= 32.0) * (0.987221613728) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 32.0 && pt <= 38.0) * (0.988053902116) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 38.0 && pt <= 44.0) * (0.946067513865) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 44.0 && pt <= 50.0) * (0.979073117004) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 50.0 && pt <= 60.0) * (1.0) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 60.0 && pt <= 70.0) * (1.0) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 90.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 150.0) * (0.0) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 10.0 && pt <= 12.0) * (0.226993865031) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 12.0 && pt <= 14.0) * (1.0) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 14.0 && pt <= 16.0) * (0.98392732355) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 16.0 && pt <= 18.0) * (1.0) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 18.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 20.0 && pt <= 26.0) * (0.975453438441) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 26.0 && pt <= 32.0) * (1.0) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 32.0 && pt <= 38.0) * (1.0) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 38.0 && pt <= 44.0) * (1.0) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 44.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 50.0 && pt <= 60.0) * (1.0) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 60.0 && pt <= 70.0) * (1.0) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 150.0) * (0.0) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 10.0 && pt <= 12.0) * (0.963779527559) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 12.0 && pt <= 14.0) * (0.718875502008) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 14.0 && pt <= 16.0) * (0.832820512821) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 16.0 && pt <= 18.0) * (0.882059800664) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 18.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 20.0 && pt <= 26.0) * (0.962575452716) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 26.0 && pt <= 32.0) * (0.931634861302) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 32.0 && pt <= 38.0) * (0.987789473684) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 38.0 && pt <= 44.0) * (0.985071075953) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 44.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 50.0 && pt <= 60.0) * (0.896896896897) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 60.0 && pt <= 70.0) * (0.806779661017) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 150.0) * (0.0) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 10.0 && pt <= 12.0) * (0.955555555556) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 12.0 && pt <= 14.0) * (1.0) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 14.0 && pt <= 16.0) * (0.847133757962) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 16.0 && pt <= 18.0) * (1.0) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 18.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 20.0 && pt <= 26.0) * (1.0) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 26.0 && pt <= 32.0) * (1.0) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 32.0 && pt <= 38.0) * (1.0) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 38.0 && pt <= 44.0) * (1.0) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 44.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 50.0 && pt <= 60.0) * (1.0) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 60.0 && pt <= 70.0) * (1.0) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 90.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 150.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 10.0 && pt <= 12.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 12.0 && pt <= 14.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 14.0 && pt <= 16.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 16.0 && pt <= 18.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 18.0 && pt <= 20.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 20.0 && pt <= 26.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 26.0 && pt <= 32.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 32.0 && pt <= 38.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 38.0 && pt <= 44.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 44.0 && pt <= 50.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 50.0 && pt <= 60.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 60.0 && pt <= 70.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 70.0 && pt <= 80.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 80.0 && pt <= 90.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 100.0 && pt <= 125.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 150.0) * (0.0)  
  }
}

##################
# Muon Tight Id  #
##################

module Efficiency MuonTightIdEfficiency {
    set InputArray MuonSmear/muons
    set OutputArray muons

      ### muon tight ID 

  set EfficiencyFormula {

   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 10.0 && pt <= 12.0) * (0.354578754579) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 12.0 && pt <= 14.0) * (0.231275720165) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 14.0 && pt <= 16.0) * (0.414252414252) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 16.0 && pt <= 18.0) * (0.525368248773) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 18.0 && pt <= 20.0) * (0.342238267148) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 20.0 && pt <= 26.0) * (0.824579296277) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 26.0 && pt <= 32.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 32.0 && pt <= 38.0) * (0.861548291633) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 38.0 && pt <= 44.0) * (0.940942237574) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 44.0 && pt <= 50.0) * (0.868943393704) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 50.0 && pt <= 60.0) * (0.950746268657) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 60.0 && pt <= 70.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 70.0 && pt <= 80.0) * (0.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 90.0 && pt <= 100.0) * (0.5) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 150.0) * (0.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 10.0 && pt <= 12.0) * (0.341763499658) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 12.0 && pt <= 14.0) * (0.588495575221) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 14.0 && pt <= 16.0) * (0.319520174482) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 16.0 && pt <= 18.0) * (0.418269230769) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 18.0 && pt <= 20.0) * (0.425850340136) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 20.0 && pt <= 26.0) * (0.870406189555) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 26.0 && pt <= 32.0) * (0.948238022278) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 32.0 && pt <= 38.0) * (0.998082327993) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 38.0 && pt <= 44.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 44.0 && pt <= 50.0) * (0.917593936462) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 50.0 && pt <= 60.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 60.0 && pt <= 70.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 100.0 && pt <= 125.0) * (0.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 150.0) * (0.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 10.0 && pt <= 12.0) * (0.341587301587) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 12.0 && pt <= 14.0) * (0.163084702908) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 14.0 && pt <= 16.0) * (0.252991452991) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 16.0 && pt <= 18.0) * (0.628997867804) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 18.0 && pt <= 20.0) * (0.634399551066) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 20.0 && pt <= 26.0) * (0.747031882254) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 26.0 && pt <= 32.0) * (0.934474939401) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 32.0 && pt <= 38.0) * (0.976565008026) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 38.0 && pt <= 44.0) * (0.931815386944) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 44.0 && pt <= 50.0) * (0.885093577014) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 50.0 && pt <= 60.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 60.0 && pt <= 70.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 90.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 100.0 && pt <= 125.0) * (0.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 150.0) * (0.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 10.0 && pt <= 12.0) * (0.28209556707) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 12.0 && pt <= 14.0) * (0.371391076115) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 14.0 && pt <= 16.0) * (0.667046101309) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 16.0 && pt <= 18.0) * (0.300974025974) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 18.0 && pt <= 20.0) * (0.785468826706) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 20.0 && pt <= 26.0) * (0.676780446692) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 26.0 && pt <= 32.0) * (0.940191064364) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 32.0 && pt <= 38.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 38.0 && pt <= 44.0) * (0.987019322292) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 44.0 && pt <= 50.0) * (0.926719895002) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 50.0 && pt <= 60.0) * (0.878981645944) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 60.0 && pt <= 70.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 150.0) * (0.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 12.0) * (0.385885885886) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 12.0 && pt <= 14.0) * (0.69306122449) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 14.0 && pt <= 16.0) * (0.520964014526) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 16.0 && pt <= 18.0) * (0.76511954993) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 18.0 && pt <= 20.0) * (0.465714285714) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 26.0) * (0.879138099902) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 26.0 && pt <= 32.0) * (0.974956884162) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 32.0 && pt <= 38.0) * (0.978031821866) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 38.0 && pt <= 44.0) * (0.968240675127) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 44.0 && pt <= 50.0) * (0.990567420947) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 60.0) * (0.93300248139) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 60.0 && pt <= 70.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 80.0 && pt <= 90.0) * (0.666666666667) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 90.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 150.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 10.0 && pt <= 12.0) * (0.130116959064) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 12.0 && pt <= 14.0) * (0.189814814815) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 14.0 && pt <= 16.0) * (0.813981042654) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 16.0 && pt <= 18.0) * (0.711711711712) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 18.0 && pt <= 20.0) * (0.787974683544) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 20.0 && pt <= 26.0) * (0.724496426251) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 26.0 && pt <= 32.0) * (0.964941118744) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 32.0 && pt <= 38.0) * (0.840448097574) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 38.0 && pt <= 44.0) * (0.951287005502) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 44.0 && pt <= 50.0) * (0.943446088795) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 50.0 && pt <= 60.0) * (0.915346121768) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 60.0 && pt <= 70.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 80.0 && pt <= 90.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 150.0) * (1.0) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 10.0 && pt <= 12.0) * (0.574324324324) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 12.0 && pt <= 14.0) * (0.8125) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 14.0 && pt <= 16.0) * (0.221468926554) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 16.0 && pt <= 18.0) * (0.545226130653) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 18.0 && pt <= 20.0) * (0.525581395349) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 20.0 && pt <= 26.0) * (0.859853190288) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 26.0 && pt <= 32.0) * (0.845528455285) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 32.0 && pt <= 38.0) * (0.924585785234) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 38.0 && pt <= 44.0) * (0.924741132408) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 44.0 && pt <= 50.0) * (0.956875993641) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 50.0 && pt <= 60.0) * (1.0) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 60.0 && pt <= 70.0) * (0.805194805195) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 90.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 150.0) * (0.0) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 10.0 && pt <= 12.0) * (0.0) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 12.0 && pt <= 14.0) * (0.657627118644) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 14.0 && pt <= 16.0) * (0.252329749104) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 16.0 && pt <= 18.0) * (0.639) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 18.0 && pt <= 20.0) * (0.463958060288) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 20.0 && pt <= 26.0) * (0.891100702576) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 26.0 && pt <= 32.0) * (0.933317201327) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 32.0 && pt <= 38.0) * (0.947060933736) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 38.0 && pt <= 44.0) * (0.994242365923) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 44.0 && pt <= 50.0) * (0.950650293788) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 50.0 && pt <= 60.0) * (0.943315508021) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 60.0 && pt <= 70.0) * (0.671361502347) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 150.0) * (0.0) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 10.0 && pt <= 12.0) * (0.49756097561) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 12.0 && pt <= 14.0) * (0.363821138211) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 14.0 && pt <= 16.0) * (0.647872340426) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 16.0 && pt <= 18.0) * (0.594957983193) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 18.0 && pt <= 20.0) * (0.78021978022) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 20.0 && pt <= 26.0) * (0.648558988182) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 26.0 && pt <= 32.0) * (0.908861907401) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 32.0 && pt <= 38.0) * (0.94157782516) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 38.0 && pt <= 44.0) * (0.898918194725) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 44.0 && pt <= 50.0) * (0.914515079193) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 50.0 && pt <= 60.0) * (0.899598393574) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 60.0 && pt <= 70.0) * (0.806779661017) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 70.0 && pt <= 80.0) * (0.509803921569) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 150.0) * (0.0) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 10.0 && pt <= 12.0) * (0.758823529412) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 12.0 && pt <= 14.0) * (0.536666666667) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 14.0 && pt <= 16.0) * (0.608974358974) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 16.0 && pt <= 18.0) * (0.547619047619) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 18.0 && pt <= 20.0) * (0.856470588235) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 20.0 && pt <= 26.0) * (0.824600520253) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 26.0 && pt <= 32.0) * (0.816951896392) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 32.0 && pt <= 38.0) * (0.990975586587) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 38.0 && pt <= 44.0) * (1.0) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 44.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 50.0 && pt <= 60.0) * (0.864795918367) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 60.0 && pt <= 70.0) * (1.0) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 70.0 && pt <= 80.0) * (1.0) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 80.0 && pt <= 90.0) * (1.0) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 90.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 150.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 10.0 && pt <= 12.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 12.0 && pt <= 14.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 14.0 && pt <= 16.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 16.0 && pt <= 18.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 18.0 && pt <= 20.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 20.0 && pt <= 26.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 26.0 && pt <= 32.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 32.0 && pt <= 38.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 38.0 && pt <= 44.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 44.0 && pt <= 50.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 50.0 && pt <= 60.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 60.0 && pt <= 70.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 70.0 && pt <= 80.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 80.0 && pt <= 90.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 90.0 && pt <= 100.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 100.0 && pt <= 125.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 125.0 && pt <= 150.0) * (0.0) +
   (abs(eta) > 2.8) * (pt > 150.0) * (0.0)  
  }
}


module JetFlavorAssociation JetFlavorAssociationPUPPI {

  set PartonInputArray Delphes/partons
  set ParticleInputArray Delphes/allParticles
  set ParticleLHEFInputArray Delphes/allParticlesLHEF
  set JetInputArray JetSmearPUPPI/jets

  set DeltaR 0.5
  set PartonPTMin 10.0
  set PartonEtaMax 4.0

}

module JetFlavorAssociation JetFlavorAssociationPUPPIAK8 {

  set PartonInputArray Delphes/partons
  set ParticleInputArray Delphes/allParticles
  set ParticleLHEFInputArray Delphes/allParticlesLHEF
  set JetInputArray JetSmearPUPPIAK8/jets

  set DeltaR 0.8
  set PartonPTMin 100.0
  set PartonEtaMax 4.0

}




module BTagging BTaggingPUPPILoose {

  set JetInputArray JetSmearPUPPI/jets

  set BitNumber 0

  add EfficiencyFormula {0}      {0.1}

  add EfficiencyFormula {5}      {1.0}

  add EfficiencyFormula {4}      {1.0}
}

module BTagging BTaggingPUPPIMedium {

  set JetInputArray JetSmearPUPPI/jets

  set BitNumber 1

  add EfficiencyFormula {0}      {0.01}

  add EfficiencyFormula {5}      {1.0}

  add EfficiencyFormula {4}      {1.0}
}

module BTagging BTaggingPUPPITight {

  set JetInputArray JetSmearPUPPI/jets

  set BitNumber 2

  add EfficiencyFormula {0}      {0.001}

  add EfficiencyFormula {5}      {1.0}

  add EfficiencyFormula {4}      {1.0}

}


module BTagging BTaggingPUPPIAK8Loose {

  set JetInputArray JetSmearPUPPI/jets

  set BitNumber 0

  add EfficiencyFormula {0}      {0.1}

  add EfficiencyFormula {5}      {1.0}

  add EfficiencyFormula {4}      {1.0}
}

module BTagging BTaggingPUPPIAK8Medium {

  set JetInputArray JetSmearPUPPI/jets

  set BitNumber 1

  add EfficiencyFormula {0}      {0.01}

  add EfficiencyFormula {5}      {1.0}

  add EfficiencyFormula {4}      {1.0}
}

module BTagging BTaggingPUPPIAK8Tight {

  set JetInputArray JetSmearPUPPI/jets

  set BitNumber 2

  add EfficiencyFormula {0}      {0.001}

  add EfficiencyFormula {5}      {1.0}

  add EfficiencyFormula {4}      {1.0}

}




#############
# tau-tagging
#############


module TauTagging TauTaggingPUPPILoose {

  set ParticleInputArray Delphes/allParticles
  set PartonInputArray Delphes/partons
  set JetInputArray JetSmearPUPPI/jets



  set DeltaR 0.5

  set TauPTMin 20.0

  set TauEtaMax 2.3

  set BitNumber 0

  add EfficiencyFormula {0}      {0.1}

  add EfficiencyFormula {15}      {1.0}
}

module TauTagging TauTaggingPUPPIMedium {

  set ParticleInputArray Delphes/allParticles
  set PartonInputArray Delphes/partons
  set JetInputArray JetSmearPUPPI/jets

  set BitNumber 1

  add EfficiencyFormula {0}      {0.01}

  add EfficiencyFormula {15}      {1.0}
}

module TauTagging TauTaggingPUPPITight {

  set ParticleInputArray Delphes/allParticles
  set PartonInputArray Delphes/partons
  set JetInputArray JetSmearPUPPI/jets

  set BitNumber 2

  add EfficiencyFormula {0}      {0.001}

  add EfficiencyFormula {15}      {1.0}

}



#################################
# Jet Fake Particle Maker Loose #
#################################

module JetFakeParticle JetFakeMakerLoose {

  set InputArray JetSmearPUPPI/jets
  set PhotonOutputArray photons
  set MuonOutputArray muons
  set ElectronOutputArray electrons
  set JetOutputArray jets

  set EfficiencyFormula {
    {11} {

         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 20.0 && pt <= 26.0) * (0.0177838577291) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 26.0 && pt <= 32.0) * (0.0236998025016) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 32.0 && pt <= 38.0) * (0.0388109000826) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 38.0 && pt <= 44.0) * (0.042957042957) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 44.0 && pt <= 50.0) * (0.0395809080326) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 50.0 && pt <= 60.0) * (0.0496515679443) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 60.0 && pt <= 70.0) * (0.0508849557522) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 70.0 && pt <= 80.0) * (0.0523385300668) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 80.0 && pt <= 90.0) * (0.0414937759336) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 90.0 && pt <= 100.0) * (0.0540145985401) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 100.0 && pt <= 125.0) * (0.0343773873186) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 125.0 && pt <= 150.0) * (0.0242805755396) +
          (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 150.0) * (0.0121509064184) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 20.0 && pt <= 26.0) * (0.0121065375303) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 26.0 && pt <= 32.0) * (0.0254120879121) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 32.0 && pt <= 38.0) * (0.0290598290598) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 38.0 && pt <= 44.0) * (0.0438972162741) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 44.0 && pt <= 50.0) * (0.0504926108374) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 50.0 && pt <= 60.0) * (0.0577586206897) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 60.0 && pt <= 70.0) * (0.0580503833516) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 70.0 && pt <= 80.0) * (0.0425055928412) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 80.0 && pt <= 90.0) * (0.0564024390244) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 90.0 && pt <= 100.0) * (0.0497592295345) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 100.0 && pt <= 125.0) * (0.0367647058824) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 125.0 && pt <= 150.0) * (0.031746031746) +
          (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 150.0) * (0.01131057383) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 20.0 && pt <= 26.0) * (0.0132585415604) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 26.0 && pt <= 32.0) * (0.02284082798) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 32.0 && pt <= 38.0) * (0.0329289428076) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 38.0 && pt <= 44.0) * (0.0367567567568) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 44.0 && pt <= 50.0) * (0.0330073349633) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 50.0 && pt <= 60.0) * (0.0694716242661) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 60.0 && pt <= 70.0) * (0.058064516129) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 70.0 && pt <= 80.0) * (0.0434782608696) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 80.0 && pt <= 90.0) * (0.0544117647059) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 90.0 && pt <= 100.0) * (0.0466237942122) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 100.0 && pt <= 125.0) * (0.0435483870968) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 125.0 && pt <= 150.0) * (0.0303623898139) +
          (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 150.0) * (0.0121342708097) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 20.0 && pt <= 26.0) * (0.00850425212606) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 26.0 && pt <= 32.0) * (0.0206847360913) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 32.0 && pt <= 38.0) * (0.0206008583691) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 38.0 && pt <= 44.0) * (0.0478309232481) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 44.0 && pt <= 50.0) * (0.035761589404) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 50.0 && pt <= 60.0) * (0.0507662835249) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 60.0 && pt <= 70.0) * (0.0455089820359) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 70.0 && pt <= 80.0) * (0.0365168539326) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 80.0 && pt <= 90.0) * (0.0404530744337) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 90.0 && pt <= 100.0) * (0.0373665480427) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 100.0 && pt <= 125.0) * (0.03) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 125.0 && pt <= 150.0) * (0.02634467618) +
          (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 150.0) * (0.0113260371596) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 26.0) * (0.0198178896626) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 26.0 && pt <= 32.0) * (0.0359602142311) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 32.0 && pt <= 38.0) * (0.0522456461962) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 38.0 && pt <= 44.0) * (0.0610778443114) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 44.0 && pt <= 50.0) * (0.0743338008415) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 60.0) * (0.100303951368) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 60.0 && pt <= 70.0) * (0.110843373494) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 70.0 && pt <= 80.0) * (0.0970042796006) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 80.0 && pt <= 90.0) * (0.0984719864177) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 90.0 && pt <= 100.0) * (0.105263157895) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 125.0) * (0.0990740740741) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 125.0 && pt <= 150.0) * (0.101136363636) +
          (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 150.0) * (0.0815170756822) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 20.0 && pt <= 26.0) * (0.0278293135436) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 26.0 && pt <= 32.0) * (0.0567901234568) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 32.0 && pt <= 38.0) * (0.0692307692308) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 38.0 && pt <= 44.0) * (0.104364326376) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 44.0 && pt <= 50.0) * (0.123255813953) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 50.0 && pt <= 60.0) * (0.161129568106) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 60.0 && pt <= 70.0) * (0.136904761905) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 70.0 && pt <= 80.0) * (0.176334106729) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 80.0 && pt <= 90.0) * (0.132394366197) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 90.0 && pt <= 100.0) * (0.162079510703) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 100.0 && pt <= 125.0) * (0.14696485623) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 125.0 && pt <= 150.0) * (0.157142857143) +
          (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 150.0) * (0.121170553269) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 20.0 && pt <= 26.0) * (0.069696969697) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 26.0 && pt <= 32.0) * (0.106944444444) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 32.0 && pt <= 38.0) * (0.142335766423) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 38.0 && pt <= 44.0) * (0.160356347439) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 44.0 && pt <= 50.0) * (0.21760391198) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 50.0 && pt <= 60.0) * (0.166975881262) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 60.0 && pt <= 70.0) * (0.18021978022) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 70.0 && pt <= 80.0) * (0.191549295775) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 80.0 && pt <= 90.0) * (0.199404761905) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 90.0 && pt <= 100.0) * (0.185667752443) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 100.0 && pt <= 125.0) * (0.143344709898) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 125.0 && pt <= 150.0) * (0.148837209302) +
          (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 150.0) * (0.107552870091) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 20.0 && pt <= 26.0) * (0.0903155603917) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 26.0 && pt <= 32.0) * (0.138138138138) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 32.0 && pt <= 38.0) * (0.142095914742) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 38.0 && pt <= 44.0) * (0.152777777778) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 44.0 && pt <= 50.0) * (0.176470588235) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 50.0 && pt <= 60.0) * (0.190661478599) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 60.0 && pt <= 70.0) * (0.169230769231) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 70.0 && pt <= 80.0) * (0.224358974359) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 80.0 && pt <= 90.0) * (0.179401993355) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 90.0 && pt <= 100.0) * (0.227488151659) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 100.0 && pt <= 125.0) * (0.171656686627) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 125.0 && pt <= 150.0) * (0.130030959752) +
          (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 150.0) * (0.121716287215) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 20.0 && pt <= 26.0) * (0.118101545254) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 26.0 && pt <= 32.0) * (0.180762852405) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 32.0 && pt <= 38.0) * (0.200435729847) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 38.0 && pt <= 44.0) * (0.213068181818) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 44.0 && pt <= 50.0) * (0.201117318436) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 50.0 && pt <= 60.0) * (0.220454545455) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 60.0 && pt <= 70.0) * (0.264788732394) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 70.0 && pt <= 80.0) * (0.21768707483) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 80.0 && pt <= 90.0) * (0.232653061224) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 90.0 && pt <= 100.0) * (0.207070707071) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 100.0 && pt <= 125.0) * (0.244215938303) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 125.0 && pt <= 150.0) * (0.234693877551) +
          (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 150.0) * (0.181705809642) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 20.0 && pt <= 26.0) * (0.13258983891) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 26.0 && pt <= 32.0) * (0.175531914894) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 32.0 && pt <= 38.0) * (0.195852534562) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 38.0 && pt <= 44.0) * (0.251515151515) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 44.0 && pt <= 50.0) * (0.201365187713) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 50.0 && pt <= 60.0) * (0.245862884161) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 60.0 && pt <= 70.0) * (0.285714285714) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 70.0 && pt <= 80.0) * (0.235537190083) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 80.0 && pt <= 90.0) * (0.322115384615) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 90.0 && pt <= 100.0) * (0.25) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 100.0 && pt <= 125.0) * (0.227692307692) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 125.0 && pt <= 150.0) * (0.27135678392) +
          (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 150.0) * (0.249530956848) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 20.0 && pt <= 26.0) * (0.0810397553517) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 26.0 && pt <= 32.0) * (0.0983899821109) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 32.0 && pt <= 38.0) * (0.154430379747) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 38.0 && pt <= 44.0) * (0.197916666667) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 44.0 && pt <= 50.0) * (0.213740458015) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 50.0 && pt <= 60.0) * (0.182584269663) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 60.0 && pt <= 70.0) * (0.234848484848) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 70.0 && pt <= 80.0) * (0.221649484536) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 80.0 && pt <= 90.0) * (0.264705882353) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 90.0 && pt <= 100.0) * (0.256578947368) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 100.0 && pt <= 125.0) * (0.204724409449) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 125.0 && pt <= 150.0) * (0.253521126761) +
          (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 150.0) * (0.210691823899) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 26.0) * (0.0553691275168) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 26.0 && pt <= 32.0) * (0.101604278075) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 32.0 && pt <= 38.0) * (0.132258064516) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 38.0 && pt <= 44.0) * (0.151020408163) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 44.0 && pt <= 50.0) * (0.158823529412) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 60.0) * (0.178707224335) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 60.0 && pt <= 70.0) * (0.149758454106) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 70.0 && pt <= 80.0) * (0.25) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 80.0 && pt <= 90.0) * (0.158415841584) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 90.0 && pt <= 100.0) * (0.158415841584) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 125.0) * (0.2) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 125.0 && pt <= 150.0) * (0.18691588785) +
          (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 150.0) * (0.16) +
          (abs(eta) > 3.0) * (pt > 10.0 && pt <= 12.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 12.0 && pt <= 14.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 14.0 && pt <= 16.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 16.0 && pt <= 18.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 18.0 && pt <= 20.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 20.0 && pt <= 26.0) * (0.00197238658777) +
          (abs(eta) > 3.0) * (pt > 26.0 && pt <= 32.0) * (0.000980872976949) +
          (abs(eta) > 3.0) * (pt > 32.0 && pt <= 38.0) * (0.00310719834283) +
          (abs(eta) > 3.0) * (pt > 38.0 && pt <= 44.0) * (0.00147058823529) +
          (abs(eta) > 3.0) * (pt > 44.0 && pt <= 50.0) * (0.00692520775623) +
          (abs(eta) > 3.0) * (pt > 50.0 && pt <= 60.0) * (0.00589101620029) +
          (abs(eta) > 3.0) * (pt > 60.0 && pt <= 70.0) * (0.00579710144928) +
          (abs(eta) > 3.0) * (pt > 70.0 && pt <= 80.0) * (0.0162162162162) +
          (abs(eta) > 3.0) * (pt > 80.0 && pt <= 90.0) * (0.021897810219) +
          (abs(eta) > 3.0) * (pt > 90.0 && pt <= 100.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 100.0 && pt <= 125.0) * (0.015873015873) +
          (abs(eta) > 3.0) * (pt > 125.0 && pt <= 150.0) * (0.0294117647059) +
          (abs(eta) > 3.0) * (pt > 150.0) * (0.0149253731343)  
  }
    {13} {

         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 20.0 && pt <= 26.0) * (0.00136798905609) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 26.0 && pt <= 32.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 32.0 && pt <= 38.0) * (0.00165152766309) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 38.0 && pt <= 44.0) * (0.001998001998) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 44.0 && pt <= 50.0) * (0.0011641443539) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 50.0 && pt <= 60.0) * (0.000871080139373) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 60.0 && pt <= 70.0) * (0.00110619469027) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 70.0 && pt <= 80.0) * (0.0022271714922) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 80.0 && pt <= 90.0) * (0.00276625172891) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 90.0 && pt <= 100.0) * (0.0029197080292) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 100.0 && pt <= 125.0) * (0.00229182582124) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 150.0) * (0.000685938265556) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 20.0 && pt <= 26.0) * (0.000484261501211) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 26.0 && pt <= 32.0) * (0.000686813186813) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 32.0 && pt <= 38.0) * (0.000854700854701) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 38.0 && pt <= 44.0) * (0.00107066381156) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 44.0 && pt <= 50.0) * (0.0012315270936) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 50.0 && pt <= 60.0) * (0.000862068965517) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 60.0 && pt <= 70.0) * (0.00328587075575) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 70.0 && pt <= 80.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 80.0 && pt <= 90.0) * (0.0015243902439) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 90.0 && pt <= 100.0) * (0.0016051364366) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 100.0 && pt <= 125.0) * (0.00147058823529) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 125.0 && pt <= 150.0) * (0.000992063492063) +
          (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 150.0) * (0.00186780118294) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 20.0 && pt <= 26.0) * (0.00101988781234) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 26.0 && pt <= 32.0) * (0.00142755174875) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 32.0 && pt <= 38.0) * (0.000866551126516) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 38.0 && pt <= 44.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 50.0 && pt <= 60.0) * (0.000978473581213) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 60.0 && pt <= 70.0) * (0.0010752688172) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 70.0 && pt <= 80.0) * (0.00127877237852) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 80.0 && pt <= 90.0) * (0.00147058823529) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 100.0 && pt <= 125.0) * (0.00241935483871) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 150.0) * (0.00158766160127) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 20.0 && pt <= 26.0) * (0.00100050025013) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 26.0 && pt <= 32.0) * (0.000713266761769) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 32.0 && pt <= 38.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 38.0 && pt <= 44.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 50.0 && pt <= 60.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 60.0 && pt <= 70.0) * (0.00119760479042) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 70.0 && pt <= 80.0) * (0.00280898876404) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 80.0 && pt <= 90.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 100.0 && pt <= 125.0) * (0.000833333333333) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 150.0) * (0.0012725884449) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 26.0) * (0.000535618639529) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 26.0 && pt <= 32.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 32.0 && pt <= 38.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 38.0 && pt <= 44.0) * (0.00119760479042) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 60.0) * (0.00202634245187) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 60.0 && pt <= 70.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 70.0 && pt <= 80.0) * (0.00142653352354) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 80.0 && pt <= 90.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 125.0) * (0.000925925925926) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 150.0) * (0.00240260854642) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 20.0 && pt <= 26.0) * (0.00077101002313) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 26.0 && pt <= 32.0) * (0.00104931794334) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 32.0 && pt <= 38.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 38.0 && pt <= 44.0) * (0.00160256410256) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 50.0 && pt <= 60.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 60.0 && pt <= 70.0) * (0.00163132137031) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 70.0 && pt <= 80.0) * (0.00392156862745) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 80.0 && pt <= 90.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 125.0 && pt <= 150.0) * (0.00167504187605) +
          (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 150.0) * (0.00578703703704) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 20.0 && pt <= 26.0) * (0.000834028356964) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 26.0 && pt <= 32.0) * (0.0) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 32.0 && pt <= 38.0) * (0.0) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 38.0 && pt <= 44.0) * (0.0) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 50.0 && pt <= 60.0) * (0.0) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 60.0 && pt <= 70.0) * (0.00194174757282) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 70.0 && pt <= 80.0) * (0.0) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 80.0 && pt <= 90.0) * (0.0) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 90.0 && pt <= 100.0) * (0.00296735905045) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 100.0 && pt <= 125.0) * (0.00144927536232) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 125.0 && pt <= 150.0) * (0.00207039337474) +
          (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 150.0) * (0.00280583613917) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 20.0 && pt <= 26.0) * (0.0) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 26.0 && pt <= 32.0) * (0.00129198966408) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 32.0 && pt <= 38.0) * (0.0) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 38.0 && pt <= 44.0) * (0.0) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 50.0 && pt <= 60.0) * (0.00349040139616) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 60.0 && pt <= 70.0) * (0.0) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 70.0 && pt <= 80.0) * (0.0) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 80.0 && pt <= 90.0) * (0.003125) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 150.0) * (0.00417362270451) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 20.0 && pt <= 26.0) * (0.0) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 26.0 && pt <= 32.0) * (0.0) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 32.0 && pt <= 38.0) * (0.0) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 38.0 && pt <= 44.0) * (0.00505050505051) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 50.0 && pt <= 60.0) * (0.0) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 60.0 && pt <= 70.0) * (0.0) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 70.0 && pt <= 80.0) * (0.0031746031746) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 80.0 && pt <= 90.0) * (0.00375939849624) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 90.0 && pt <= 100.0) * (0.00434782608696) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 125.0 && pt <= 150.0) * (0.00369003690037) +
          (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 150.0) * (0.00431654676259) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 20.0 && pt <= 26.0) * (0.00123456790123) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 26.0 && pt <= 32.0) * (0.0) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 32.0 && pt <= 38.0) * (0.0) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 38.0 && pt <= 44.0) * (0.0) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 50.0 && pt <= 60.0) * (0.00232558139535) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 60.0 && pt <= 70.0) * (0.0) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 70.0 && pt <= 80.0) * (0.0) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 80.0 && pt <= 90.0) * (0.00480769230769) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 150.0) * (0.00254452926209) +
          (abs(eta) > 2.8) * (pt > 10.0 && pt <= 12.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 12.0 && pt <= 14.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 14.0 && pt <= 16.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 16.0 && pt <= 18.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 18.0 && pt <= 20.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 20.0 && pt <= 26.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 26.0 && pt <= 32.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 32.0 && pt <= 38.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 38.0 && pt <= 44.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 44.0 && pt <= 50.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 50.0 && pt <= 60.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 60.0 && pt <= 70.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 70.0 && pt <= 80.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 80.0 && pt <= 90.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 90.0 && pt <= 100.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 100.0 && pt <= 125.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 150.0) * (0.0)  
  }
    {22} {

         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 20.0 && pt <= 26.0) * (0.0077519379845) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 26.0 && pt <= 32.0) * (0.00789993416722) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 32.0 && pt <= 38.0) * (0.0198183319571) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 38.0 && pt <= 44.0) * (0.016983016983) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 44.0 && pt <= 50.0) * (0.0151338766007) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 50.0 && pt <= 60.0) * (0.00783972125436) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 60.0 && pt <= 70.0) * (0.0121681415929) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 70.0 && pt <= 80.0) * (0.0100222717149) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 80.0 && pt <= 90.0) * (0.00553250345781) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 90.0 && pt <= 100.0) * (0.00583941605839) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 100.0 && pt <= 125.0) * (0.00381970970206) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 125.0 && pt <= 150.0) * (0.0089928057554) +
          (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 150.0) * (0.00342969132778) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 20.0 && pt <= 26.0) * (0.00871670702179) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 26.0 && pt <= 32.0) * (0.0130494505495) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 32.0 && pt <= 38.0) * (0.0136752136752) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 38.0 && pt <= 44.0) * (0.00963597430407) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 44.0 && pt <= 50.0) * (0.0160098522167) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 50.0 && pt <= 60.0) * (0.0146551724138) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 60.0 && pt <= 70.0) * (0.0109529025192) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 70.0 && pt <= 80.0) * (0.0111856823266) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 80.0 && pt <= 90.0) * (0.00762195121951) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 90.0 && pt <= 100.0) * (0.00321027287319) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 100.0 && pt <= 125.0) * (0.00514705882353) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 125.0 && pt <= 150.0) * (0.00694444444444) +
          (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 150.0) * (0.00332053543634) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 20.0 && pt <= 26.0) * (0.0117287098419) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 26.0 && pt <= 32.0) * (0.0164168451106) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 32.0 && pt <= 38.0) * (0.0190641247834) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 38.0 && pt <= 44.0) * (0.0183783783784) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 44.0 && pt <= 50.0) * (0.0146699266504) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 50.0 && pt <= 60.0) * (0.0225048923679) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 60.0 && pt <= 70.0) * (0.0182795698925) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 70.0 && pt <= 80.0) * (0.0115089514066) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 80.0 && pt <= 90.0) * (0.00882352941176) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 90.0 && pt <= 100.0) * (0.0112540192926) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 100.0 && pt <= 125.0) * (0.00483870967742) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 125.0 && pt <= 150.0) * (0.00195886385896) +
          (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 150.0) * (0.00476298480381) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 20.0 && pt <= 26.0) * (0.0180090045023) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 26.0 && pt <= 32.0) * (0.0235378031384) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 32.0 && pt <= 38.0) * (0.0334763948498) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 38.0 && pt <= 44.0) * (0.0266963292547) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 44.0 && pt <= 50.0) * (0.0344370860927) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 50.0 && pt <= 60.0) * (0.0325670498084) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 60.0 && pt <= 70.0) * (0.0311377245509) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 70.0 && pt <= 80.0) * (0.0266853932584) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 80.0 && pt <= 90.0) * (0.0194174757282) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 90.0 && pt <= 100.0) * (0.00355871886121) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 100.0 && pt <= 125.0) * (0.0125) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 125.0 && pt <= 150.0) * (0.00768386388584) +
          (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 150.0) * (0.00547213031306) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 26.0) * (0.0385645420461) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 26.0 && pt <= 32.0) * (0.0413159908187) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 32.0 && pt <= 38.0) * (0.0650779101742) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 38.0 && pt <= 44.0) * (0.0502994011976) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 44.0 && pt <= 50.0) * (0.070126227209) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 60.0) * (0.0658561296859) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 60.0 && pt <= 70.0) * (0.0674698795181) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 70.0 && pt <= 80.0) * (0.0599144079886) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 80.0 && pt <= 90.0) * (0.0424448217317) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 90.0 && pt <= 100.0) * (0.0607287449393) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 125.0) * (0.0314814814815) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 125.0 && pt <= 150.0) * (0.0227272727273) +
          (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 150.0) * (0.0202505577484) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 20.0 && pt <= 26.0) * (0.0227882037534) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 26.0 && pt <= 32.0) * (0.0363636363636) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 32.0 && pt <= 38.0) * (0.0249433106576) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 38.0 && pt <= 44.0) * (0.0253164556962) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 44.0 && pt <= 50.0) * (0.0272572402044) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60.0) * (0.0243013365735) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70.0) * (0.0202898550725) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80.0) * (0.0297202797203) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90.0) * (0.0100200400802) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100.0) * (0.0133333333333) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 100.0 && pt <= 125.0) * (0.0102857142857) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 125.0 && pt <= 150.0) * (0.0103397341211) +
          (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 150.0) * (0.00239316239316) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 20.0 && pt <= 26.0) * (0.0112023898432) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 26.0 && pt <= 32.0) * (0.0103305785124) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 32.0 && pt <= 38.0) * (0.00259067357513) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 38.0 && pt <= 44.0) * (0.00480769230769) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 44.0 && pt <= 50.0) * (0.00376647834275) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 50.0 && pt <= 60.0) * (0.00414937759336) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 60.0 && pt <= 70.0) * (0.00171526586621) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 70.0 && pt <= 80.0) * (0.0) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 80.0 && pt <= 90.0) * (0.0) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 150.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 26.0) * (0.00563153660499) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 26.0 && pt <= 32.0) * (0.00346820809249) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 32.0 && pt <= 38.0) * (0.00151975683891) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 38.0 && pt <= 44.0) * (0.0020618556701) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 150.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 20.0 && pt <= 26.0) * (0.000973709834469) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 26.0 && pt <= 32.0) * (0.0012987012987) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 32.0 && pt <= 38.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 38.0 && pt <= 44.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 50.0 && pt <= 60.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 60.0 && pt <= 70.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 70.0 && pt <= 80.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 80.0 && pt <= 90.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 150.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 26.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 26.0 && pt <= 32.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 32.0 && pt <= 38.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 38.0 && pt <= 44.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 60.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 60.0 && pt <= 70.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 70.0 && pt <= 80.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 80.0 && pt <= 90.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 150.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 10.0 && pt <= 12.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 12.0 && pt <= 14.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 14.0 && pt <= 16.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 16.0 && pt <= 18.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 18.0 && pt <= 20.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 20.0 && pt <= 26.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 26.0 && pt <= 32.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 32.0 && pt <= 38.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 38.0 && pt <= 44.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 44.0 && pt <= 50.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 50.0 && pt <= 60.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 60.0 && pt <= 70.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 70.0 && pt <= 80.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 80.0 && pt <= 90.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 90.0 && pt <= 100.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 100.0 && pt <= 125.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 150.0) * (0.0)  
  }
  }

}


#################################
# Jet Fake Particle Maker Medium #
#################################

module JetFakeParticle JetFakeMakerMedium {

  set InputArray JetSmearPUPPI/jets
  set PhotonOutputArray photons
  set MuonOutputArray muons
  set ElectronOutputArray electrons
  set JetOutputArray jets

  set EfficiencyFormula {
    {11} {

         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 20.0 && pt <= 26.0) * (0.00273597811218) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 26.0 && pt <= 32.0) * (0.00658327847268) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 32.0 && pt <= 38.0) * (0.00495458298927) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 38.0 && pt <= 44.0) * (0.004995004995) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 44.0 && pt <= 50.0) * (0.0023282887078) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 50.0 && pt <= 60.0) * (0.00522648083624) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 60.0 && pt <= 70.0) * (0.00553097345133) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 70.0 && pt <= 80.0) * (0.00668151447661) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 80.0 && pt <= 90.0) * (0.00276625172891) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 90.0 && pt <= 100.0) * (0.00729927007299) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 100.0 && pt <= 125.0) * (0.00381970970206) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 125.0 && pt <= 150.0) * (0.00089928057554) +
          (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 150.0) * (0.00127388535032) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 20.0 && pt <= 26.0) * (0.000968523002421) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 26.0 && pt <= 32.0) * (0.00343406593407) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 32.0 && pt <= 38.0) * (0.0034188034188) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 38.0 && pt <= 44.0) * (0.00428265524625) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 44.0 && pt <= 50.0) * (0.00246305418719) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 50.0 && pt <= 60.0) * (0.00689655172414) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 60.0 && pt <= 70.0) * (0.00766703176342) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 70.0 && pt <= 80.0) * (0.00335570469799) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 80.0 && pt <= 90.0) * (0.00762195121951) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 90.0 && pt <= 100.0) * (0.00802568218299) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 100.0 && pt <= 125.0) * (0.00220588235294) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 125.0 && pt <= 150.0) * (0.00297619047619) +
          (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 150.0) * (0.000518833661928) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 20.0 && pt <= 26.0) * (0.00203977562468) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 26.0 && pt <= 32.0) * (0.00428265524625) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 32.0 && pt <= 38.0) * (0.00173310225303) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 38.0 && pt <= 44.0) * (0.00216216216216) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 44.0 && pt <= 50.0) * (0.00122249388753) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 50.0 && pt <= 60.0) * (0.0117416829746) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 60.0 && pt <= 70.0) * (0.00645161290323) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 70.0 && pt <= 80.0) * (0.00383631713555) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 80.0 && pt <= 90.0) * (0.00147058823529) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 90.0 && pt <= 100.0) * (0.0048231511254) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 100.0 && pt <= 125.0) * (0.00241935483871) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 125.0 && pt <= 150.0) * (0.00195886385896) +
          (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 150.0) * (0.00113404400091) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 20.0 && pt <= 26.0) * (0.00100050025013) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 26.0 && pt <= 32.0) * (0.00285306704708) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 32.0 && pt <= 38.0) * (0.00343347639485) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 38.0 && pt <= 44.0) * (0.00222469410456) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 44.0 && pt <= 50.0) * (0.00264900662252) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 50.0 && pt <= 60.0) * (0.00383141762452) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 60.0 && pt <= 70.0) * (0.00359281437126) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 70.0 && pt <= 80.0) * (0.00140449438202) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 80.0 && pt <= 90.0) * (0.00323624595469) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 90.0 && pt <= 100.0) * (0.00355871886121) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 100.0 && pt <= 125.0) * (0.00166666666667) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 125.0 && pt <= 150.0) * (0.00439077936334) +
          (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 150.0) * (0.000381776533469) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 26.0) * (0.00267809319764) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 26.0 && pt <= 32.0) * (0.00459066564652) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 32.0 && pt <= 38.0) * (0.00733272227314) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 38.0 && pt <= 44.0) * (0.00838323353293) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 44.0 && pt <= 50.0) * (0.00841514726508) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 60.0) * (0.0081053698075) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 60.0 && pt <= 70.0) * (0.0120481927711) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 70.0 && pt <= 80.0) * (0.00998573466476) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 80.0 && pt <= 90.0) * (0.00339558573854) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 90.0 && pt <= 100.0) * (0.00202429149798) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 125.0) * (0.00555555555556) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 125.0 && pt <= 150.0) * (0.00909090909091) +
          (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 150.0) * (0.0082375150163) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 20.0 && pt <= 26.0) * (0.00185528756957) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 26.0 && pt <= 32.0) * (0.00740740740741) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 32.0 && pt <= 38.0) * (0.00615384615385) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 38.0 && pt <= 44.0) * (0.00759013282732) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 44.0 && pt <= 50.0) * (0.00232558139535) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 50.0 && pt <= 60.0) * (0.0116279069767) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 60.0 && pt <= 70.0) * (0.00396825396825) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 70.0 && pt <= 80.0) * (0.0162412993039) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 80.0 && pt <= 90.0) * (0.0056338028169) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 90.0 && pt <= 100.0) * (0.0183486238532) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 100.0 && pt <= 125.0) * (0.0143769968051) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 125.0 && pt <= 150.0) * (0.0142857142857) +
          (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 150.0) * (0.00960219478738) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 20.0 && pt <= 26.0) * (0.00808080808081) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 26.0 && pt <= 32.0) * (0.0180555555556) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 32.0 && pt <= 38.0) * (0.0127737226277) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 38.0 && pt <= 44.0) * (0.0133630289532) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 44.0 && pt <= 50.0) * (0.0171149144254) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 50.0 && pt <= 60.0) * (0.0204081632653) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 60.0 && pt <= 70.0) * (0.00879120879121) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 70.0 && pt <= 80.0) * (0.0197183098592) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 80.0 && pt <= 90.0) * (0.0178571428571) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 90.0 && pt <= 100.0) * (0.00651465798046) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 100.0 && pt <= 125.0) * (0.0204778156997) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 125.0 && pt <= 150.0) * (0.0046511627907) +
          (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 150.0) * (0.00906344410876) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 20.0 && pt <= 26.0) * (0.00652883569097) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 26.0 && pt <= 32.0) * (0.015015015015) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 32.0 && pt <= 38.0) * (0.0106571936057) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 38.0 && pt <= 44.0) * (0.0185185185185) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 44.0 && pt <= 50.0) * (0.0147058823529) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 50.0 && pt <= 60.0) * (0.0136186770428) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 60.0 && pt <= 70.0) * (0.0128205128205) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 70.0 && pt <= 80.0) * (0.0192307692308) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 80.0 && pt <= 90.0) * (0.00664451827243) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 90.0 && pt <= 100.0) * (0.0236966824645) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 100.0 && pt <= 125.0) * (0.0059880239521) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 125.0 && pt <= 150.0) * (0.0030959752322) +
          (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 150.0) * (0.00175131348511) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 20.0 && pt <= 26.0) * (0.0165562913907) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 26.0 && pt <= 32.0) * (0.0364842454395) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 32.0 && pt <= 38.0) * (0.0348583877996) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 38.0 && pt <= 44.0) * (0.0198863636364) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 44.0 && pt <= 50.0) * (0.0139664804469) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 50.0 && pt <= 60.0) * (0.0318181818182) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 60.0 && pt <= 70.0) * (0.0169014084507) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 70.0 && pt <= 80.0) * (0.0238095238095) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 80.0 && pt <= 90.0) * (0.0244897959184) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 90.0 && pt <= 100.0) * (0.020202020202) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 100.0 && pt <= 125.0) * (0.0411311053985) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 125.0 && pt <= 150.0) * (0.0136054421769) +
          (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 150.0) * (0.0160692212608) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 20.0 && pt <= 26.0) * (0.0198265179678) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 26.0 && pt <= 32.0) * (0.0230496453901) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 32.0 && pt <= 38.0) * (0.0345622119816) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 38.0 && pt <= 44.0) * (0.0363636363636) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 44.0 && pt <= 50.0) * (0.0477815699659) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 50.0 && pt <= 60.0) * (0.0307328605201) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 60.0 && pt <= 70.0) * (0.0357142857143) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 70.0 && pt <= 80.0) * (0.0330578512397) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 80.0 && pt <= 90.0) * (0.0144230769231) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 90.0 && pt <= 100.0) * (0.015625) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 100.0 && pt <= 125.0) * (0.0307692307692) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 125.0 && pt <= 150.0) * (0.0150753768844) +
          (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 150.0) * (0.0168855534709) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 20.0 && pt <= 26.0) * (0.0183486238532) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 26.0 && pt <= 32.0) * (0.0161001788909) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 32.0 && pt <= 38.0) * (0.0278481012658) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 38.0 && pt <= 44.0) * (0.03125) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 44.0 && pt <= 50.0) * (0.030534351145) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 50.0 && pt <= 60.0) * (0.00842696629213) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 60.0 && pt <= 70.0) * (0.0378787878788) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 70.0 && pt <= 80.0) * (0.0103092783505) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 80.0 && pt <= 90.0) * (0.0352941176471) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 90.0 && pt <= 100.0) * (0.0328947368421) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 100.0 && pt <= 125.0) * (0.0196850393701) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 125.0 && pt <= 150.0) * (0.00704225352113) +
          (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 150.0) * (0.0345911949686) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 26.0) * (0.00838926174497) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 26.0 && pt <= 32.0) * (0.00267379679144) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 32.0 && pt <= 38.0) * (0.0129032258065) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 38.0 && pt <= 44.0) * (0.0122448979592) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 44.0 && pt <= 50.0) * (0.0352941176471) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 60.0) * (0.0228136882129) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 60.0 && pt <= 70.0) * (0.0193236714976) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 70.0 && pt <= 80.0) * (0.0263157894737) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 80.0 && pt <= 90.0) * (0.029702970297) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 90.0 && pt <= 100.0) * (0.029702970297) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 125.0) * (0.0451612903226) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 125.0 && pt <= 150.0) * (0.00934579439252) +
          (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 150.0) * (0.0457142857143) +
          (abs(eta) > 3.0) * (pt > 10.0 && pt <= 12.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 12.0 && pt <= 14.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 14.0 && pt <= 16.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 16.0 && pt <= 18.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 18.0 && pt <= 20.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 20.0 && pt <= 26.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 26.0 && pt <= 32.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 32.0 && pt <= 38.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 38.0 && pt <= 44.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 44.0 && pt <= 50.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 50.0 && pt <= 60.0) * (0.00294550810015) +
          (abs(eta) > 3.0) * (pt > 60.0 && pt <= 70.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 70.0 && pt <= 80.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 80.0 && pt <= 90.0) * (0.00729927007299) +
          (abs(eta) > 3.0) * (pt > 90.0 && pt <= 100.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 100.0 && pt <= 125.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 150.0) * (0.0)  
  }
    {13} {

         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 20.0 && pt <= 26.0) * (0.00136798905609) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 26.0 && pt <= 32.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 32.0 && pt <= 38.0) * (0.000825763831544) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 38.0 && pt <= 44.0) * (0.001998001998) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 44.0 && pt <= 50.0) * (0.0011641443539) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 50.0 && pt <= 60.0) * (0.000871080139373) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 60.0 && pt <= 70.0) * (0.00110619469027) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 70.0 && pt <= 80.0) * (0.0022271714922) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 80.0 && pt <= 90.0) * (0.00276625172891) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 90.0 && pt <= 100.0) * (0.0029197080292) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 100.0 && pt <= 125.0) * (0.00229182582124) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 150.0) * (0.000587947084762) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 20.0 && pt <= 26.0) * (0.000484261501211) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 26.0 && pt <= 32.0) * (0.000686813186813) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 32.0 && pt <= 38.0) * (0.000854700854701) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 38.0 && pt <= 44.0) * (0.00107066381156) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 50.0 && pt <= 60.0) * (0.000862068965517) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 60.0 && pt <= 70.0) * (0.00328587075575) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 70.0 && pt <= 80.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 80.0 && pt <= 90.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 90.0 && pt <= 100.0) * (0.0016051364366) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 100.0 && pt <= 125.0) * (0.00147058823529) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 125.0 && pt <= 150.0) * (0.000992063492063) +
          (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 150.0) * (0.00166026771817) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 20.0 && pt <= 26.0) * (0.00101988781234) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 26.0 && pt <= 32.0) * (0.00142755174875) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 32.0 && pt <= 38.0) * (0.000866551126516) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 38.0 && pt <= 44.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 50.0 && pt <= 60.0) * (0.000978473581213) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 60.0 && pt <= 70.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 70.0 && pt <= 80.0) * (0.00127877237852) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 80.0 && pt <= 90.0) * (0.00147058823529) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 100.0 && pt <= 125.0) * (0.00161290322581) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 150.0) * (0.00158766160127) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 20.0 && pt <= 26.0) * (0.00100050025013) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 26.0 && pt <= 32.0) * (0.000713266761769) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 32.0 && pt <= 38.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 38.0 && pt <= 44.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 50.0 && pt <= 60.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 60.0 && pt <= 70.0) * (0.00119760479042) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 70.0 && pt <= 80.0) * (0.00280898876404) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 80.0 && pt <= 90.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 100.0 && pt <= 125.0) * (0.000833333333333) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 150.0) * (0.0012725884449) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 26.0) * (0.000535618639529) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 26.0 && pt <= 32.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 32.0 && pt <= 38.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 38.0 && pt <= 44.0) * (0.00119760479042) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 60.0) * (0.00202634245187) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 60.0 && pt <= 70.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 70.0 && pt <= 80.0) * (0.00142653352354) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 80.0 && pt <= 90.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 125.0) * (0.000925925925926) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 150.0) * (0.00205937875408) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 20.0 && pt <= 26.0) * (0.00077101002313) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 26.0 && pt <= 32.0) * (0.00104931794334) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 32.0 && pt <= 38.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 38.0 && pt <= 44.0) * (0.00160256410256) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 50.0 && pt <= 60.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 60.0 && pt <= 70.0) * (0.00163132137031) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 70.0 && pt <= 80.0) * (0.00392156862745) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 80.0 && pt <= 90.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 125.0 && pt <= 150.0) * (0.00167504187605) +
          (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 150.0) * (0.00424382716049) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 20.0 && pt <= 26.0) * (0.0) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 26.0 && pt <= 32.0) * (0.0) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 32.0 && pt <= 38.0) * (0.0) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 38.0 && pt <= 44.0) * (0.0) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 50.0 && pt <= 60.0) * (0.0) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 60.0 && pt <= 70.0) * (0.00194174757282) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 70.0 && pt <= 80.0) * (0.0) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 80.0 && pt <= 90.0) * (0.0) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 100.0 && pt <= 125.0) * (0.00144927536232) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 125.0 && pt <= 150.0) * (0.00207039337474) +
          (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 150.0) * (0.00280583613917) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 20.0 && pt <= 26.0) * (0.0) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 26.0 && pt <= 32.0) * (0.00129198966408) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 32.0 && pt <= 38.0) * (0.0) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 38.0 && pt <= 44.0) * (0.0) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 50.0 && pt <= 60.0) * (0.00174520069808) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 60.0 && pt <= 70.0) * (0.0) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 70.0 && pt <= 80.0) * (0.0) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 80.0 && pt <= 90.0) * (0.003125) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 150.0) * (0.0025041736227) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 20.0 && pt <= 26.0) * (0.0) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 26.0 && pt <= 32.0) * (0.0) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 32.0 && pt <= 38.0) * (0.0) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 38.0 && pt <= 44.0) * (0.00505050505051) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 50.0 && pt <= 60.0) * (0.0) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 60.0 && pt <= 70.0) * (0.0) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 70.0 && pt <= 80.0) * (0.0) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 80.0 && pt <= 90.0) * (0.00375939849624) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 90.0 && pt <= 100.0) * (0.00434782608696) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 125.0 && pt <= 150.0) * (0.00369003690037) +
          (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 150.0) * (0.00287769784173) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 20.0 && pt <= 26.0) * (0.00123456790123) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 26.0 && pt <= 32.0) * (0.0) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 32.0 && pt <= 38.0) * (0.0) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 38.0 && pt <= 44.0) * (0.0) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 50.0 && pt <= 60.0) * (0.00232558139535) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 60.0 && pt <= 70.0) * (0.0) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 70.0 && pt <= 80.0) * (0.0) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 80.0 && pt <= 90.0) * (0.00480769230769) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 150.0) * (0.00254452926209) +
          (abs(eta) > 2.8) * (pt > 10.0 && pt <= 12.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 12.0 && pt <= 14.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 14.0 && pt <= 16.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 16.0 && pt <= 18.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 18.0 && pt <= 20.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 20.0 && pt <= 26.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 26.0 && pt <= 32.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 32.0 && pt <= 38.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 38.0 && pt <= 44.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 44.0 && pt <= 50.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 50.0 && pt <= 60.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 60.0 && pt <= 70.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 70.0 && pt <= 80.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 80.0 && pt <= 90.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 90.0 && pt <= 100.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 100.0 && pt <= 125.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 150.0) * (0.0)  
  }
    {22} {

         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 20.0 && pt <= 26.0) * (0.00638394892841) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 26.0 && pt <= 32.0) * (0.00592495062541) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 32.0 && pt <= 38.0) * (0.0156895127993) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 38.0 && pt <= 44.0) * (0.013986013986) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 44.0 && pt <= 50.0) * (0.0104772991851) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 50.0 && pt <= 60.0) * (0.00435540069686) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 60.0 && pt <= 70.0) * (0.00774336283186) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 70.0 && pt <= 80.0) * (0.00445434298441) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 80.0 && pt <= 90.0) * (0.00276625172891) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 90.0 && pt <= 100.0) * (0.0029197080292) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 100.0 && pt <= 125.0) * (0.000763941940413) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 125.0 && pt <= 150.0) * (0.00719424460432) +
          (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 150.0) * (0.00264576188143) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 20.0 && pt <= 26.0) * (0.00774818401937) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 26.0 && pt <= 32.0) * (0.0116758241758) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 32.0 && pt <= 38.0) * (0.0136752136752) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 38.0 && pt <= 44.0) * (0.00856531049251) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 44.0 && pt <= 50.0) * (0.012315270936) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 50.0 && pt <= 60.0) * (0.00775862068966) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 60.0 && pt <= 70.0) * (0.00876232201533) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 70.0 && pt <= 80.0) * (0.00671140939597) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 80.0 && pt <= 90.0) * (0.0030487804878) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 90.0 && pt <= 100.0) * (0.0016051364366) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 100.0 && pt <= 125.0) * (0.00367647058824) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 125.0 && pt <= 150.0) * (0.00595238095238) +
          (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 150.0) * (0.0029054685068) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 20.0 && pt <= 26.0) * (0.00968893421724) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 26.0 && pt <= 32.0) * (0.0149892933619) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 32.0 && pt <= 38.0) * (0.0181975736568) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 38.0 && pt <= 44.0) * (0.0108108108108) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 44.0 && pt <= 50.0) * (0.0110024449878) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 50.0 && pt <= 60.0) * (0.0166340508806) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 60.0 && pt <= 70.0) * (0.0139784946237) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 70.0 && pt <= 80.0) * (0.00511508951407) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 80.0 && pt <= 90.0) * (0.00147058823529) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 90.0 && pt <= 100.0) * (0.008038585209) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 100.0 && pt <= 125.0) * (0.00322580645161) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 125.0 && pt <= 150.0) * (0.00195886385896) +
          (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 150.0) * (0.00408255840327) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 20.0 && pt <= 26.0) * (0.016008004002) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 26.0 && pt <= 32.0) * (0.0192582025678) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 32.0 && pt <= 38.0) * (0.0223175965665) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 38.0 && pt <= 44.0) * (0.0233592880979) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 44.0 && pt <= 50.0) * (0.0238410596026) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 50.0 && pt <= 60.0) * (0.0239463601533) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 60.0 && pt <= 70.0) * (0.0203592814371) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 70.0 && pt <= 80.0) * (0.0196629213483) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 80.0 && pt <= 90.0) * (0.00809061488673) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 90.0 && pt <= 100.0) * (0.0017793594306) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 100.0 && pt <= 125.0) * (0.00916666666667) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 125.0 && pt <= 150.0) * (0.00439077936334) +
          (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 150.0) * (0.0036905064902) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 26.0) * (0.0353508302089) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 26.0 && pt <= 32.0) * (0.0359602142311) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 32.0 && pt <= 38.0) * (0.0540788267644) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 38.0 && pt <= 44.0) * (0.037125748503) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 44.0 && pt <= 50.0) * (0.0504908835905) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 60.0) * (0.048632218845) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 60.0 && pt <= 70.0) * (0.0542168674699) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 70.0 && pt <= 80.0) * (0.0413694721826) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 80.0 && pt <= 90.0) * (0.025466893039) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 90.0 && pt <= 100.0) * (0.0485829959514) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 125.0) * (0.0212962962963) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 125.0 && pt <= 150.0) * (0.0204545454545) +
          (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 150.0) * (0.0152737257594) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 20.0 && pt <= 26.0) * (0.0154155495979) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 26.0 && pt <= 32.0) * (0.0218181818182) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 32.0 && pt <= 38.0) * (0.0136054421769) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 38.0 && pt <= 44.0) * (0.0168776371308) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 44.0 && pt <= 50.0) * (0.015332197615) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60.0) * (0.0145808019441) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70.0) * (0.0144927536232) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80.0) * (0.020979020979) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90.0) * (0.00801603206413) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100.0) * (0.00666666666667) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 100.0 && pt <= 125.0) * (0.00457142857143) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 125.0 && pt <= 150.0) * (0.0103397341211) +
          (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 150.0) * (0.0017094017094) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 20.0 && pt <= 26.0) * (0.00597460791636) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 26.0 && pt <= 32.0) * (0.0051652892562) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 32.0 && pt <= 38.0) * (0.00259067357513) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 38.0 && pt <= 44.0) * (0.00480769230769) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 44.0 && pt <= 50.0) * (0.00188323917137) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 50.0 && pt <= 60.0) * (0.00138312586445) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 60.0 && pt <= 70.0) * (0.0) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 70.0 && pt <= 80.0) * (0.0) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 80.0 && pt <= 90.0) * (0.0) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 150.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 26.0) * (0.00160901045857) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 26.0 && pt <= 32.0) * (0.00346820809249) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 32.0 && pt <= 38.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 38.0 && pt <= 44.0) * (0.0020618556701) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 150.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 20.0 && pt <= 26.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 26.0 && pt <= 32.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 32.0 && pt <= 38.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 38.0 && pt <= 44.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 50.0 && pt <= 60.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 60.0 && pt <= 70.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 70.0 && pt <= 80.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 80.0 && pt <= 90.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 150.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 26.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 26.0 && pt <= 32.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 32.0 && pt <= 38.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 38.0 && pt <= 44.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 60.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 60.0 && pt <= 70.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 70.0 && pt <= 80.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 80.0 && pt <= 90.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 150.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 10.0 && pt <= 12.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 12.0 && pt <= 14.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 14.0 && pt <= 16.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 16.0 && pt <= 18.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 18.0 && pt <= 20.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 20.0 && pt <= 26.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 26.0 && pt <= 32.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 32.0 && pt <= 38.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 38.0 && pt <= 44.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 44.0 && pt <= 50.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 50.0 && pt <= 60.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 60.0 && pt <= 70.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 70.0 && pt <= 80.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 80.0 && pt <= 90.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 90.0 && pt <= 100.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 100.0 && pt <= 125.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 150.0) * (0.0)  
  }
  }

}

#################################
# Jet Fake Particle Maker Tight #
#################################

module JetFakeParticle JetFakeMakerTight {

  set InputArray JetSmearPUPPI/jets
  set PhotonOutputArray photons
  set MuonOutputArray muons
  set ElectronOutputArray electrons
  set JetOutputArray jets

  set EfficiencyFormula {
    {11} {

         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 20.0 && pt <= 26.0) * (0.000455996352029) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 26.0 && pt <= 32.0) * (0.0019749835418) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 32.0 && pt <= 38.0) * (0.000825763831544) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 38.0 && pt <= 44.0) * (0.000999000999001) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 50.0 && pt <= 60.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 60.0 && pt <= 70.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 70.0 && pt <= 80.0) * (0.0022271714922) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 80.0 && pt <= 90.0) * (0.00138312586445) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 90.0 && pt <= 100.0) * (0.0029197080292) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 150.0) * (0.000391964723175) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 20.0 && pt <= 26.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 26.0 && pt <= 32.0) * (0.000686813186813) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 32.0 && pt <= 38.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 38.0 && pt <= 44.0) * (0.00107066381156) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 50.0 && pt <= 60.0) * (0.000862068965517) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 60.0 && pt <= 70.0) * (0.00109529025192) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 70.0 && pt <= 80.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 80.0 && pt <= 90.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 90.0 && pt <= 100.0) * (0.0016051364366) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 150.0) * (0.000103766732386) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 20.0 && pt <= 26.0) * (0.00050994390617) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 26.0 && pt <= 32.0) * (0.000713775874375) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 32.0 && pt <= 38.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 38.0 && pt <= 44.0) * (0.00216216216216) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 50.0 && pt <= 60.0) * (0.000978473581213) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 60.0 && pt <= 70.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 70.0 && pt <= 80.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 80.0 && pt <= 90.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 90.0 && pt <= 100.0) * (0.0016077170418) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 150.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 20.0 && pt <= 26.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 26.0 && pt <= 32.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 32.0 && pt <= 38.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 38.0 && pt <= 44.0) * (0.00111234705228) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 50.0 && pt <= 60.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 60.0 && pt <= 70.0) * (0.00239520958084) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 70.0 && pt <= 80.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 80.0 && pt <= 90.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 125.0 && pt <= 150.0) * (0.00109769484083) +
          (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 150.0) * (0.00012725884449) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 26.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 26.0 && pt <= 32.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 32.0 && pt <= 38.0) * (0.00183318056829) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 38.0 && pt <= 44.0) * (0.00359281437126) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 44.0 && pt <= 50.0) * (0.00140252454418) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 60.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 60.0 && pt <= 70.0) * (0.00240963855422) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 70.0 && pt <= 80.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 80.0 && pt <= 90.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 125.0 && pt <= 150.0) * (0.00113636363636) +
          (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 150.0) * (0.000343229792346) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 20.0 && pt <= 26.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 26.0 && pt <= 32.0) * (0.00123456790123) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 32.0 && pt <= 38.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 38.0 && pt <= 44.0) * (0.00189753320683) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 50.0 && pt <= 60.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 60.0 && pt <= 70.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 70.0 && pt <= 80.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 80.0 && pt <= 90.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 100.0 && pt <= 125.0) * (0.00159744408946) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 150.0) * (0.000914494741655) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 20.0 && pt <= 26.0) * (0.0) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 26.0 && pt <= 32.0) * (0.00138888888889) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 32.0 && pt <= 38.0) * (0.00547445255474) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 38.0 && pt <= 44.0) * (0.0) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 50.0 && pt <= 60.0) * (0.0) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 60.0 && pt <= 70.0) * (0.0) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 70.0 && pt <= 80.0) * (0.0056338028169) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 80.0 && pt <= 90.0) * (0.0) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 100.0 && pt <= 125.0) * (0.00170648464164) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 150.0) * (0.000604229607251) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 20.0 && pt <= 26.0) * (0.00108813928183) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 26.0 && pt <= 32.0) * (0.0) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 32.0 && pt <= 38.0) * (0.0) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 38.0 && pt <= 44.0) * (0.0) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 44.0 && pt <= 50.0) * (0.00294117647059) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 50.0 && pt <= 60.0) * (0.0) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 60.0 && pt <= 70.0) * (0.0) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 70.0 && pt <= 80.0) * (0.0) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 80.0 && pt <= 90.0) * (0.0) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 150.0) * (0.0) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 20.0 && pt <= 26.0) * (0.00331125827815) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 26.0 && pt <= 32.0) * (0.0016583747927) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 32.0 && pt <= 38.0) * (0.00871459694989) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 38.0 && pt <= 44.0) * (0.00284090909091) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 44.0 && pt <= 50.0) * (0.00279329608939) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 50.0 && pt <= 60.0) * (0.0) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 60.0 && pt <= 70.0) * (0.0) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 70.0 && pt <= 80.0) * (0.0) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 80.0 && pt <= 90.0) * (0.00408163265306) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 100.0 && pt <= 125.0) * (0.00771208226221) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 150.0) * (0.00247218788628) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 20.0 && pt <= 26.0) * (0.0) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 26.0 && pt <= 32.0) * (0.00531914893617) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 32.0 && pt <= 38.0) * (0.00230414746544) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 38.0 && pt <= 44.0) * (0.0030303030303) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 44.0 && pt <= 50.0) * (0.00682593856655) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 50.0 && pt <= 60.0) * (0.00236406619385) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 60.0 && pt <= 70.0) * (0.0) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 70.0 && pt <= 80.0) * (0.00413223140496) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 80.0 && pt <= 90.0) * (0.0) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 150.0) * (0.00187617260788) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 20.0 && pt <= 26.0) * (0.0) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 26.0 && pt <= 32.0) * (0.00357781753131) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 32.0 && pt <= 38.0) * (0.0) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 38.0 && pt <= 44.0) * (0.0) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 44.0 && pt <= 50.0) * (0.00763358778626) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 50.0 && pt <= 60.0) * (0.00280898876404) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 60.0 && pt <= 70.0) * (0.0) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 70.0 && pt <= 80.0) * (0.0) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 80.0 && pt <= 90.0) * (0.0) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 150.0) * (0.0) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 26.0) * (0.0) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 26.0 && pt <= 32.0) * (0.0) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 32.0 && pt <= 38.0) * (0.00322580645161) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 38.0 && pt <= 44.0) * (0.0) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 60.0) * (0.00380228136882) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 60.0 && pt <= 70.0) * (0.0048309178744) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 70.0 && pt <= 80.0) * (0.0) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 80.0 && pt <= 90.0) * (0.0) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 90.0 && pt <= 100.0) * (0.00990099009901) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 125.0) * (0.00645161290323) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 150.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 10.0 && pt <= 12.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 12.0 && pt <= 14.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 14.0 && pt <= 16.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 16.0 && pt <= 18.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 18.0 && pt <= 20.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 20.0 && pt <= 26.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 26.0 && pt <= 32.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 32.0 && pt <= 38.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 38.0 && pt <= 44.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 44.0 && pt <= 50.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 50.0 && pt <= 60.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 60.0 && pt <= 70.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 70.0 && pt <= 80.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 80.0 && pt <= 90.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 90.0 && pt <= 100.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 100.0 && pt <= 125.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 150.0) * (0.0)  
  }
    {13} {

         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 20.0 && pt <= 26.0) * (0.000455996352029) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 26.0 && pt <= 32.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 32.0 && pt <= 38.0) * (0.000825763831544) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 38.0 && pt <= 44.0) * (0.000999000999001) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 50.0 && pt <= 60.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 60.0 && pt <= 70.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 70.0 && pt <= 80.0) * (0.0011135857461) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 80.0 && pt <= 90.0) * (0.00138312586445) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 90.0 && pt <= 100.0) * (0.0014598540146) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 100.0 && pt <= 125.0) * (0.000763941940413) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 150.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 20.0 && pt <= 26.0) * (0.000484261501211) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 26.0 && pt <= 32.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 32.0 && pt <= 38.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 38.0 && pt <= 44.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 50.0 && pt <= 60.0) * (0.000862068965517) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 60.0 && pt <= 70.0) * (0.00109529025192) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 70.0 && pt <= 80.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 80.0 && pt <= 90.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 150.0) * (0.000103766732386) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 20.0 && pt <= 26.0) * (0.00050994390617) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 26.0 && pt <= 32.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 32.0 && pt <= 38.0) * (0.000866551126516) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 38.0 && pt <= 44.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 50.0 && pt <= 60.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 60.0 && pt <= 70.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 70.0 && pt <= 80.0) * (0.00127877237852) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 80.0 && pt <= 90.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 150.0) * (0.000453617600363) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 20.0 && pt <= 26.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 26.0 && pt <= 32.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 32.0 && pt <= 38.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 38.0 && pt <= 44.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 50.0 && pt <= 60.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 60.0 && pt <= 70.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 70.0 && pt <= 80.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 80.0 && pt <= 90.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 150.0) * (0.00012725884449) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 26.0) * (0.000535618639529) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 26.0 && pt <= 32.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 32.0 && pt <= 38.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 38.0 && pt <= 44.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 60.0) * (0.00202634245187) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 60.0 && pt <= 70.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 70.0 && pt <= 80.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 80.0 && pt <= 90.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 150.0) * (0.000343229792346) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 20.0 && pt <= 26.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 26.0 && pt <= 32.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 32.0 && pt <= 38.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 38.0 && pt <= 44.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 50.0 && pt <= 60.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 60.0 && pt <= 70.0) * (0.00163132137031) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 70.0 && pt <= 80.0) * (0.00392156862745) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 80.0 && pt <= 90.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 125.0 && pt <= 150.0) * (0.00167504187605) +
          (abs(eta) > 1.5 && abs(eta) <= 1.76) * (pt > 150.0) * (0.00115740740741) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 20.0 && pt <= 26.0) * (0.0) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 26.0 && pt <= 32.0) * (0.0) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 32.0 && pt <= 38.0) * (0.0) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 38.0 && pt <= 44.0) * (0.0) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 50.0 && pt <= 60.0) * (0.0) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 60.0 && pt <= 70.0) * (0.0) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 70.0 && pt <= 80.0) * (0.0) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 80.0 && pt <= 90.0) * (0.0) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 100.0 && pt <= 125.0) * (0.00144927536232) +
         (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 1.76 && abs(eta) <= 2.02) * (pt > 150.0) * (0.0) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 20.0 && pt <= 26.0) * (0.0) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 26.0 && pt <= 32.0) * (0.0) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 32.0 && pt <= 38.0) * (0.0) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 38.0 && pt <= 44.0) * (0.0) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 50.0 && pt <= 60.0) * (0.0) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 60.0 && pt <= 70.0) * (0.0) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 70.0 && pt <= 80.0) * (0.0) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 80.0 && pt <= 90.0) * (0.003125) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 2.02 && abs(eta) <= 2.28) * (pt > 150.0) * (0.0) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 20.0 && pt <= 26.0) * (0.0) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 26.0 && pt <= 32.0) * (0.0) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 32.0 && pt <= 38.0) * (0.0) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 38.0 && pt <= 44.0) * (0.00252525252525) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 50.0 && pt <= 60.0) * (0.0) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 60.0 && pt <= 70.0) * (0.0) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 70.0 && pt <= 80.0) * (0.0) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 80.0 && pt <= 90.0) * (0.00375939849624) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 125.0 && pt <= 150.0) * (0.00369003690037) +
          (abs(eta) > 2.28 && abs(eta) <= 2.54) * (pt > 150.0) * (0.0) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 20.0 && pt <= 26.0) * (0.00123456790123) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 26.0 && pt <= 32.0) * (0.0) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 32.0 && pt <= 38.0) * (0.0) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 38.0 && pt <= 44.0) * (0.0) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 50.0 && pt <= 60.0) * (0.00232558139535) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 60.0 && pt <= 70.0) * (0.0) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 70.0 && pt <= 80.0) * (0.0) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 80.0 && pt <= 90.0) * (0.0) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 2.54 && abs(eta) <= 2.8) * (pt > 150.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 10.0 && pt <= 12.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 12.0 && pt <= 14.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 14.0 && pt <= 16.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 16.0 && pt <= 18.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 18.0 && pt <= 20.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 20.0 && pt <= 26.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 26.0 && pt <= 32.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 32.0 && pt <= 38.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 38.0 && pt <= 44.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 44.0 && pt <= 50.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 50.0 && pt <= 60.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 60.0 && pt <= 70.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 70.0 && pt <= 80.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 80.0 && pt <= 90.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 90.0 && pt <= 100.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 100.0 && pt <= 125.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 150.0) * (0.0)  
  }
    {22} {

         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 20.0 && pt <= 26.0) * (0.00273597811218) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 26.0 && pt <= 32.0) * (0.00394996708361) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 32.0 && pt <= 38.0) * (0.0107349298101) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 38.0 && pt <= 44.0) * (0.00599400599401) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 44.0 && pt <= 50.0) * (0.0046565774156) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 50.0 && pt <= 60.0) * (0.00261324041812) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 60.0 && pt <= 70.0) * (0.00221238938053) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 70.0 && pt <= 80.0) * (0.0022271714922) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 80.0 && pt <= 90.0) * (0.00138312586445) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 125.0 && pt <= 150.0) * (0.00269784172662) +
          (abs(eta) > 0.0 && abs(eta) <= 0.3) * (pt > 150.0) * (0.00186183243508) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 20.0 && pt <= 26.0) * (0.0043583535109) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 26.0 && pt <= 32.0) * (0.00618131868132) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 32.0 && pt <= 38.0) * (0.00512820512821) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 38.0 && pt <= 44.0) * (0.00428265524625) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 44.0 && pt <= 50.0) * (0.00492610837438) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 50.0 && pt <= 60.0) * (0.00431034482759) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 60.0 && pt <= 70.0) * (0.00219058050383) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 70.0 && pt <= 80.0) * (0.00223713646532) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 80.0 && pt <= 90.0) * (0.0015243902439) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 90.0 && pt <= 100.0) * (0.0016051364366) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 100.0 && pt <= 125.0) * (0.000735294117647) +
         (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 125.0 && pt <= 150.0) * (0.00396825396825) +
          (abs(eta) > 0.3 && abs(eta) <= 0.6) * (pt > 150.0) * (0.0021791013801) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 20.0 && pt <= 26.0) * (0.00560938296787) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 26.0 && pt <= 32.0) * (0.00999286224126) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 32.0 && pt <= 38.0) * (0.0112651646447) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 38.0 && pt <= 44.0) * (0.00756756756757) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 44.0 && pt <= 50.0) * (0.00611246943765) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 50.0 && pt <= 60.0) * (0.00587084148728) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 60.0 && pt <= 70.0) * (0.00752688172043) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 70.0 && pt <= 80.0) * (0.00255754475703) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 80.0 && pt <= 90.0) * (0.00147058823529) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 90.0 && pt <= 100.0) * (0.0016077170418) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 100.0 && pt <= 125.0) * (0.00161290322581) +
         (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 125.0 && pt <= 150.0) * (0.000979431929481) +
          (abs(eta) > 0.6 && abs(eta) <= 0.9) * (pt > 150.0) * (0.00226808800181) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 20.0 && pt <= 26.0) * (0.00900450225113) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 26.0 && pt <= 32.0) * (0.0128388017118) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 32.0 && pt <= 38.0) * (0.00944206008584) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 38.0 && pt <= 44.0) * (0.0122358175751) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 44.0 && pt <= 50.0) * (0.0132450331126) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 50.0 && pt <= 60.0) * (0.0124521072797) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 60.0 && pt <= 70.0) * (0.00479041916168) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 70.0 && pt <= 80.0) * (0.00702247191011) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 80.0 && pt <= 90.0) * (0.00485436893204) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 90.0 && pt <= 100.0) * (0.0017793594306) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 100.0 && pt <= 125.0) * (0.00416666666667) +
         (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 125.0 && pt <= 150.0) * (0.00219538968167) +
          (abs(eta) > 0.9 && abs(eta) <= 1.2) * (pt > 150.0) * (0.00229065920081) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 26.0) * (0.0208891269416) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 26.0 && pt <= 32.0) * (0.0183626625861) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 32.0 && pt <= 38.0) * (0.0348304307974) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 38.0 && pt <= 44.0) * (0.0239520958084) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 44.0 && pt <= 50.0) * (0.0266479663394) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 60.0) * (0.0253292806484) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 60.0 && pt <= 70.0) * (0.0349397590361) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 70.0 && pt <= 80.0) * (0.0171184022825) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 80.0 && pt <= 90.0) * (0.00679117147708) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 90.0 && pt <= 100.0) * (0.0323886639676) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 125.0) * (0.00833333333333) +
         (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 125.0 && pt <= 150.0) * (0.00909090909091) +
          (abs(eta) > 1.2 && abs(eta) <= 1.5) * (pt > 150.0) * (0.00772267032778) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 20.0 && pt <= 26.0) * (0.0100536193029) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 26.0 && pt <= 32.0) * (0.00727272727273) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 32.0 && pt <= 38.0) * (0.00793650793651) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 38.0 && pt <= 44.0) * (0.0042194092827) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 44.0 && pt <= 50.0) * (0.00340715502555) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60.0) * (0.0109356014581) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70.0) * (0.00579710144928) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80.0) * (0.013986013986) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90.0) * (0.00801603206413) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100.0) * (0.00666666666667) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 125.0 && pt <= 150.0) * (0.00738552437223) +
          (abs(eta) > 1.5 && abs(eta) <= 1.8) * (pt > 150.0) * (0.00034188034188) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 20.0 && pt <= 26.0) * (0.00298730395818) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 26.0 && pt <= 32.0) * (0.00309917355372) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 32.0 && pt <= 38.0) * (0.0) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 38.0 && pt <= 44.0) * (0.0) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 50.0 && pt <= 60.0) * (0.0) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 60.0 && pt <= 70.0) * (0.0) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 70.0 && pt <= 80.0) * (0.0) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 80.0 && pt <= 90.0) * (0.0) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 1.8 && abs(eta) <= 2.1) * (pt > 150.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 26.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 26.0 && pt <= 32.0) * (0.00115606936416) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 32.0 && pt <= 38.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 38.0 && pt <= 44.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 2.1 && abs(eta) <= 2.4) * (pt > 150.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 20.0 && pt <= 26.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 26.0 && pt <= 32.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 32.0 && pt <= 38.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 38.0 && pt <= 44.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 50.0 && pt <= 60.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 60.0 && pt <= 70.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 70.0 && pt <= 80.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 80.0 && pt <= 90.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 2.4 && abs(eta) <= 2.7) * (pt > 150.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 10.0 && pt <= 12.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 12.0 && pt <= 14.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 14.0 && pt <= 16.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 16.0 && pt <= 18.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 18.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 26.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 26.0 && pt <= 32.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 32.0 && pt <= 38.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 38.0 && pt <= 44.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 44.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 60.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 60.0 && pt <= 70.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 70.0 && pt <= 80.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 80.0 && pt <= 90.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 90.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 2.7 && abs(eta) <= 3.0) * (pt > 150.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 10.0 && pt <= 12.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 12.0 && pt <= 14.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 14.0 && pt <= 16.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 16.0 && pt <= 18.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 18.0 && pt <= 20.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 20.0 && pt <= 26.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 26.0 && pt <= 32.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 32.0 && pt <= 38.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 38.0 && pt <= 44.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 44.0 && pt <= 50.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 50.0 && pt <= 60.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 60.0 && pt <= 70.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 70.0 && pt <= 80.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 80.0 && pt <= 90.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 90.0 && pt <= 100.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 100.0 && pt <= 125.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 3.0) * (pt > 150.0) * (0.0)  
  }
  }

}


############################
# Photon fake merger loose
############################

module Merger PhotonFakeMergerLoose {
# add InputArray InputArray
  add InputArray PhotonLooseID/photons
  add InputArray JetFakeMakerLoose/photons
  set OutputArray photons
}

############################
# Photon fake merger medium
############################

module Merger PhotonFakeMergerMedium {
# add InputArray InputArray
  add InputArray PhotonMediumID/photons
  add InputArray JetFakeMakerMedium/photons
  set OutputArray photons
}

############################
# Photon fake merger tight
############################

module Merger PhotonFakeMergerTight {
# add InputArray InputArray
  add InputArray PhotonTightID/photons
  add InputArray JetFakeMakerTight/photons
  set OutputArray photons
}

############################
# Electron fake merger loose
############################

module Merger ElectronFakeMergerLoose {
# add InputArray InputArray
  add InputArray ElectronLooseEfficiency/electrons
  add InputArray JetFakeMakerLoose/electrons
  set OutputArray electrons
}

############################
# Electron fake merger medium
############################

module Merger ElectronFakeMergerMedium {
# add InputArray InputArray
  add InputArray ElectronMediumEfficiency/electrons
  add InputArray JetFakeMakerMedium/electrons
  set OutputArray electrons
}

############################
# Electron fake merger tight
############################

module Merger ElectronFakeMergerTight {
# add InputArray InputArray
  add InputArray ElectronTightEfficiency/electrons
  add InputArray JetFakeMakerTight/electrons
  set OutputArray electrons
}


############################
# Muon fake merger loose
############################

module Merger MuonFakeMergerLoose {
# add InputArray InputArray
  add InputArray MuonLooseIdEfficiency/muons
  add InputArray JetFakeMakerLoose/muons
  set OutputArray muons
}

############################
# Muon fake merger medium
############################

module Merger MuonFakeMergerMedium {
# add InputArray InputArray
  add InputArray MuonMediumIdEfficiency/muons
  add InputArray JetFakeMakerMedium/muons
  set OutputArray muons
}

############################
# Muon fake merger tight
############################

module Merger MuonFakeMergerTight {
# add InputArray InputArray
  add InputArray MuonMediumIdEfficiency/muons
  add InputArray JetFakeMakerTight/muons
  set OutputArray muons
}




###############################################################################################################
# StatusPidFilter: this module removes all generated particles except electrons, muons, taus, and status == 3 #
###############################################################################################################

module StatusPidFilter GenParticleFilter {

    set InputArray Delphes/allParticles
    set OutputArray filteredParticles
    set PTMin 0.0

}


####################
# ROOT tree writer
####################

module TreeWriter TreeWriter {

# add Branch InputArray BranchName BranchClass
  #add Branch GenParticleFilter/filteredParticles Particle GenParticle
  add Branch Delphes/allParticles Particle GenParticle
  add Branch PileUpMerger/vertices Vertex Vertex

  add Branch GenJetFinder/jets GenJet Jet
  add Branch GenJetFinderAK8/jetsAK8 GenJetAK8 Jet
  add Branch GenMissingET/momentum GenMissingET MissingET

#  add Branch HCal/eflowTracks EFlowTrack Track
#  add Branch ECal/eflowPhotons EFlowPhoton Tower
#  add Branch HCal/eflowNeutralHadrons EFlowNeutralHadron Tower

  add Branch RunPUPPI/PuppiParticles ParticleFlowCandidate ParticleFlowCandidate

  add Branch PhotonSmear/photons Photon Photon
  add Branch PhotonFakeMergerLoose/photons PhotonLoose Photon
  add Branch PhotonFakeMergerMedium/photons PhotonMedium Photon
  add Branch PhotonFakeMergerTight/photons PhotonTight Photon

  add Branch ElectronSmear/electrons Electron Electron
  add Branch ElectronFakeMergerLoose/electrons ElectronLoose Electron
  add Branch ElectronFakeMergerMedium/electrons ElectronMedium Electron
  add Branch ElectronFakeMergerTight/electrons ElectronTight Electron

  add Branch MuonSmear/muons Muon Muon
  add Branch MuonFakeMergerLoose/muons MuonLoose Muon
  add Branch MuonFakeMergerMedium/muons MuonMedium Muon
  add Branch MuonFakeMergerTight/muons MuonTight Muon

  add Branch JetSmearPUPPI/jets JetPUPPI Jet
  add Branch JetSmearPUPPIAK8/jets JetPUPPIAK8 Jet

  add Branch JetLooseID/jets JetLoose Jet
  add Branch JetTightID/jets JetTight Jet

  add Branch Rho/rho Rho Rho
  add Branch PuppiMissingET/momentum PuppiMissingET MissingET
  add Branch GenPileUpMissingET/momentum GenPileUpMissingET MissingET
  add Branch ScalarHT/energy ScalarHT ScalarHT

}
