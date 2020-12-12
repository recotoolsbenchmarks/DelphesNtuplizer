
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
   (abs(eta) > 3.0 && abs(eta) <= 5.0) * (pt > 200.0 && pt <= 400.0) * (0.63) +
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
   (abs(eta) > 0.0 && abs(eta) <= 1.3) * (pt > 20.0 && pt <= 50.0) * (0.063782) +
   (abs(eta) > 0.0 && abs(eta) <= 1.3) * (pt > 50.0 && pt <= 100.0) * (0.110779) +
   (abs(eta) > 0.0 && abs(eta) <= 1.3) * (pt > 100.0 && pt <= 200.0) * (0.116716) +
   (abs(eta) > 0.0 && abs(eta) <= 1.3) * (pt > 200.0 && pt <= 400.0) * (0.186809) +
   (abs(eta) > 0.0 && abs(eta) <= 1.3) * (pt > 400.0 && pt <= 14000.0) * (0.059356) +
   (abs(eta) > 1.3 && abs(eta) <= 2.5) * (pt > 20.0 && pt <= 50.0) * (0.085486) +
   (abs(eta) > 1.3 && abs(eta) <= 2.5) * (pt > 50.0 && pt <= 100.0) * (0.117623) +
   (abs(eta) > 1.3 && abs(eta) <= 2.5) * (pt > 100.0 && pt <= 200.0) * (0.152267) +
   (abs(eta) > 1.3 && abs(eta) <= 2.5) * (pt > 200.0 && pt <= 400.0) * (0.000001) +
   (abs(eta) > 1.3 && abs(eta) <= 2.5) * (pt > 400.0 && pt <= 14000.0) * (0.164762) +
   (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 50.0) * (0.075379) +
   (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 100.0) * (0.117915) +
   (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 200.0) * (0.129564) +
   (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 200.0 && pt <= 400.0) * (0.112097) +
   (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 400.0 && pt <= 14000.0) * (0.999633) +
   (abs(eta) > 3.0 && abs(eta) <= 5.0) * (pt > 20.0 && pt <= 50.0) * (0.156740) +
   (abs(eta) > 3.0 && abs(eta) <= 5.0) * (pt > 50.0 && pt <= 100.0) * (0.000001) +
   (abs(eta) > 3.0 && abs(eta) <= 5.0) * (pt > 100.0 && pt <= 200.0) * (0.997596) +
   (abs(eta) > 3.0 && abs(eta) <= 5.0) * (pt > 200.0 && pt <= 400.0) * (0.000001) +
   (abs(eta) > 3.0 && abs(eta) <= 5.0) * (pt > 400.0 && pt <= 14000.0) * (0.673551)  
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

   (abs(eta) > 0.0 && abs(eta) <= 0.65) * (pt > 20.0 && pt <= 35.0) * (0.694807130451) +
   (abs(eta) > 0.0 && abs(eta) <= 0.65) * (pt > 35.0 && pt <= 50.0) * (0.904983829541) +
   (abs(eta) > 0.0 && abs(eta) <= 0.65) * (pt > 50.0 && pt <= 75.0) * (0.947461345908) +
   (abs(eta) > 0.0 && abs(eta) <= 0.65) * (pt > 75.0 && pt <= 100.0) * (0.977742538174) +
   (abs(eta) > 0.0 && abs(eta) <= 0.65) * (pt > 100.0 && pt <= 150.0) * (0.99080735194) +
   (abs(eta) > 0.0 && abs(eta) <= 0.65) * (pt > 150.0 && pt <= 200.0) * (0.990189113115) +
   (abs(eta) > 0.0 && abs(eta) <= 0.65) * (pt > 200.0 && pt <= 300.0) * (0.998201841376) +
   (abs(eta) > 0.0 && abs(eta) <= 0.65) * (pt > 300.0 && pt <= 400.0) * (0.997569144339) +
   (abs(eta) > 0.0 && abs(eta) <= 0.65) * (pt > 400.0) * (0.998595478438) +
   (abs(eta) > 0.65 && abs(eta) <= 1.3) * (pt > 20.0 && pt <= 35.0) * (0.754043394336) +
   (abs(eta) > 0.65 && abs(eta) <= 1.3) * (pt > 35.0 && pt <= 50.0) * (0.906688346402) +
   (abs(eta) > 0.65 && abs(eta) <= 1.3) * (pt > 50.0 && pt <= 75.0) * (0.968047152663) +
   (abs(eta) > 0.65 && abs(eta) <= 1.3) * (pt > 75.0 && pt <= 100.0) * (0.982557856057) +
   (abs(eta) > 0.65 && abs(eta) <= 1.3) * (pt > 100.0 && pt <= 150.0) * (0.993551171026) +
   (abs(eta) > 0.65 && abs(eta) <= 1.3) * (pt > 150.0 && pt <= 200.0) * (0.998464420923) +
   (abs(eta) > 0.65 && abs(eta) <= 1.3) * (pt > 200.0 && pt <= 300.0) * (1.0) +
   (abs(eta) > 0.65 && abs(eta) <= 1.3) * (pt > 300.0 && pt <= 400.0) * (0.999771449037) +
   (abs(eta) > 0.65 && abs(eta) <= 1.3) * (pt > 400.0) * (1.0) +
   (abs(eta) > 1.3 && abs(eta) <= 1.9) * (pt > 20.0 && pt <= 35.0) * (0.673949197816) +
   (abs(eta) > 1.3 && abs(eta) <= 1.9) * (pt > 35.0 && pt <= 50.0) * (0.760336471446) +
   (abs(eta) > 1.3 && abs(eta) <= 1.9) * (pt > 50.0 && pt <= 75.0) * (0.827989935851) +
   (abs(eta) > 1.3 && abs(eta) <= 1.9) * (pt > 75.0 && pt <= 100.0) * (0.792151429041) +
   (abs(eta) > 1.3 && abs(eta) <= 1.9) * (pt > 100.0 && pt <= 150.0) * (0.807914341774) +
   (abs(eta) > 1.3 && abs(eta) <= 1.9) * (pt > 150.0 && pt <= 200.0) * (0.854415509259) +
   (abs(eta) > 1.3 && abs(eta) <= 1.9) * (pt > 200.0 && pt <= 300.0) * (0.798419793631) +
   (abs(eta) > 1.3 && abs(eta) <= 1.9) * (pt > 300.0 && pt <= 400.0) * (0.738030209091) +
   (abs(eta) > 1.3 && abs(eta) <= 1.9) * (pt > 400.0) * (0.715819975739) +
   (abs(eta) > 1.9 && abs(eta) <= 2.5) * (pt > 20.0 && pt <= 35.0) * (0.480518380612) +
   (abs(eta) > 1.9 && abs(eta) <= 2.5) * (pt > 35.0 && pt <= 50.0) * (0.688434612367) +
   (abs(eta) > 1.9 && abs(eta) <= 2.5) * (pt > 50.0 && pt <= 75.0) * (0.736953934646) +
   (abs(eta) > 1.9 && abs(eta) <= 2.5) * (pt > 75.0 && pt <= 100.0) * (0.741652611194) +
   (abs(eta) > 1.9 && abs(eta) <= 2.5) * (pt > 100.0 && pt <= 150.0) * (0.728457834286) +
   (abs(eta) > 1.9 && abs(eta) <= 2.5) * (pt > 150.0 && pt <= 200.0) * (0.673561444908) +
   (abs(eta) > 1.9 && abs(eta) <= 2.5) * (pt > 200.0 && pt <= 300.0) * (0.654171821892) +
   (abs(eta) > 1.9 && abs(eta) <= 2.5) * (pt > 300.0 && pt <= 400.0) * (0.612942612943) +
   (abs(eta) > 1.9 && abs(eta) <= 2.5) * (pt > 400.0) * (0.549180327869) +
   (abs(eta) > 2.5 && abs(eta) <= 2.75) * (pt > 20.0 && pt <= 35.0) * (0.65856188631) +
   (abs(eta) > 2.5 && abs(eta) <= 2.75) * (pt > 35.0 && pt <= 50.0) * (0.816545080434) +
   (abs(eta) > 2.5 && abs(eta) <= 2.75) * (pt > 50.0 && pt <= 75.0) * (0.750517675823) +
   (abs(eta) > 2.5 && abs(eta) <= 2.75) * (pt > 75.0 && pt <= 100.0) * (0.702645252716) +
   (abs(eta) > 2.5 && abs(eta) <= 2.75) * (pt > 100.0 && pt <= 150.0) * (0.73118729097) +
   (abs(eta) > 2.5 && abs(eta) <= 2.75) * (pt > 150.0 && pt <= 200.0) * (0.618518901054) +
   (abs(eta) > 2.5 && abs(eta) <= 2.75) * (pt > 200.0 && pt <= 300.0) * (0.593644781145) +
   (abs(eta) > 2.5 && abs(eta) <= 2.75) * (pt > 300.0 && pt <= 400.0) * (0.377259036145) +
   (abs(eta) > 2.5 && abs(eta) <= 2.75) * (pt > 400.0) * (0.444444444444) +
   (abs(eta) > 2.75 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 35.0) * (0.495393237082) +
   (abs(eta) > 2.75 && abs(eta) <= 3.0) * (pt > 35.0 && pt <= 50.0) * (0.810540221094) +
   (abs(eta) > 2.75 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 75.0) * (0.698979056414) +
   (abs(eta) > 2.75 && abs(eta) <= 3.0) * (pt > 75.0 && pt <= 100.0) * (0.727066174382) +
   (abs(eta) > 2.75 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 150.0) * (0.69433893602) +
   (abs(eta) > 2.75 && abs(eta) <= 3.0) * (pt > 150.0 && pt <= 200.0) * (0.64) +
   (abs(eta) > 2.75 && abs(eta) <= 3.0) * (pt > 200.0 && pt <= 300.0) * (0.501824817518) +
   (abs(eta) > 2.75 && abs(eta) <= 3.0) * (pt > 300.0 && pt <= 400.0) * (0.666666666667) +
   (abs(eta) > 2.75 && abs(eta) <= 3.0) * (pt > 400.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 20.0 && pt <= 35.0) * (0.0990122148822) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 35.0 && pt <= 50.0) * (0.336202510399) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 50.0 && pt <= 75.0) * (0.412796208531) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 75.0 && pt <= 100.0) * (0.315018219677) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 100.0 && pt <= 150.0) * (0.261485441372) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 150.0 && pt <= 200.0) * (0.277266483516) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 200.0 && pt <= 300.0) * (0.347826086957) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 300.0 && pt <= 400.0) * (0.5) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 400.0) * (1.0) +
   (abs(eta) > 4.0 && abs(eta) <= 5.0) * (pt > 20.0 && pt <= 35.0) * (0.465783972125) +
   (abs(eta) > 4.0 && abs(eta) <= 5.0) * (pt > 35.0 && pt <= 50.0) * (0.264168401573) +
   (abs(eta) > 4.0 && abs(eta) <= 5.0) * (pt > 50.0 && pt <= 75.0) * (0.279279279279) +
   (abs(eta) > 4.0 && abs(eta) <= 5.0) * (pt > 75.0 && pt <= 100.0) * (0.170168067227) +
   (abs(eta) > 4.0 && abs(eta) <= 5.0) * (pt > 100.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 4.0 && abs(eta) <= 5.0) * (pt > 150.0 && pt <= 200.0) * (1.0) +
   (abs(eta) > 4.0 && abs(eta) <= 5.0) * (pt > 200.0 && pt <= 300.0) * (1.0) +
   (abs(eta) > 4.0 && abs(eta) <= 5.0) * (pt > 300.0 && pt <= 400.0) * (1.0) +
   (abs(eta) > 4.0 && abs(eta) <= 5.0) * (pt > 400.0) * (1.0) +
   (abs(eta) > 5.0) * (pt > 20.0 && pt <= 35.0) * (1.0) +
   (abs(eta) > 5.0) * (pt > 35.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 5.0) * (pt > 50.0 && pt <= 75.0) * (1.0) +
   (abs(eta) > 5.0) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 5.0) * (pt > 100.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 5.0) * (pt > 150.0 && pt <= 200.0) * (1.0) +
   (abs(eta) > 5.0) * (pt > 200.0 && pt <= 300.0) * (1.0) +
   (abs(eta) > 5.0) * (pt > 300.0 && pt <= 400.0) * (1.0) +
   (abs(eta) > 5.0) * (pt > 400.0) * (1.0)  
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

   (abs(eta) > 0.0 && abs(eta) <= 0.65) * (pt > 20.0 && pt <= 35.0) * (0.694807130451) +
   (abs(eta) > 0.0 && abs(eta) <= 0.65) * (pt > 35.0 && pt <= 50.0) * (0.903033433356) +
   (abs(eta) > 0.0 && abs(eta) <= 0.65) * (pt > 50.0 && pt <= 75.0) * (0.947461345908) +
   (abs(eta) > 0.0 && abs(eta) <= 0.65) * (pt > 75.0 && pt <= 100.0) * (0.977742538174) +
   (abs(eta) > 0.0 && abs(eta) <= 0.65) * (pt > 100.0 && pt <= 150.0) * (0.989015656909) +
   (abs(eta) > 0.0 && abs(eta) <= 0.65) * (pt > 150.0 && pt <= 200.0) * (0.990189113115) +
   (abs(eta) > 0.0 && abs(eta) <= 0.65) * (pt > 200.0 && pt <= 300.0) * (0.996267341684) +
   (abs(eta) > 0.0 && abs(eta) <= 0.65) * (pt > 300.0 && pt <= 400.0) * (0.997569144339) +
   (abs(eta) > 0.0 && abs(eta) <= 0.65) * (pt > 400.0) * (0.995550980028) +
   (abs(eta) > 0.65 && abs(eta) <= 1.3) * (pt > 20.0 && pt <= 35.0) * (0.754043394336) +
   (abs(eta) > 0.65 && abs(eta) <= 1.3) * (pt > 35.0 && pt <= 50.0) * (0.904444068317) +
   (abs(eta) > 0.65 && abs(eta) <= 1.3) * (pt > 50.0 && pt <= 75.0) * (0.968047152663) +
   (abs(eta) > 0.65 && abs(eta) <= 1.3) * (pt > 75.0 && pt <= 100.0) * (0.982557856057) +
   (abs(eta) > 0.65 && abs(eta) <= 1.3) * (pt > 100.0 && pt <= 150.0) * (0.993551171026) +
   (abs(eta) > 0.65 && abs(eta) <= 1.3) * (pt > 150.0 && pt <= 200.0) * (0.998464420923) +
   (abs(eta) > 0.65 && abs(eta) <= 1.3) * (pt > 200.0 && pt <= 300.0) * (0.99706517903) +
   (abs(eta) > 0.65 && abs(eta) <= 1.3) * (pt > 300.0 && pt <= 400.0) * (0.999771449037) +
   (abs(eta) > 0.65 && abs(eta) <= 1.3) * (pt > 400.0) * (1.0) +
   (abs(eta) > 1.3 && abs(eta) <= 1.9) * (pt > 20.0 && pt <= 35.0) * (0.673949197816) +
   (abs(eta) > 1.3 && abs(eta) <= 1.9) * (pt > 35.0 && pt <= 50.0) * (0.757400809394) +
   (abs(eta) > 1.3 && abs(eta) <= 1.9) * (pt > 50.0 && pt <= 75.0) * (0.827989935851) +
   (abs(eta) > 1.3 && abs(eta) <= 1.9) * (pt > 75.0 && pt <= 100.0) * (0.78800403936) +
   (abs(eta) > 1.3 && abs(eta) <= 1.9) * (pt > 100.0 && pt <= 150.0) * (0.805230240971) +
   (abs(eta) > 1.3 && abs(eta) <= 1.9) * (pt > 150.0 && pt <= 200.0) * (0.854415509259) +
   (abs(eta) > 1.3 && abs(eta) <= 1.9) * (pt > 200.0 && pt <= 300.0) * (0.798419793631) +
   (abs(eta) > 1.3 && abs(eta) <= 1.9) * (pt > 300.0 && pt <= 400.0) * (0.738030209091) +
   (abs(eta) > 1.3 && abs(eta) <= 1.9) * (pt > 400.0) * (0.713178573614) +
   (abs(eta) > 1.9 && abs(eta) <= 2.5) * (pt > 20.0 && pt <= 35.0) * (0.480518380612) +
   (abs(eta) > 1.9 && abs(eta) <= 2.5) * (pt > 35.0 && pt <= 50.0) * (0.688434612367) +
   (abs(eta) > 1.9 && abs(eta) <= 2.5) * (pt > 50.0 && pt <= 75.0) * (0.736953934646) +
   (abs(eta) > 1.9 && abs(eta) <= 2.5) * (pt > 75.0 && pt <= 100.0) * (0.741652611194) +
   (abs(eta) > 1.9 && abs(eta) <= 2.5) * (pt > 100.0 && pt <= 150.0) * (0.728457834286) +
   (abs(eta) > 1.9 && abs(eta) <= 2.5) * (pt > 150.0 && pt <= 200.0) * (0.673561444908) +
   (abs(eta) > 1.9 && abs(eta) <= 2.5) * (pt > 200.0 && pt <= 300.0) * (0.649020862665) +
   (abs(eta) > 1.9 && abs(eta) <= 2.5) * (pt > 300.0 && pt <= 400.0) * (0.612942612943) +
   (abs(eta) > 1.9 && abs(eta) <= 2.5) * (pt > 400.0) * (0.549180327869) +
   (abs(eta) > 2.5 && abs(eta) <= 2.75) * (pt > 20.0 && pt <= 35.0) * (0.65856188631) +
   (abs(eta) > 2.5 && abs(eta) <= 2.75) * (pt > 35.0 && pt <= 50.0) * (0.816545080434) +
   (abs(eta) > 2.5 && abs(eta) <= 2.75) * (pt > 50.0 && pt <= 75.0) * (0.750517675823) +
   (abs(eta) > 2.5 && abs(eta) <= 2.75) * (pt > 75.0 && pt <= 100.0) * (0.702645252716) +
   (abs(eta) > 2.5 && abs(eta) <= 2.75) * (pt > 100.0 && pt <= 150.0) * (0.73118729097) +
   (abs(eta) > 2.5 && abs(eta) <= 2.75) * (pt > 150.0 && pt <= 200.0) * (0.618518901054) +
   (abs(eta) > 2.5 && abs(eta) <= 2.75) * (pt > 200.0 && pt <= 300.0) * (0.593644781145) +
   (abs(eta) > 2.5 && abs(eta) <= 2.75) * (pt > 300.0 && pt <= 400.0) * (0.377259036145) +
   (abs(eta) > 2.5 && abs(eta) <= 2.75) * (pt > 400.0) * (0.444444444444) +
   (abs(eta) > 2.75 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 35.0) * (0.495393237082) +
   (abs(eta) > 2.75 && abs(eta) <= 3.0) * (pt > 35.0 && pt <= 50.0) * (0.810540221094) +
   (abs(eta) > 2.75 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 75.0) * (0.698979056414) +
   (abs(eta) > 2.75 && abs(eta) <= 3.0) * (pt > 75.0 && pt <= 100.0) * (0.727066174382) +
   (abs(eta) > 2.75 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 150.0) * (0.69433893602) +
   (abs(eta) > 2.75 && abs(eta) <= 3.0) * (pt > 150.0 && pt <= 200.0) * (0.64) +
   (abs(eta) > 2.75 && abs(eta) <= 3.0) * (pt > 200.0 && pt <= 300.0) * (0.501824817518) +
   (abs(eta) > 2.75 && abs(eta) <= 3.0) * (pt > 300.0 && pt <= 400.0) * (0.666666666667) +
   (abs(eta) > 2.75 && abs(eta) <= 3.0) * (pt > 400.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 20.0 && pt <= 35.0) * (0.0958182724667) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 35.0 && pt <= 50.0) * (0.330949346174) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 50.0 && pt <= 75.0) * (0.412796208531) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 75.0 && pt <= 100.0) * (0.315018219677) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 100.0 && pt <= 150.0) * (0.261485441372) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 150.0 && pt <= 200.0) * (0.277266483516) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 200.0 && pt <= 300.0) * (0.347826086957) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 300.0 && pt <= 400.0) * (0.5) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 400.0) * (1.0) +
   (abs(eta) > 4.0 && abs(eta) <= 5.0) * (pt > 20.0 && pt <= 35.0) * (0.326048780488) +
   (abs(eta) > 4.0 && abs(eta) <= 5.0) * (pt > 35.0 && pt <= 50.0) * (0.132084200786) +
   (abs(eta) > 4.0 && abs(eta) <= 5.0) * (pt > 50.0 && pt <= 75.0) * (0.0930930930931) +
   (abs(eta) > 4.0 && abs(eta) <= 5.0) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 4.0 && abs(eta) <= 5.0) * (pt > 100.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 4.0 && abs(eta) <= 5.0) * (pt > 150.0 && pt <= 200.0) * (1.0) +
   (abs(eta) > 4.0 && abs(eta) <= 5.0) * (pt > 200.0 && pt <= 300.0) * (1.0) +
   (abs(eta) > 4.0 && abs(eta) <= 5.0) * (pt > 300.0 && pt <= 400.0) * (1.0) +
   (abs(eta) > 4.0 && abs(eta) <= 5.0) * (pt > 400.0) * (1.0) +
   (abs(eta) > 5.0) * (pt > 20.0 && pt <= 35.0) * (1.0) +
   (abs(eta) > 5.0) * (pt > 35.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 5.0) * (pt > 50.0 && pt <= 75.0) * (1.0) +
   (abs(eta) > 5.0) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 5.0) * (pt > 100.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 5.0) * (pt > 150.0 && pt <= 200.0) * (1.0) +
   (abs(eta) > 5.0) * (pt > 200.0 && pt <= 300.0) * (1.0) +
   (abs(eta) > 5.0) * (pt > 300.0 && pt <= 400.0) * (1.0) +
   (abs(eta) > 5.0) * (pt > 400.0) * (1.0)  
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
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 20.0) * (0.044573) +
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 50.0) * (0.000001) +
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 100.0) * (0.007629) +
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 150.0) * (0.006910) +
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 150.0 && pt <= 14000.0) * (0.0097990) +
   (abs(eta) > 1.5 && abs(eta) <= 3.0) * (pt > 10.0 && pt <= 20.0) * (0.087093) +
   (abs(eta) > 1.5 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 50.0) * (0.000001) +
   (abs(eta) > 1.5 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 100.0) * (0.000001) +
   (abs(eta) > 1.5 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 150.0) * (0.013357) +
   (abs(eta) > 1.5 && abs(eta) <= 3.0) * (pt > 150.0 && pt <= 14000.0) * (0.039037)  
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

   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 10.0 && pt <= 15.0) * (0.818126979339) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 15.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 20.0 && pt <= 35.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 35.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 50.0 && pt <= 75.0) * (0.997521701182) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 100.0 && pt <= 125.0) * (0.972638101565) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 125.0 && pt <= 150.0) * (0.976663244399) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 150.0) * (0.974520231134) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 15.0) * (1.0) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 15.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 35.0) * (1.0) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 35.0 && pt <= 50.0) * (0.995749432389) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 75.0) * (0.979181148326) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 75.0 && pt <= 100.0) * (0.980673551863) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 125.0) * (0.942359109329) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 125.0 && pt <= 150.0) * (0.963418741101) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 150.0) * (0.962511280388) +
   (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 10.0 && pt <= 15.0) * (0.071384397915) +
   (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 15.0 && pt <= 20.0) * (0.25160697888) +
   (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 20.0 && pt <= 35.0) * (0.228459742657) +
   (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 35.0 && pt <= 50.0) * (0.167420531527) +
   (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 50.0 && pt <= 75.0) * (0.115267627568) +
   (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 75.0 && pt <= 100.0) * (0.0481554123471) +
   (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 100.0 && pt <= 125.0) * (0.0643756698821) +
   (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 125.0 && pt <= 150.0) * (0.0724028446801) +
   (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 150.0) * (0.0453412461048) +
   (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 10.0 && pt <= 15.0) * (1.0) +
   (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 15.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 35.0) * (0.0221088435374) +
   (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 35.0 && pt <= 50.0) * (0.0118148148148) +
   (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 75.0) * (0.0179496233522) +
   (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 75.0 && pt <= 100.0) * (0.0145833333333) +
   (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 150.0) * (0.021631006006) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 10.0 && pt <= 15.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 15.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 20.0 && pt <= 35.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 35.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 50.0 && pt <= 75.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 150.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 10.0 && pt <= 15.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 15.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 20.0 && pt <= 35.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 35.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 50.0 && pt <= 75.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 150.0) * (1.0)  
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

   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 10.0 && pt <= 15.0) * (0.79843225084) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 15.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 20.0 && pt <= 35.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 35.0 && pt <= 50.0) * (0.990918320626) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 50.0 && pt <= 75.0) * (0.97987289514) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 75.0 && pt <= 100.0) * (0.980245407434) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 100.0 && pt <= 125.0) * (0.964680263729) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 125.0 && pt <= 150.0) * (0.943163006478) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 150.0) * (0.965657772174) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 15.0) * (1.0) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 15.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 35.0) * (1.0) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 35.0 && pt <= 50.0) * (0.984185616338) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 75.0) * (0.972985856584) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 75.0 && pt <= 100.0) * (0.97153927333) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 125.0) * (0.929718769749) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 125.0 && pt <= 150.0) * (0.941000705827) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 150.0) * (0.951830808488) +
   (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 10.0 && pt <= 15.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 15.0 && pt <= 20.0) * (0.209510155317) +
   (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 20.0 && pt <= 35.0) * (0.165117878691) +
   (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 35.0 && pt <= 50.0) * (0.111216023041) +
   (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 50.0 && pt <= 75.0) * (0.0777500574467) +
   (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 75.0 && pt <= 100.0) * (0.0322720970917) +
   (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 100.0 && pt <= 125.0) * (0.0548538752128) +
   (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 125.0 && pt <= 150.0) * (0.0437511206742) +
   (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 150.0) * (0.024721867248) +
   (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 10.0 && pt <= 15.0) * (1.0) +
   (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 15.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 35.0) * (0.0224202356999) +
   (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 35.0 && pt <= 50.0) * (0.0120150659134) +
   (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 75.0) * (0.00908522236454) +
   (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 150.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 10.0 && pt <= 15.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 15.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 20.0 && pt <= 35.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 35.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 50.0 && pt <= 75.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 150.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 10.0 && pt <= 15.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 15.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 20.0 && pt <= 35.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 35.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 50.0 && pt <= 75.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 150.0) * (1.0)  
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

   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 10.0 && pt <= 15.0) * (0.653730285309) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 15.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 20.0 && pt <= 35.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 35.0 && pt <= 50.0) * (0.959562411547) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 50.0 && pt <= 75.0) * (0.929533894018) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 75.0 && pt <= 100.0) * (0.914386078962) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 100.0 && pt <= 125.0) * (0.894709818211) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 125.0 && pt <= 150.0) * (0.915110630455) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 150.0) * (0.911341265491) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 15.0) * (1.0) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 15.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 35.0) * (1.0) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 35.0 && pt <= 50.0) * (0.985152109668) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 75.0) * (0.935784126498) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 75.0 && pt <= 100.0) * (0.90712506685) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 125.0) * (0.883476140346) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 125.0 && pt <= 150.0) * (0.911468833329) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 150.0) * (0.896381642915) +
   (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 10.0 && pt <= 15.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 15.0 && pt <= 20.0) * (0.0281241981011) +
   (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 20.0 && pt <= 35.0) * (0.0713389617875) +
   (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 35.0 && pt <= 50.0) * (0.0548597906148) +
   (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 50.0 && pt <= 75.0) * (0.0363946861118) +
   (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 75.0 && pt <= 100.0) * (0.0123107778998) +
   (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 100.0 && pt <= 125.0) * (0.0351727527502) +
   (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 125.0 && pt <= 150.0) * (0.0292507492507) +
   (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 150.0) * (0.00830263105992) +
   (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 10.0 && pt <= 15.0) * (1.0) +
   (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 15.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 35.0) * (1.0) +
   (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 35.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 75.0) * (1.0) +
   (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 150.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 10.0 && pt <= 15.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 15.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 20.0 && pt <= 35.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 35.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 50.0 && pt <= 75.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 150.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 10.0 && pt <= 15.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 15.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 20.0 && pt <= 35.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 35.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 50.0 && pt <= 75.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 150.0) * (1.0)  
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
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 150.0) * (0.99) +
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
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 20.0) * (0.000001) +
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 50.0) * (0.000001) +
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 100.0) * (0.000001) +
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 150.0) * (0.013302) +
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 150.0 && pt <= 14000.0) * (0.017425) +
   (abs(eta) > 1.5 && abs(eta) <= 3.0) * (pt > 10.0 && pt <= 20.0) * (0.000001) +
   (abs(eta) > 1.5 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 50.0) * (0.000001) +
   (abs(eta) > 1.5 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 100.0) * (0.000001) +
   (abs(eta) > 1.5 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 150.0) * (0.009484) +
   (abs(eta) > 1.5 && abs(eta) <= 3.0) * (pt > 150.0 && pt <= 14000.0) * (0.000001)  
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

   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 10.0 && pt <= 15.0) * (0.872727272727) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 15.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 20.0 && pt <= 35.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 35.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 50.0 && pt <= 75.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 100.0 && pt <= 125.0) * (0.851254480287) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 150.0) * (1.0) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 15.0) * (0.623094454983) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 15.0 && pt <= 20.0) * (0.96191794102) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 35.0) * (0.917443588704) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 35.0 && pt <= 50.0) * (0.917551802094) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 75.0) * (0.880114860014) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 150.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 10.0 && pt <= 15.0) * (0.273201856148) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 15.0 && pt <= 20.0) * (0.414117647059) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 20.0 && pt <= 35.0) * (0.680297987003) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 35.0 && pt <= 50.0) * (0.756365809922) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 50.0 && pt <= 75.0) * (0.767734420501) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 100.0 && pt <= 125.0) * (0.5) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 150.0) * (1.0) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 10.0 && pt <= 15.0) * (0.399357171555) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 15.0 && pt <= 20.0) * (0.779381443299) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 20.0 && pt <= 35.0) * (0.809345939847) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 35.0 && pt <= 50.0) * (0.768261676711) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 50.0 && pt <= 75.0) * (0.951492537313) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 100.0 && pt <= 125.0) * (0.5) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 150.0) * (1.0) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 10.0 && pt <= 15.0) * (0.412177985948) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 15.0 && pt <= 20.0) * (0.576397515528) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 20.0 && pt <= 35.0) * (0.801927437642) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 35.0 && pt <= 50.0) * (0.780804889418) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 50.0 && pt <= 75.0) * (0.872368421053) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 75.0 && pt <= 100.0) * (0.666666666667) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 150.0) * (1.0) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 10.0 && pt <= 15.0) * (0.914765906363) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 15.0 && pt <= 20.0) * (0.78337236534) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 20.0 && pt <= 35.0) * (0.832051214034) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 35.0 && pt <= 50.0) * (0.852301563481) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 50.0 && pt <= 75.0) * (0.933333333333) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 150.0) * (1.0) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 10.0 && pt <= 15.0) * (0.937704918033) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 15.0 && pt <= 20.0) * (0.840082766775) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 20.0 && pt <= 35.0) * (0.89237092913) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 35.0 && pt <= 50.0) * (0.954000291342) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 50.0 && pt <= 75.0) * (0.939695550351) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 150.0) * (1.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 10.0 && pt <= 15.0) * (0.983788058521) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 15.0 && pt <= 20.0) * (0.836786684783) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 20.0 && pt <= 35.0) * (0.925430613513) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 35.0 && pt <= 50.0) * (0.973871818382) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 50.0 && pt <= 75.0) * (1.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 150.0) * (1.0) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 10.0 && pt <= 15.0) * (0.67380952381) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 15.0 && pt <= 20.0) * (0.970629370629) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 35.0) * (0.953695110634) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 35.0 && pt <= 50.0) * (0.980806953463) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 75.0) * (0.862637362637) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 150.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 10.0 && pt <= 15.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 15.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 20.0 && pt <= 35.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 35.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 50.0 && pt <= 75.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 150.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 10.0 && pt <= 15.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 15.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 20.0 && pt <= 35.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 35.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 50.0 && pt <= 75.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 150.0) * (1.0)  
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

   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 10.0 && pt <= 15.0) * (0.620622568093) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 15.0 && pt <= 20.0) * (0.968460673185) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 20.0 && pt <= 35.0) * (0.902061326655) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 35.0 && pt <= 50.0) * (0.964965343415) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 50.0 && pt <= 75.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 75.0 && pt <= 100.0) * (0.970779220779) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 100.0 && pt <= 125.0) * (0.971254480287) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 150.0) * (1.0) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 15.0) * (0.368701007839) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 15.0 && pt <= 20.0) * (0.859996460291) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 35.0) * (0.82942959493) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 35.0 && pt <= 50.0) * (0.878894802527) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 75.0) * (0.8031068688) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 75.0 && pt <= 100.0) * (0.793927327028) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 150.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 10.0 && pt <= 15.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 15.0 && pt <= 20.0) * (0.3168) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 20.0 && pt <= 35.0) * (0.475502008032) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 35.0 && pt <= 50.0) * (0.608137432188) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 50.0 && pt <= 75.0) * (0.650910364146) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 100.0 && pt <= 125.0) * (0.5) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 150.0) * (1.0) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 10.0 && pt <= 15.0) * (0.231742146063) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 15.0 && pt <= 20.0) * (0.732558139535) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 20.0 && pt <= 35.0) * (0.689974293059) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 35.0 && pt <= 50.0) * (0.647570197311) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 50.0 && pt <= 75.0) * (0.953271028037) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 75.0 && pt <= 100.0) * (0.506944444444) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 100.0 && pt <= 125.0) * (0.5) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 150.0) * (1.0) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 10.0 && pt <= 15.0) * (0.324191567745) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 15.0 && pt <= 20.0) * (0.499804151978) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 20.0 && pt <= 35.0) * (0.74995426132) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 35.0 && pt <= 50.0) * (0.680467202002) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 50.0 && pt <= 75.0) * (0.878145695364) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 75.0 && pt <= 100.0) * (0.675799086758) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 150.0) * (1.0) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 10.0 && pt <= 15.0) * (0.805213103205) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 15.0 && pt <= 20.0) * (0.702915681639) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 20.0 && pt <= 35.0) * (0.795428180894) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 35.0 && pt <= 50.0) * (0.810778598974) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 50.0 && pt <= 75.0) * (0.868860759494) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 150.0) * (1.0) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 10.0 && pt <= 15.0) * (0.702702702703) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 15.0 && pt <= 20.0) * (0.744385026738) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 20.0 && pt <= 35.0) * (0.857488431759) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 35.0 && pt <= 50.0) * (0.927244304766) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 50.0 && pt <= 75.0) * (0.877049180328) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 150.0) * (1.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 10.0 && pt <= 15.0) * (0.510045100451) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 15.0 && pt <= 20.0) * (0.807528409091) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 20.0 && pt <= 35.0) * (0.947815313267) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 35.0 && pt <= 50.0) * (0.980604759349) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 50.0 && pt <= 75.0) * (1.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 150.0) * (1.0) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 10.0 && pt <= 15.0) * (0.308027210884) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 15.0 && pt <= 20.0) * (0.921782762692) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 35.0) * (0.944184471013) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 35.0 && pt <= 50.0) * (0.996443787795) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 75.0) * (0.862637362637) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 150.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 10.0 && pt <= 15.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 15.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 20.0 && pt <= 35.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 35.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 50.0 && pt <= 75.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 150.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 10.0 && pt <= 15.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 15.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 20.0 && pt <= 35.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 35.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 50.0 && pt <= 75.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 150.0) * (1.0)  
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

   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 10.0 && pt <= 15.0) * (0.181124219292) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 15.0 && pt <= 20.0) * (0.321201429831) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 20.0 && pt <= 35.0) * (0.456118383608) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 35.0 && pt <= 50.0) * (0.679247457906) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 50.0 && pt <= 75.0) * (0.735178731455) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 75.0 && pt <= 100.0) * (0.762755102041) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 100.0 && pt <= 125.0) * (0.681003584229) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 150.0) * (1.0) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 15.0) * (0.130089492088) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 15.0 && pt <= 20.0) * (0.398260705333) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 35.0) * (0.420344378527) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 35.0 && pt <= 50.0) * (0.609265916129) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 75.0) * (0.47000191681) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 75.0 && pt <= 100.0) * (0.510563380282) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 125.0) * (0.506172839506) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 150.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 10.0 && pt <= 15.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 15.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 20.0 && pt <= 35.0) * (0.0626693766938) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 35.0 && pt <= 50.0) * (0.148760330579) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 50.0 && pt <= 75.0) * (0.241486068111) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 150.0) * (1.0) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 10.0 && pt <= 15.0) * (0.0614275414564) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 15.0 && pt <= 20.0) * (0.139689578714) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 20.0 && pt <= 35.0) * (0.234364727153) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 35.0 && pt <= 50.0) * (0.210624721975) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 50.0 && pt <= 75.0) * (0.342857142857) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 75.0 && pt <= 100.0) * (0.506944444444) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 150.0) * (1.0) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 10.0 && pt <= 15.0) * (0.115687992989) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 15.0 && pt <= 20.0) * (0.141076314989) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 20.0 && pt <= 35.0) * (0.245497931912) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 35.0 && pt <= 50.0) * (0.240909674933) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 50.0 && pt <= 75.0) * (0.340757238307) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 150.0) * (1.0) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 10.0 && pt <= 15.0) * (0.425539836188) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 15.0 && pt <= 20.0) * (0.092684954281) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 20.0 && pt <= 35.0) * (0.462075999093) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 35.0 && pt <= 50.0) * (0.54372014571) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 50.0 && pt <= 75.0) * (0.673469387755) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 75.0 && pt <= 100.0) * (0.675) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 150.0) * (1.0) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 10.0 && pt <= 15.0) * (0.279569892473) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 15.0 && pt <= 20.0) * (0.261724415794) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 20.0 && pt <= 35.0) * (0.447675086108) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 35.0 && pt <= 50.0) * (0.721723790323) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 50.0 && pt <= 75.0) * (0.692352941176) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 75.0 && pt <= 100.0) * (0.666666666667) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 150.0) * (1.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 10.0 && pt <= 15.0) * (1.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 15.0 && pt <= 20.0) * (0.424253731343) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 20.0 && pt <= 35.0) * (0.639762379055) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 35.0 && pt <= 50.0) * (0.708891101523) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 50.0 && pt <= 75.0) * (0.772479564033) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 150.0) * (1.0) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 10.0 && pt <= 15.0) * (0.0816738816739) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 15.0 && pt <= 20.0) * (0.42628992629) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 35.0) * (0.605955540146) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 35.0 && pt <= 50.0) * (0.765297774869) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 75.0) * (0.888260254597) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 150.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 10.0 && pt <= 15.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 15.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 20.0 && pt <= 35.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 35.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 50.0 && pt <= 75.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 150.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 10.0 && pt <= 15.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 15.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 20.0 && pt <= 35.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 35.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 50.0 && pt <= 75.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 4.0) * (pt > 150.0) * (1.0)  
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
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 20.0) * (1.000) +
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 50.0) * (1.000) +
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 100.0) * (1.000) +
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 150.0) * (0.985) +
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 150.0 && pt <= 14000.0) * (0.985) +
   (abs(eta) > 1.5 && abs(eta) <= 2.8) * (pt > 10.0 && pt <= 20.0) * (0.996) +
   (abs(eta) > 1.5 && abs(eta) <= 2.8) * (pt > 20.0 && pt <= 50.0) * (1.002) +
   (abs(eta) > 1.5 && abs(eta) <= 2.8) * (pt > 50.0 && pt <= 100.0) * (0.996) +
   (abs(eta) > 1.5 && abs(eta) <= 2.8) * (pt > 100.0 && pt <= 150.0) * (0.992) +
   (abs(eta) > 1.5 && abs(eta) <= 2.8) * (pt > 150.0 && pt <= 14000.0) * (0.992)  
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
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 100.0) * (0.000001) +
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 150.0) * (0.017425) +
   (abs(eta) > 0.0 && abs(eta) <= 1.5) * (pt > 150.0 && pt <= 14000.0) * (0.017425) +
   (abs(eta) > 1.5 && abs(eta) <= 2.8) * (pt > 10.0 && pt <= 20.0) * (0.008370) +
   (abs(eta) > 1.5 && abs(eta) <= 2.8) * (pt > 20.0 && pt <= 50.0) * (0.000001) +
   (abs(eta) > 1.5 && abs(eta) <= 2.8) * (pt > 50.0 && pt <= 100.0) * (0.008728) +
   (abs(eta) > 1.5 && abs(eta) <= 2.8) * (pt > 100.0 && pt <= 150.0) * (0.012709) +
   (abs(eta) > 1.5 && abs(eta) <= 2.8) * (pt > 150.0 && pt <= 14000.0) * (0.029131)  
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

   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 10.0 && pt <= 15.0) * (0.989461005721) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 15.0 && pt <= 20.0) * (0.977777777778) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 20.0 && pt <= 35.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 35.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 50.0 && pt <= 75.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 150.0) * (1.0) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 15.0) * (0.862651975684) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 15.0 && pt <= 20.0) * (0.983391608392) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 35.0) * (0.996004560628) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 35.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 75.0) * (1.0) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 150.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 10.0 && pt <= 15.0) * (0.837531107869) +
   (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 15.0 && pt <= 20.0) * (0.981856573082) +
   (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 20.0 && pt <= 35.0) * (0.97260740866) +
   (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 35.0 && pt <= 50.0) * (0.996758406071) +
   (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 50.0 && pt <= 75.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 150.0) * (1.0) +
   (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 10.0 && pt <= 15.0) * (0.960046222736) +
   (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 15.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 20.0 && pt <= 35.0) * (0.991392176372) +
   (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 35.0 && pt <= 50.0) * (0.998071093019) +
   (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 50.0 && pt <= 75.0) * (0.951228928066) +
   (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 150.0) * (1.0) +
   (abs(eta) > 2.8) * (pt > 10.0 && pt <= 15.0) * (1.0) +
   (abs(eta) > 2.8) * (pt > 15.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 2.8) * (pt > 20.0 && pt <= 35.0) * (1.0) +
   (abs(eta) > 2.8) * (pt > 35.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 2.8) * (pt > 50.0 && pt <= 75.0) * (1.0) +
   (abs(eta) > 2.8) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 2.8) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 2.8) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 2.8) * (pt > 150.0) * (1.0)  
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

   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 10.0 && pt <= 15.0) * (0.875453342271) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 15.0 && pt <= 20.0) * (0.977382875606) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 20.0 && pt <= 35.0) * (0.994739311783) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 35.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 50.0 && pt <= 75.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 150.0) * (1.0) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 15.0) * (0.787021791768) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 15.0 && pt <= 20.0) * (0.981881022204) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 35.0) * (0.989312260631) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 35.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 75.0) * (1.0) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 150.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 10.0 && pt <= 15.0) * (0.823547327431) +
   (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 15.0 && pt <= 20.0) * (0.991432068543) +
   (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 20.0 && pt <= 35.0) * (0.979399591211) +
   (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 35.0 && pt <= 50.0) * (0.988652197279) +
   (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 50.0 && pt <= 75.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 150.0) * (1.0) +
   (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 10.0 && pt <= 15.0) * (0.885950413223) +
   (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 15.0 && pt <= 20.0) * (0.964591391142) +
   (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 20.0 && pt <= 35.0) * (0.986670275837) +
   (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 35.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 50.0 && pt <= 75.0) * (0.95512102679) +
   (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 150.0) * (1.0) +
   (abs(eta) > 2.8) * (pt > 10.0 && pt <= 15.0) * (1.0) +
   (abs(eta) > 2.8) * (pt > 15.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 2.8) * (pt > 20.0 && pt <= 35.0) * (1.0) +
   (abs(eta) > 2.8) * (pt > 35.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 2.8) * (pt > 50.0 && pt <= 75.0) * (1.0) +
   (abs(eta) > 2.8) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 2.8) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 2.8) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 2.8) * (pt > 150.0) * (1.0)  
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

   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 10.0 && pt <= 15.0) * (0.317582765165) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 15.0 && pt <= 20.0) * (0.454045129334) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 20.0 && pt <= 35.0) * (0.922017030335) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 35.0 && pt <= 50.0) * (0.939245527402) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 50.0 && pt <= 75.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 75.0 && pt <= 100.0) * (0.890681003584) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 150.0) * (1.0) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 15.0) * (0.475858636137) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 15.0 && pt <= 20.0) * (0.464233812369) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 35.0) * (0.893488157989) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 35.0 && pt <= 50.0) * (0.963299388887) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 75.0) * (0.954442412466) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 75.0 && pt <= 100.0) * (0.936503496503) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 125.0) * (0.75) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 150.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 10.0 && pt <= 15.0) * (0.342249916546) +
   (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 15.0 && pt <= 20.0) * (0.590744101633) +
   (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 20.0 && pt <= 35.0) * (0.873663845224) +
   (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 35.0 && pt <= 50.0) * (0.939061099031) +
   (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 50.0 && pt <= 75.0) * (0.948077280327) +
   (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 150.0) * (1.0) +
   (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 10.0 && pt <= 15.0) * (0.582169256298) +
   (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 15.0 && pt <= 20.0) * (0.61577672996) +
   (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 20.0 && pt <= 35.0) * (0.856233489504) +
   (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 35.0 && pt <= 50.0) * (0.96494863259) +
   (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 50.0 && pt <= 75.0) * (0.876962899051) +
   (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 150.0) * (1.0) +
   (abs(eta) > 2.8) * (pt > 10.0 && pt <= 15.0) * (1.0) +
   (abs(eta) > 2.8) * (pt > 15.0 && pt <= 20.0) * (1.0) +
   (abs(eta) > 2.8) * (pt > 20.0 && pt <= 35.0) * (1.0) +
   (abs(eta) > 2.8) * (pt > 35.0 && pt <= 50.0) * (1.0) +
   (abs(eta) > 2.8) * (pt > 50.0 && pt <= 75.0) * (1.0) +
   (abs(eta) > 2.8) * (pt > 75.0 && pt <= 100.0) * (1.0) +
   (abs(eta) > 2.8) * (pt > 100.0 && pt <= 125.0) * (1.0) +
   (abs(eta) > 2.8) * (pt > 125.0 && pt <= 150.0) * (1.0) +
   (abs(eta) > 2.8) * (pt > 150.0) * (1.0)  
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

         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 20.0 && pt <= 35.0) * (0.0212765957447) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 35.0 && pt <= 50.0) * (0.0398035230352) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 50.0 && pt <= 75.0) * (0.055546662398) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 75.0 && pt <= 100.0) * (0.0497624971726) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 100.0 && pt <= 125.0) * (0.036800486618) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 125.0 && pt <= 150.0) * (0.0279667422525) +
          (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 150.0) * (0.0116054951815) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 35.0) * (0.0212722227931) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 35.0 && pt <= 50.0) * (0.0470719051149) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 75.0) * (0.0714285714286) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 75.0 && pt <= 100.0) * (0.0622329059829) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 125.0) * (0.0589451913133) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 125.0 && pt <= 150.0) * (0.0564304461942) +
          (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 150.0) * (0.0345768374165) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 20.0 && pt <= 35.0) * (0.0434782608696) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 35.0 && pt <= 50.0) * (0.104430379747) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 50.0 && pt <= 75.0) * (0.158682634731) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 75.0 && pt <= 100.0) * (0.147225368063) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 100.0 && pt <= 125.0) * (0.14696485623) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 125.0 && pt <= 150.0) * (0.157142857143) +
          (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 150.0) * (0.121170553269) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 20.0 && pt <= 35.0) * (0.094227923039) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 35.0 && pt <= 50.0) * (0.178145087236) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 50.0 && pt <= 75.0) * (0.181587837838) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 75.0 && pt <= 100.0) * (0.184405940594) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 100.0 && pt <= 125.0) * (0.143344709898) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 125.0 && pt <= 150.0) * (0.148837209302) +
          (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 150.0) * (0.107552870091) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 20.0 && pt <= 35.0) * (0.116675546084) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 35.0 && pt <= 50.0) * (0.155321188878) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 50.0 && pt <= 75.0) * (0.189719626168) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 75.0 && pt <= 100.0) * (0.202127659574) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 100.0 && pt <= 125.0) * (0.171656686627) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 125.0 && pt <= 150.0) * (0.130030959752) +
          (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 150.0) * (0.121716287215) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 20.0 && pt <= 35.0) * (0.15019988578) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 35.0 && pt <= 50.0) * (0.2071197411) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 50.0 && pt <= 75.0) * (0.237179487179) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 75.0 && pt <= 100.0) * (0.219798657718) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 100.0 && pt <= 125.0) * (0.244215938303) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 125.0 && pt <= 150.0) * (0.234693877551) +
          (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 150.0) * (0.181705809642) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 20.0 && pt <= 35.0) * (0.153941651148) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 35.0 && pt <= 50.0) * (0.226438188494) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 50.0 && pt <= 75.0) * (0.263157894737) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 75.0 && pt <= 100.0) * (0.26833976834) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 100.0 && pt <= 125.0) * (0.227692307692) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 125.0 && pt <= 150.0) * (0.27135678392) +
          (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 150.0) * (0.249530956848) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 20.0 && pt <= 35.0) * (0.0988700564972) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 35.0 && pt <= 50.0) * (0.191374663073) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 50.0 && pt <= 75.0) * (0.204735376045) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 75.0 && pt <= 100.0) * (0.255980861244) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 100.0 && pt <= 125.0) * (0.204724409449) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 125.0 && pt <= 150.0) * (0.253521126761) +
          (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 150.0) * (0.210691823899) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 35.0) * (0.0820829655781) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 35.0 && pt <= 50.0) * (0.14768683274) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 75.0) * (0.182149362477) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 75.0 && pt <= 100.0) * (0.174545454545) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 125.0) * (0.2) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 125.0 && pt <= 150.0) * (0.18691588785) +
          (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 150.0) * (0.16) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 20.0 && pt <= 35.0) * (0.00442804428044) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 35.0 && pt <= 50.0) * (0.00824175824176) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 50.0 && pt <= 75.0) * (0.0114613180516) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 75.0 && pt <= 100.0) * (0.0157480314961) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 100.0 && pt <= 125.0) * (0.017094017094) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 125.0 && pt <= 150.0) * (0.0298507462687) +
          (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 150.0) * (0.015625) +
          (abs(eta) > 4.0) * (pt > 10.0 && pt <= 15.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 15.0 && pt <= 20.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 20.0 && pt <= 35.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 35.0 && pt <= 50.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 50.0 && pt <= 75.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 75.0 && pt <= 100.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 100.0 && pt <= 125.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 150.0) * (0.0)  
  }
    {13} {

         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 20.0 && pt <= 35.0) * (0.000858696689247) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 35.0 && pt <= 50.0) * (0.000846883468835) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 50.0 && pt <= 75.0) * (0.00144069153194) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 75.0 && pt <= 100.0) * (0.0013571590138) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 100.0 && pt <= 125.0) * (0.00152068126521) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 125.0 && pt <= 150.0) * (0.000377928949358) +
          (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 150.0) * (0.00143530859135) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 35.0) * (0.000719350529236) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 35.0 && pt <= 50.0) * (0.000370644922165) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 75.0) * (0.00125) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 75.0 && pt <= 100.0) * (0.000534188034188) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 125.0) * (0.00172354360565) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 150.0) * (0.00155902004454) +
         (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 20.0 && pt <= 35.0) * (0.000483481063658) +
         (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 35.0 && pt <= 50.0) * (0.000291460215681) +
         (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 50.0 && pt <= 75.0) * (0.000830105146652) +
         (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 75.0 && pt <= 100.0) * (0.00126315789474) +
         (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 100.0 && pt <= 125.0) * (0.000578703703704) +
         (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 125.0 && pt <= 150.0) * (0.0015923566879) +
          (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 150.0) * (0.00437723836053) +
         (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 20.0 && pt <= 35.0) * (0.000417623721027) +
         (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 35.0 && pt <= 50.0) * (0.000799680127949) +
         (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 50.0 && pt <= 75.0) * (0.00158353127474) +
         (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 75.0 && pt <= 100.0) * (0.00196721311475) +
         (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 125.0 && pt <= 150.0) * (0.00157480314961) +
          (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 150.0) * (0.00428396572827) +
          (abs(eta) > 2.8) * (pt > 10.0 && pt <= 15.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 15.0 && pt <= 20.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 20.0 && pt <= 35.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 35.0 && pt <= 50.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 50.0 && pt <= 75.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 75.0 && pt <= 100.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 100.0 && pt <= 125.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 150.0) * (0.0)  
  }
    {22} {

         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 20.0 && pt <= 35.0) * (0.010686003244) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 35.0 && pt <= 50.0) * (0.0162601626016) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 50.0 && pt <= 75.0) * (0.0126460701137) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 75.0 && pt <= 100.0) * (0.00723818140692) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 100.0 && pt <= 125.0) * (0.00456204379562) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 125.0 && pt <= 150.0) * (0.00718065003779) +
          (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 150.0) * (0.00352675825302) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 35.0) * (0.0285684924468) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 35.0 && pt <= 50.0) * (0.0400296515938) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 75.0) * (0.0433928571429) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 75.0 && pt <= 100.0) * (0.0267094017094) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 125.0) * (0.0179248534988) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 125.0 && pt <= 150.0) * (0.011811023622) +
          (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 150.0) * (0.0102449888641) +
         (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 20.0 && pt <= 35.0) * (0.0172095007823) +
         (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 35.0 && pt <= 50.0) * (0.012606122974) +
         (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 50.0 && pt <= 75.0) * (0.0120511559272) +
         (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 75.0 && pt <= 100.0) * (0.00637659414854) +
         (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 100.0 && pt <= 125.0) * (0.00470957613815) +
         (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 125.0 && pt <= 150.0) * (0.00499286733238) +
          (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 150.0) * (0.00128228613299) +
         (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 35.0) * (0.0011964107677) +
         (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 35.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 75.0) * (0.0) +
         (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 75.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 150.0) * (0.0) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 20.0 && pt <= 35.0) * (0.0) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 35.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 50.0 && pt <= 75.0) * (0.0) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 75.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 150.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 10.0 && pt <= 15.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 15.0 && pt <= 20.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 20.0 && pt <= 35.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 35.0 && pt <= 50.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 50.0 && pt <= 75.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 75.0 && pt <= 100.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 100.0 && pt <= 125.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 150.0) * (0.0)  
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

         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 20.0 && pt <= 35.0) * (0.00314855452724) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 35.0 && pt <= 50.0) * (0.00321815718157) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 50.0 && pt <= 75.0) * (0.00656315031215) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 75.0 && pt <= 100.0) * (0.00520244288622) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 100.0 && pt <= 125.0) * (0.00273722627737) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 125.0 && pt <= 150.0) * (0.00226757369615) +
          (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 150.0) * (0.000861185154808) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 35.0) * (0.00277463775563) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 35.0 && pt <= 50.0) * (0.0051890289103) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 75.0) * (0.0075) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 75.0 && pt <= 100.0) * (0.00293803418803) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 125.0) * (0.00344708721131) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 125.0 && pt <= 150.0) * (0.00524934383202) +
          (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 150.0) * (0.0032293986637) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 20.0 && pt <= 35.0) * (0.00403406544151) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 35.0 && pt <= 50.0) * (0.00632911392405) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 50.0 && pt <= 75.0) * (0.00898203592814) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 75.0 && pt <= 100.0) * (0.0135900339751) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 100.0 && pt <= 125.0) * (0.0143769968051) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 125.0 && pt <= 150.0) * (0.0142857142857) +
          (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 150.0) * (0.00960219478738) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 20.0 && pt <= 35.0) * (0.0133201776024) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 35.0 && pt <= 50.0) * (0.0128558310376) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 50.0 && pt <= 75.0) * (0.0160472972973) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 75.0 && pt <= 100.0) * (0.0136138613861) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 100.0 && pt <= 125.0) * (0.0204778156997) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 125.0 && pt <= 150.0) * (0.0046511627907) +
          (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 150.0) * (0.00906344410876) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 20.0 && pt <= 35.0) * (0.0101225359616) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 35.0 && pt <= 50.0) * (0.0153403643337) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 50.0 && pt <= 75.0) * (0.014953271028) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 75.0 && pt <= 100.0) * (0.0136778115502) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 100.0 && pt <= 125.0) * (0.0059880239521) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 125.0 && pt <= 150.0) * (0.0030959752322) +
          (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 150.0) * (0.00175131348511) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 20.0 && pt <= 35.0) * (0.0279840091376) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 35.0 && pt <= 50.0) * (0.017259978425) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 50.0 && pt <= 75.0) * (0.025641025641) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 75.0 && pt <= 100.0) * (0.0218120805369) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 100.0 && pt <= 125.0) * (0.0411311053985) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 125.0 && pt <= 150.0) * (0.0136054421769) +
          (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 150.0) * (0.0160692212608) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 20.0 && pt <= 35.0) * (0.0223463687151) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 35.0 && pt <= 50.0) * (0.0416156670747) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 50.0 && pt <= 75.0) * (0.0327485380117) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 75.0 && pt <= 100.0) * (0.019305019305) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 100.0 && pt <= 125.0) * (0.0307692307692) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 125.0 && pt <= 150.0) * (0.0150753768844) +
          (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 150.0) * (0.0168855534709) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 20.0 && pt <= 35.0) * (0.0183615819209) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 35.0 && pt <= 50.0) * (0.0309973045822) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 50.0 && pt <= 75.0) * (0.0194986072423) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 75.0 && pt <= 100.0) * (0.0287081339713) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 100.0 && pt <= 125.0) * (0.0196850393701) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 125.0 && pt <= 150.0) * (0.00704225352113) +
          (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 150.0) * (0.0345911949686) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 35.0) * (0.00794351279788) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 35.0 && pt <= 50.0) * (0.017793594306) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 75.0) * (0.0218579234973) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 75.0 && pt <= 100.0) * (0.0290909090909) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 125.0) * (0.0451612903226) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 125.0 && pt <= 150.0) * (0.00934579439252) +
          (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 150.0) * (0.0457142857143) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 20.0 && pt <= 35.0) * (0.0) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 35.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 50.0 && pt <= 75.0) * (0.00286532951289) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 75.0 && pt <= 100.0) * (0.00393700787402) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 150.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 10.0 && pt <= 15.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 15.0 && pt <= 20.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 20.0 && pt <= 35.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 35.0 && pt <= 50.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 50.0 && pt <= 75.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 75.0 && pt <= 100.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 100.0 && pt <= 125.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 150.0) * (0.0)  
  }
    {13} {

         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 20.0 && pt <= 35.0) * (0.000763285945998) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 35.0 && pt <= 50.0) * (0.000677506775068) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 50.0 && pt <= 75.0) * (0.00144069153194) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 75.0 && pt <= 100.0) * (0.00113096584483) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 100.0 && pt <= 125.0) * (0.00152068126521) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 125.0 && pt <= 150.0) * (0.000377928949358) +
          (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 150.0) * (0.00131228214066) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 35.0) * (0.000719350529236) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 35.0 && pt <= 50.0) * (0.000370644922165) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 75.0) * (0.00107142857143) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 75.0 && pt <= 100.0) * (0.000534188034188) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 125.0) * (0.00137883488452) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 150.0) * (0.00144766146993) +
         (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 20.0 && pt <= 35.0) * (0.000322320709106) +
         (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 35.0 && pt <= 50.0) * (0.000291460215681) +
         (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 50.0 && pt <= 75.0) * (0.000830105146652) +
         (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 75.0 && pt <= 100.0) * (0.000842105263158) +
         (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 100.0 && pt <= 125.0) * (0.000578703703704) +
         (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 125.0 && pt <= 150.0) * (0.0015923566879) +
          (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 150.0) * (0.00338241146041) +
         (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 20.0 && pt <= 35.0) * (0.000417623721027) +
         (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 35.0 && pt <= 50.0) * (0.000799680127949) +
         (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 50.0 && pt <= 75.0) * (0.000791765637371) +
         (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 75.0 && pt <= 100.0) * (0.00196721311475) +
         (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 125.0 && pt <= 150.0) * (0.00157480314961) +
          (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 150.0) * (0.0030599755202) +
          (abs(eta) > 2.8) * (pt > 10.0 && pt <= 15.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 15.0 && pt <= 20.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 20.0 && pt <= 35.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 35.0 && pt <= 50.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 50.0 && pt <= 75.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 75.0 && pt <= 100.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 100.0 && pt <= 125.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 150.0) * (0.0)  
  }
    {22} {

         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 20.0 && pt <= 35.0) * (0.00896860986547) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 35.0 && pt <= 50.0) * (0.0133807588076) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 50.0 && pt <= 75.0) * (0.00800384184409) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 75.0 && pt <= 100.0) * (0.00316670436553) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 100.0 && pt <= 125.0) * (0.00243309002433) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 125.0 && pt <= 150.0) * (0.00604686318972) +
          (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 150.0) * (0.00291162599959) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 35.0) * (0.0248689754393) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 35.0 && pt <= 50.0) * (0.0296515937732) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 75.0) * (0.0328571428571) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 75.0 && pt <= 100.0) * (0.0162927350427) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 125.0) * (0.0124095139607) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 125.0 && pt <= 150.0) * (0.00962379702537) +
          (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 150.0) * (0.00768374164811) +
         (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 20.0 && pt <= 35.0) * (0.010382591381) +
         (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 35.0 && pt <= 50.0) * (0.00771803447389) +
         (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 50.0 && pt <= 75.0) * (0.00737825873094) +
         (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 75.0 && pt <= 100.0) * (0.00450112528132) +
         (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 100.0 && pt <= 125.0) * (0.00209314495029) +
         (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 125.0 && pt <= 150.0) * (0.00499286733238) +
          (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 150.0) * (0.000915918666422) +
         (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 35.0) * (0.000398803589232) +
         (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 35.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 75.0) * (0.0) +
         (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 75.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 150.0) * (0.0) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 20.0 && pt <= 35.0) * (0.0) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 35.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 50.0 && pt <= 75.0) * (0.0) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 75.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 150.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 10.0 && pt <= 15.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 15.0 && pt <= 20.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 20.0 && pt <= 35.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 35.0 && pt <= 50.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 50.0 && pt <= 75.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 75.0 && pt <= 100.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 100.0 && pt <= 125.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 150.0) * (0.0)  
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

         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 20.0 && pt <= 35.0) * (0.000667875202748) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 35.0 && pt <= 50.0) * (0.000508130081301) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 50.0 && pt <= 75.0) * (0.000480230510645) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 75.0 && pt <= 100.0) * (0.0013571590138) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 150.0) * (0.000205044084478) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 35.0) * (0.000102764361319) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 35.0 && pt <= 50.0) * (0.00148257968866) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 75.0) * (0.000892857142857) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 75.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 125.0 && pt <= 150.0) * (0.00087489063867) +
          (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 150.0) * (0.000167037861915) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 20.0 && pt <= 35.0) * (0.000448229493501) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 35.0 && pt <= 50.0) * (0.000791139240506) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 50.0 && pt <= 75.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 75.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 100.0 && pt <= 125.0) * (0.00159744408946) +
         (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 1.5 && abs(eta) <= 1.71428571429) * (pt > 150.0) * (0.000914494741655) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 20.0 && pt <= 35.0) * (0.0014800197336) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 35.0 && pt <= 50.0) * (0.000918273645546) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 50.0 && pt <= 75.0) * (0.000844594594595) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 75.0 && pt <= 100.0) * (0.00123762376238) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 100.0 && pt <= 125.0) * (0.00170648464164) +
         (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 1.71428571429 && abs(eta) <= 1.92857142857) * (pt > 150.0) * (0.000604229607251) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 20.0 && pt <= 35.0) * (0.000532765050613) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 35.0 && pt <= 50.0) * (0.000958772770853) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 50.0 && pt <= 75.0) * (0.0) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 75.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 1.92857142857 && abs(eta) <= 2.14285714286) * (pt > 150.0) * (0.0) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 20.0 && pt <= 35.0) * (0.00399771559109) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 35.0 && pt <= 50.0) * (0.00323624595469) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 50.0 && pt <= 75.0) * (0.0) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 75.0 && pt <= 100.0) * (0.00167785234899) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 100.0 && pt <= 125.0) * (0.00771208226221) +
         (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 2.14285714286 && abs(eta) <= 2.35714285714) * (pt > 150.0) * (0.00247218788628) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 20.0 && pt <= 35.0) * (0.00248292985723) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 35.0 && pt <= 50.0) * (0.00367197062424) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 50.0 && pt <= 75.0) * (0.00116959064327) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 75.0 && pt <= 100.0) * (0.0019305019305) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 2.35714285714 && abs(eta) <= 2.57142857143) * (pt > 150.0) * (0.00187617260788) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 20.0 && pt <= 35.0) * (0.00141242937853) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 35.0 && pt <= 50.0) * (0.00269541778976) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 50.0 && pt <= 75.0) * (0.00139275766017) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 75.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 2.57142857143 && abs(eta) <= 2.78571428571) * (pt > 150.0) * (0.0) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 35.0) * (0.0) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 35.0 && pt <= 50.0) * (0.0017793594306) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 75.0) * (0.00364298724954) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 75.0 && pt <= 100.0) * (0.00363636363636) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 125.0) * (0.00645161290323) +
         (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 2.78571428571 && abs(eta) <= 3.0) * (pt > 150.0) * (0.0) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 20.0 && pt <= 35.0) * (0.0) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 35.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 50.0 && pt <= 75.0) * (0.0) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 75.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 150.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 10.0 && pt <= 15.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 15.0 && pt <= 20.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 20.0 && pt <= 35.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 35.0 && pt <= 50.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 50.0 && pt <= 75.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 75.0 && pt <= 100.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 100.0 && pt <= 125.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 150.0) * (0.0)  
  }
    {13} {

         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 20.0 && pt <= 35.0) * (0.000286232229749) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 35.0 && pt <= 50.0) * (0.000169376693767) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 50.0 && pt <= 75.0) * (0.000480230510645) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 75.0 && pt <= 100.0) * (0.000452386337933) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 100.0 && pt <= 125.0) * (0.000304136253041) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 150.0) * (0.000123026450687) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 35.0) * (0.000205528722639) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 35.0 && pt <= 50.0) * (0.000185322461082) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 75.0) * (0.000535714285714) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 75.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 150.0) * (0.000278396436526) +
         (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 20.0 && pt <= 35.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 35.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 50.0 && pt <= 75.0) * (0.000553403431101) +
         (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 75.0 && pt <= 100.0) * (0.000842105263158) +
         (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 100.0 && pt <= 125.0) * (0.000578703703704) +
         (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 125.0 && pt <= 150.0) * (0.000796178343949) +
          (abs(eta) > 1.5 && abs(eta) <= 2.15) * (pt > 150.0) * (0.000596896140072) +
         (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 20.0 && pt <= 35.0) * (0.000208811860514) +
         (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 35.0 && pt <= 50.0) * (0.000399840063974) +
         (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 50.0 && pt <= 75.0) * (0.000395882818686) +
         (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 75.0 && pt <= 100.0) * (0.000655737704918) +
         (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 125.0 && pt <= 150.0) * (0.00157480314961) +
          (abs(eta) > 2.15 && abs(eta) <= 2.8) * (pt > 150.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 10.0 && pt <= 15.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 15.0 && pt <= 20.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 20.0 && pt <= 35.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 35.0 && pt <= 50.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 50.0 && pt <= 75.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 75.0 && pt <= 100.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 100.0 && pt <= 125.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 2.8) * (pt > 150.0) * (0.0)  
  }
    {22} {

         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 20.0 && pt <= 35.0) * (0.00457971567599) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 35.0 && pt <= 50.0) * (0.00711382113821) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 50.0 && pt <= 75.0) * (0.00304145990075) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 75.0 && pt <= 100.0) * (0.0013571590138) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 100.0 && pt <= 125.0) * (0.000608272506083) +
         (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 125.0 && pt <= 150.0) * (0.00302343159486) +
          (abs(eta) > 0.0 && abs(eta) <= 0.75) * (pt > 150.0) * (0.00205044084478) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 35.0) * (0.014489774946) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 35.0 && pt <= 50.0) * (0.0174203113417) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 75.0) * (0.0166071428571) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 75.0 && pt <= 100.0) * (0.00801282051282) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 100.0 && pt <= 125.0) * (0.00517063081696) +
         (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 125.0 && pt <= 150.0) * (0.00437445319335) +
          (abs(eta) > 0.75 && abs(eta) <= 1.5) * (pt > 150.0) * (0.00406458797327) +
         (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 20.0 && pt <= 35.0) * (0.00497795477173) +
         (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 35.0 && pt <= 50.0) * (0.00205814252637) +
         (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 50.0 && pt <= 75.0) * (0.00418101328087) +
         (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 75.0 && pt <= 100.0) * (0.00412603150788) +
         (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 125.0 && pt <= 150.0) * (0.00356633380884) +
          (abs(eta) > 1.5 && abs(eta) <= 2.25) * (pt > 150.0) * (0.000183183733284) +
         (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 35.0) * (0.0) +
         (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 35.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 75.0) * (0.0) +
         (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 75.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 2.25 && abs(eta) <= 3.0) * (pt > 150.0) * (0.0) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 10.0 && pt <= 15.0) * (0.0) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 15.0 && pt <= 20.0) * (0.0) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 20.0 && pt <= 35.0) * (0.0) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 35.0 && pt <= 50.0) * (0.0) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 50.0 && pt <= 75.0) * (0.0) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 75.0 && pt <= 100.0) * (0.0) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 100.0 && pt <= 125.0) * (0.0) +
         (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 150.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 10.0 && pt <= 15.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 15.0 && pt <= 20.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 20.0 && pt <= 35.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 35.0 && pt <= 50.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 50.0 && pt <= 75.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 75.0 && pt <= 100.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 100.0 && pt <= 125.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 125.0 && pt <= 150.0) * (0.0) +
          (abs(eta) > 4.0) * (pt > 150.0) * (0.0)  
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
