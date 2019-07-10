DelphesNtuplizer
=============

This package allows you to produce to flat Ntuples from Delphes PhaseII samples.

Table of contents
=================
  * [Clone and initialisation](#clone-and-initilisation)
  * [Produce validation Delphes samples](#producing-delphes)
  * [Produce Delphes Flat trees](#producing-flatrees)


Clone and initialisation
========================

If you do not attempt to contribute to the repository, simply clone it:
```
git clone git@github.com:recotoolsbenchmarks/DelphesNtuplizer.git
```

If you aim at contributing to the repository, you need to fork this repository (via the fork button) and then clone the forked repository:
```
git clone git@github.com:YOURGITUSERNAME/DelphesNtuplizer.git
cd DelphesNtuplizer
git remote add upstream git@github.com:recotoolsbenchmarks/DelphesNtuplizer.git
```
You can then regularly update your fork via:
```
git fetch upstream && git merge upstream/master
```

If you want to submit a new feature to ```recotoolsbenchmarks/DelphesNtuplizer``` you have to do it via pull-request (PR):
So, first commit and push your changes to ```YOURGITUSERNAME/DelphesNtuplizer``` and then make a PR via the github interface. 

This package requires Delphes to be installed, and CMSSW for gcc, ROOT, FWLite and other dependencies:

```
cmsrel CMSSW_10_0_5
cd CMSSW_10_0_5
cmsenv
cd ..
git clone https://github.com/delphes/delphes.git
cd delphes
cp libDelphes.so ..
git checkout tags/3.4.2pre17
./configure
sed -i -e 's/c++0x/c++1y/g' Makefile
make -j 10
```
Make a dummy test that Delphes runs properly on GEN-SIM samples (once few events have been processed you can stop processing with CTRL+C):

```
./DelphesCMSFWLite cards/gen_card.tcl test_gensim.root /eos/cms/store/relval/CMSSW_10_3_0/RelValZMM_13/GEN-SIM-RECO/PU25ns_103X_upgrade2018_realistic_v7-v1/10000/50FEF759-6699-7344-82CE-894E8A724442.root
```

Produce validation Delphes samples 
===================================

For producing new validation samples you need to set-up the proper environment

```
cd CMSSW_10_0_5
cmsenv
cd ../delphes
```

The Delphes card stored in this repository is configured so that all needed information for validation is available in the Delphes output (gen particles, PF candidates etc ...)
As a result, the output size increases substantially compared to normal Delphes samples. 

To produce Delphes validation samples run this command (by changing the appropiate input GEN-SIM file of interest): 

```
./DelphesCMSFWLite ../cards/CMS_PhaseII_200PU_v04VAL.tcl delphes_input.root /eos/cms/store/relval/CMSSW_10_3_0/RelValZMM_13/GEN-SIM-RECO/PU25ns_103X_upgrade2018_realistic_v7-v1/10000/50FEF759-6699-7344-82CE-894E8A724442.root
```

Produce Delphes flat trees
==========================

Set up the proper environment:

```
cd CMSSW_10_0_5
cmsenv
cd ..
```

The following command will produce a flat Ntuple, with 10 events.

``` 
python bin/Ntuplizer.py -i delphes/delphes.root -o flat_tree.root -n 10
```

