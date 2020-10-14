# DimuonSpectrum2011 MC&Data

- zh_CN  [简体中文](https://github.com/PhyM73/DimuonSpectrum2011MC/blob/datamaster/README.zh_CN.md)

This is a simple analysis example to compute the spectrum of two muon final state with CMS Open Data. The goal is to reproduce the results histogram in this paper, [Searching in CMS Open Data for Dimuon Resonances with Substantial Transverse Momentum](https://arxiv.org/abs/1902.04222). This repository also made reference to another GitHub project, [DimounSpectrum2011](https://github.com/cms-opendata-analyses/DimuonSpectrum2011)



## Usage

You can run this code in [CMS Open Data VM](http://opendata.web.cern.ch/VM/CMS/2010). If you have not installed the CMSSW area do the following:

```bash
cmsrel CMSSW_5_3_32
```

You can also run this code in a [CMS Open Data Docker](http://opendata.cern.ch/docs/cms-guide-docker) container if you are interested in that:

```bash
docker run --name dimu2011mc -it cmsopendata/cmssw_5_3_32 bash
```

If you already have the CMSSW on your VM or run the docker container, start directly with:

```bash
cd CMSSW_5_3_32/src
cmsenv
```

For this example, you need to create an additional directory, you can call it `WorkDir` or choose another name. The name doesn't matter.
Go to this directory, and download the example code. (the default branch is `datamaster`)

```bash
mkdir WorkDir
cd WorkDir
git clone git://github.com/PhyM73/DimuonSpectrum2011MC.git
```

Go to the example directory, and compile with `scram b`. 

```bash
cd DimuonSpectrum2011MC
scram b
```

The input files are defined in the configuration file `demoanalyzer_cfg.py` and are already in the `datasets` directory. In order not to overwrite the existing DoubleMu2011.root, to which you can compare your output, rename it before you run.

Now run the configuration file. 

```bash
cmsRun demoanalyzer_cfg.py
```

The output of the example is a root file containing several histograms, by default DoubleMu2011.root with 10000 input events (small subset of data). These can be looked at using a Root Browser.

For more detailed information, read the comments in src/DimuonSpectrum2011MC.cc.



### Different branches

There are some differences between the analysis workflow for collision data and MC samples. The MC samples don't have any trigger path but do need to perform some analysis to compute efficiency. Hence the example has different branches for different datasets and analysis platforms, i.e. local computer or condor. Below are the descriptions of different branches. 

- `datamaster` is for code testing and simple analysis, This branch process 10000 events from 2011a datasets by default and this is the default branch.
- `datacondor` is for condor submitting. This branch is for processing all the events found in the path `./data` which should be created in the input executable for condor submitting. You can use this branch to analyze all the events from 2011a datasets by setting a proper number of jobs.
- `mcmaster` is for code testing and simple analysis, This branch process 10000 events from Monte Carlo datasets by default.
- `mccondor` is for condor submitting. This branch process all the events found in the path `./data` which should be created in the input executable for condor submitting. You can use this branch to analyze all the events from MC datasets by setting a proper number of jobs.




## Files description

`datasets` include files below

- `datasets/double`  [DoubleMu primary dataset in AOD format from RunA of 2011](http://opendata.cern.ch/record/17)

- `datasets/single`  [SingleMu primary dataset in AOD format from RunA of 2011](http://opendata.cern.ch/record/32)
- `datasets/mc` 
   - [Simulated dataset DYJetsToLL_M-10To50_TuneZ2_7TeV-pythia6 in AODSIM format for 2011 collision data](http://opendata.cern.ch/record/1393) 
   - [Simulated dataset DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola in AODSIM format for 2011 collision data](http://opendata.cern.ch/record/1395)

- `datasets/double100` Files in `datasets/double` divided into every 100 indexes.
- `datasets/mc100` Files in `datasets/mc` divided into every 100 indexes.

- [CMS list of validated runs for primary datasets of 2011 data taking](http://opendata.cern.ch/record/1001)

You can obtain the luminosity information from [CMS luminosity information for 2011 CMS open data](http://opendata.cern.ch/record/1051) as needed.



