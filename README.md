# Stage-2 Emulator wrapper

This repository is a first step in integrating the standalone C++ emulator in a python code.
It is possible to use a configuration file (`config.yaml`) and easily change the input configuration parameters of the Stage-2 emulator.

To create the binding between Python and C++ code, a special library has been used: [cppyy](https://cppyy.readthedocs.io/en/latest/). Some short example about how import C++ libraries, how to call class into Python, etc are reported in [this](https://github.com/mchiusi/python-bindings/tree/main) repository on my profile.

For the moment, there is no `Makefile` to compile the code. The command that should be used is the following. The shared library has to be inserted in the python file `run_emulator.py`.
```
g++ -std=c++14 -shared -o lib_configuration.so L1Trigger/L1THGCal/src/backend_emulator/*.cc -g -fPIC -I.
```

In the main python file `run_emulator.py` and in `data_handle/event.py` data are prepared and injected. 

Currently I'm studying the back-end mapping (TC to S1, S2 links and frames) to understand how the data has to be prepared starting from a simple simulated dataset. XML files [here](https://gitlab.cern.ch/hgcal-tpg/mapping/-/tree/master) and copied in the `config_files` directory. The geometry is read from the XML file in `data_handle/geometry.py`.


### Installation
To run the Stage-2 emulator using the Python wrapper, some packages need to be installed in your environment. Here, I'll list the most important ones:
I suggest creating a `conda` environment for easy inclusion and activation of the following packages.
```
conda create --name <env_name> python=3.7
conda install -c conda-forge cppyy 

conda install -c conda-forge pyyaml

pip install matplotlib
pip install scipy
pip install awkward
pip install uproot

conda install gxx
```

### Useful commands
```
# run the emulator without producing any plots
python run_emulator.py

# processing 2 events and producing plots post unpacking and post seeding
python run_emulator.py -n 2 --plot --col

# processing 100 events and plotting the seeding efficiency vs eta, pT and energy threshold
python run_emulator.py -n 100 --thr_seed

# to specified the particles (photons or pions) and pileup (0 or 200)
python run_emulator.py -n 100 --pileup PU0 --particles photons --thr_seed
```
Other options can be displayed using `python run_emulator.py --help`.
Single particle plots and efficiency plots can be found [here](https://mchiusi.web.cern.ch/Emulator_seeding/).

Some emulator parameters can be easily changes via the config `yaml` file. For instance, the seeding threshold (in GeV) or the seeding window width.

### Data samples
I'm producing new data samples using Geometry V16 following the recipe outlined [here](https://twiki.cern.ch/twiki/bin/viewauth/CMS/HGCALTriggerPrimitivesSimulation) and considering [Phase2Fall22campaign](https://cmsweb.cern.ch/das/request?view=list&limit=50&instance=prod%2Fglobal&input=dataset+dataset%3D%2F*%2FPhase2Fall22DRMiniAOD-*125X*%2FGEN-SIM-DIGI-RAW-MINIAOD). Public samples (photons and pions) that I'm using are stored in this eos folder:
```
/eos/home-m/mchiusi/s2_emulator
```
The location can be change via the `yaml` configuration file.
These samples have been produced using threshold algorithm in the CE-E and STC for CE-H and CE-S. A simulation of best choice algorithm (BS) for the CE-E is implemented in the Python wrapper.
