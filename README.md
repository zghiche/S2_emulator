# Stage-2 Emulator wrapper

This repository is a first step in integrating the standalone C++ emulator in a python code.
It is possible to use a configuration file (`config.yaml`) and easily change the input configuration parameters of the Stage-2 emulator.

To create the binding between Python and C++ code, a special library has been used: [cppyy](https://cppyy.readthedocs.io/en/latest/). Some short example about how import C++ libraries, how to call class into Python, etc are reported in [this](https://github.com/mchiusi/python-bindings/tree/main) repository on my profile.

For the moment, there is no `Makefile` to compile the code. The command that should be used is the following. The shared library has to be inserted in the python file `run_configuration.py`.
```
g++ -std=c++14 -shared -o lib_configuration.so L1Trigger/L1THGCal/src/backend_emulator/*.cc -fPIC -I.
```
