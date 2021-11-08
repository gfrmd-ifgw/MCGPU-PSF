# MC-GPU PSF

---


This repository contains the tools described in the technical note:
**"Technical Note: MC-GPU breast dosimetry validations with other Monte Carlo codes and Phase Space File implementation"**
RT Massera, RM Thomson and A Tomal.

Each folder contains specific codes and implementations (with further informations of the original sources).


Contens:
1. [MCGPU_files](MCGPU_files/README.md): the modified MC-GPU files to include phase space files
2. [example](example/README.md): a simulation for testing the phase space files
3. [PSFBin2IAEA](PSFBin2IAEA/README.md):  Fortran routine to convert raw binary files to IAEA format

Instructions:
- Compile MC-GPU binary [1](MCGPU_files/) and place on example folder [2](example/). 
- Run the example following the instructions within the folder.
- After the simulation is finished, three files will be generated in examples/images: mcgpu_image_x_psf.raw (x=1,2,3).
- Compile the psfconvert binary [3](PSFBin2IAEA/).
- Run psfbin2fort.py script following the instructions to convert the PSF in binary to Fortran compatible.
- Copy the psfconvert executable to the directory containing the PSF binary (the one generated with the Python script) and run
- If successful, two files will be generated: filename.IAEAheader and  filename.IAEAphsp.
- Those PSF now can be loaded in other MC codes.

Each code has a specific license, indicated in each folder. Otherwise, if not indicated, the contributions of this work are released under the LICENSE file in the main directory.
