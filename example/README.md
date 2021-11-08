# Example

This is an example of simulation to generate the PSF files for a mammography projection. Notice that the data is stored in the "image/" folder as a raw binary data.
A script is used to convert the file to a PSF following the IAEA format.
The cross sections are adapted from PENELOPE 2018.
To run the simulation, compile the code and paste the executable in this directory. With the command line type:
./MC-GPU_v1.5b.x simulation.in | tee simulation.out
