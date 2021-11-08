# PSFIAEAGen
This repository uses an adapted penEasy (v. 2019) fortran routine to write Phase Space Files in IAEA format.
We use the IAEA phase space file routines (https://www-nds.iaea.org/phsp/phsp.htmlx), see LICENSE and IAEA_PHSP_README.TXT for more information.

The file PhaseSpaceFileConverter.F contains the Fortran routines for the PSF generation. The original disclaimer: \
!********************************************************************************************* \
!*penEasy copyright and disclaimer \
!*penEasy \
!*Copyright (c) 2003-2019 \
!*Universitat Politecnica de Catalunya \
!*Permission to use, copy, modify and re-distribute copies of this software, or parts of it,\
!*and its documentation for any purpose is hereby granted without fee, provided that the \
!*above copyright notice appears in all copies and that both that copyright notice and this \
!*permission notice appear in all supporting documentation. The Universitat Politecnica \
!*de Catalunya makes no representations about the suitability of this software for any \
!*purpose. It is provided \as is" without express or implied warranty. \
!**********************************************************************************************

How to compile: run the files "compile_script.bat" for windows or "compile_script_linux.sh" for linux. A fortran compiler is necessary.\
We also provide the compiled versions (static builds).\
A python script called "psfbin2fort.py" is also provided and serves as "glue" between the binary and the fortran routine. (Numpy and Scipy packages are required).



