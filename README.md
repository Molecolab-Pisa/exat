EXAT - EXcitonic Analysis Tool
--------------------------------------------

  EXAT is a tool developed by the Molecolab group at the
  Department of Chemistry, University of Pisa, Italy.
  For the other Molecolab tools please visit the website
  http://molecolab.dcci.unipi.it/software.html

  EXAT is written in python. It is designed to compute the 
  excitonic properties of a multichromophoric system and 
  to simulate absorption and circular dichroism spectra.

  The program itself is divided into two executables,
  and a few modules; the executables are:

  - `exat.py`:     the program that computes the excitonic properties

  - `spectrum.py`: the program that computes the linear absorption and
                 circular dichroism spectra from the output of exat.py

Disclaimer and copyright:
-------------------------
  EXAT is Copyright (C) 2014-2024 L. Cupellini, S. Jurinovich,
  C.A. Guido, and B. Mennucci

  The terms for using, copying, modifying, and distributing EXAT 
  are specified in the file LICENSE.
      
Contacts:
-------------------------
  Lorenzo Cupellini

  lorenzo.cupellini@unipi.it

  Dipartimento di Chimica e Chimica Industriale

  Via G. Moruzzi 13, I-56124 Pisa PI, Italy

Documentation:
-------------------------
  The main documentation for EXAT usage is in doc/
  There is also a limited help for commands:
```
    exat.py -h
    spectrum.py -h
```

  Some examples of usage can be found in examples/

System Requirements:
-------------------------
  The program has been tested with the following configuration
  on linux systems:

*    `python` version 2.7.2, 2.7.17, 3.7.7
*    `numpy` version 1.9, 1.16, 1.21
*    `scipy` version 1.7.3 (_*_)
*    `matplotlib` version 2.2.3, 3.2.2 (_*_)(_**_) 

  (_*_)  Only needed to run `spectrum.py`
  (_**_) Only needed for live spectrum plotting 

  Please, run tests to check your system compatiblity. 
  Note that the GUI will only work with python3
  
Installation and Testing:
-------------------------
  Once verified the system requirementsi, the executables `exat.py` 
  and `spectrum.py` can be directly run. 

  You may add the `src/` directory to your PATH, or link `exat.py`
  and spectrum.py in your `bin/` directory.

  To run tests go to `tests/` and run
  ```
  $ make
  ```

  Check carefully the `POSSIBLE FAILURE` indications. 
  A sign change in the exciton coefficients is not a problem.

Program's citation:
-------------------------
  Please cite this program as:

  S. Jurinovich, L. Cupellini, C.A. Guido, B. Mennucci,
  "EXAT - EXcitonic Analysis Tool", 
  J. Comput. Chem. 2018, 39, 279-286, DOI: 10.1002/jcc.25118 

Reference papers:
-------------------------
  Please cite the following papers.

  S. Jurinovich, G. Pescitelli, L. Di Bari, B. Mennucci,
  PCCP, 2014, 16, 16407-16418 (doi: 10.1039/c3cp55428g)

  S. Jurinovich, C. A. Guido, T. Bruhn, G. Pescitelli, B. Mennucci,
  ChemComm, 2015, 51, 10498-10501 (doi: 10.1039/C5CC03167B)

