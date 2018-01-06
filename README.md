SyncDFT
======
A **SyncDFT**, or 'Synchronized DFT' is a DFT where number of samples used for the DFT is as close as possible to number of samples of one period of the base frequency of the signal. (or a multiple of it)

These Matlab programs show how preprocessing of the signal based on Zero Crossing Detection (ZCD) can be used to find the base frequency and achieve a more meaningful und usable frequency spectrum in some cases.

A SyncDFT of a Saw wave signal at 500 Hz looks ![like this](https://www.dropbox.com/s/kssf4342pz6629n/SyncDFT.PNG?dl=0 "SyncDFT"), 
while a FFT with 1024 points looks ![like this](https://www.dropbox.com/s/tjxf5oy151711fr/UnSyncDFT.PNG?dl=0 "FFT").


## Download
* [Version 0.1](https://github.com/f-heil/SyncDFT/archive/master.zip)
* ```$ git clone https://github.com/f-heil/SyncDFT.git
...```

## Usage
The programs here are meant to be used as examples to test whether a ZCD-based SyncDFT will be useful, or even usable at all.

### Third party libraries
If you want graphs to be generated in TikZ:
![matlab2tikz](https://github.com/matlab2tikz/matlab2tikz "matlab2tikz").

This option is commented out by default in the programs to minimize confusion.

## How-to use this code
Points 2.x are only necessary if matlab2tikz is used:

1   : Download the code
  2.1 : Make sure you have met the dependencies for matlab2tikz
  2.2 : Put matlab2tikz in the program folder, or use 'addpath' to tell Matlab where it is.
3   : Run the program in Matlab
4   : All files will be generated in the program folder, and automatically created sub-folders.

## License 
* MIT, see [LICENSE](https://github.com/f-heil/SyncDFT/blob/master/LICENSE) file

## Version 
* Version 0.1

## Contact
#### Developer
* e-mail: heil_florian@gmx.de 
