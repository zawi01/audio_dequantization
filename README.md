
Audio Dequantization Using (Co)Sparse (Non)Convex Methods
========================================================================
[Pavel Záviška](https://orcid.org/0000-0003-2221-2058), [Pavel Rajmic](https://orcid.org/0000-0002-8381-4442), and [Ondřej Mokrý](https://orcid.org/0000-0003-1806-5809)
------------------------------------------------------------------------

This readme file describes the MATLAB toolbox accompanying the article from the title.  
Published version of the article is available at https://ieeexplore.ieee.org/document/9414637.  
Postprint is also available at https://arxiv.org/abs/2010.16386.


### Requirements
The code has been developed in MATLAB version R2019b and it relies on the LTFAT toolbox (version 2.4.0 was used).
To reproduce all the results from the paper, evaluate the dequantized signals using the PEMO-Q ([paper](https://ieeexplore.ieee.org/document/1709880), [website](https://www.hoertech.de/de/produkte/pemo-q.html)) metric. 
Nevertheless, for running the dequantization tests, this software is not essential.

### Quick Tutorial
To use this declipping toolbox, download all the files, add them to the MATLAB path and make sure that the LTFAT toolbox is properly installed.

The toolbox is organized as follows:
   - "Algorithms" folder contains implementations of all ten algorithms used in the experiments. 
   - "Sounds" folder contains wav-files used for testing.
   - "Tools" folder contains support functions for the dequantization algorithms, quantizing the signal, etc.

The root folder contains two main files. 

The m-file "dequantization_main.m" is designed to run one dequantization experiment with selected settings and parameters.
It is possible to select the testing audio file (`audio_file`), level of quantization in bits per sample (`param.wordlength`), and dequantization algorithm (`param.algorithm`).
Other options, such as frame settings and general options of the algorithms are also adjustable. 
Except for the SPADQ algorithms, the algorithm parameters are set directly in the respective m-files.
Note that default values are the values used for the experiments in the paper.

The other main file "dequantization_whole_database.m" serves to easily reproduce the results from the paper.
One can select algorithms (`alg_idxs`), sounds (`sound_idxs`), and word lengths (`wordlengths`) and run the experiments at once. 
It is also possible to enable or disable storing the dequantized signals (`STORE_DEQ_SOUNDS`) or computing and storing values of the objective function (`STORE_OBJ_PROCESS`) and SDR (`STORE_dSDR_PROCESS`) in each iteration. 

To reproduce the &Delta;SDR values from paper, enable the `paramsolver.dec_term` parameter, which terminates the algorithms after a SDR drop.
On the other hand, to reproduce the PEMO-Q ODG values, let this parameter disabled and run full 500 iterations. 

### How to cite this toolbox
Please cite the following paper:

P. Záviška, P. Rajmic and O. Mokrý, 
"Audio Dequantization Using (Co)Sparse (Non)Convex Methods," 
in *2021 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP)*, 
2021, pp. 701-705, doi: 10.1109/ICASSP39728.2021.9414637.

### License
The code of this toolbox is distributed under the terms of the GNU Public License version 3 (http://www.gnu.org/licenses/gpl.txt).

--------------------------------------------------
Pavel Záviška, Brno University of Technology, 2021