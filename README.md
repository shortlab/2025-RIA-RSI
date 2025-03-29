# 2025-RIA-RSI
Data repository for Aurelien's RSI paper on our fiber RIA spectroscopy facility

Data supporting the experiments can be found in the "Experiments" folder. It includes the code used to run the instruments (Arduino TTL trigger box, Fiber temperature monitoring code), which complements information given in the article to reproduce this experimental setup. This folder also contains data gathered from the literature that can be used to perform spectrum decomposition over the different known defect absorption bands. The "Experimental results for article" folder contains compressed folder with the experimental data (transmitted spectra and temperature data), stamps to identify key moments of each experiments in the spectrometer data, the analysis codes in jupyter notebooks, and python modules that can be used for gaussian decomposition, spectrum smoothing and noise region removal.

The "Simulations" folder contains the OpenMC calculation (jupyter notebook and NIST XCOM data) used to evaluate the actual dose rate by the fiber samples, compared to the nominal dose rate given for an empty Gammacell chamber, and the COMSOL model used to evaluate the maximum difference between the temperature measured by the sensor and the fiber sample temperature.

The supporting data repository can be found at https://doi.org/10.5281/zenodo.15107009. It contains some of the data files that are too heavy for a GitHub repository.
