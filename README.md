# Adam-FastICA

Adam-FastICA is a high-density surface electromyography (HD-sEMG) decomposition algorithm that integrates the FastICA framework with the Adam optimization strategy to improve the convergence efficiency and stability of motor unit (MU) decomposition.

This repository provides the reference MATLAB implementation used in the corresponding research work, including both offline decomposition and online adaptive decomposition.

> **Note**
>
> This repository was originally prepared for the peer-review process and has been made publicly available following acceptance of the corresponding manuscript. The code and documentation may be further updated to improve clarity, usability, and reproducibility.

## Requirements

* MATLAB R2018b or later
* Signal Processing Toolbox
* Statistics and Machine Learning Toolbox

## Input Data Format

The input HD-sEMG signal should be organized as

```text
channels × samples
```

For example:

```matlab
data = randn(64, 120000);
```

where:

* rows correspond to recording channels;
* columns correspond to time samples.

## Quick Start

### 1. Prepare the data

Load your HD-sEMG recording into MATLAB.

```matlab
load('your_data.mat');
```

The training and testing datasets should be stored as matrices with dimensions:

```text
channels × samples
```

### 2. Set parameters

```matlab
FS = 2048;
win_size_ms = 200;
niter = 200;
```

### 3. Offline decomposition

```matlab
[signal, parameters, signalprocess] = adam_fastica(niter, train_data, FS);
```

### 4. Online decomposition

```matlab
[ipt, spikes] = online_decomp(test_data, win_size_ms, signalprocess, FS);
```

## Example

```matlab
%% Load data
load('your_data.mat');

%% Parameters
FS = 2048;
win_size_ms = 200;
niter = 200;

%% Offline decomposition
[signal, parameters, signalprocess] = adam_fastica( ...
    niter, train_data, FS);

%% Online decomposition
[ipt, spikes] = online_decomp( ...
    test_data, win_size_ms, signalprocess, FS);
```

## Notes

* The implementation is provided for academic research and reproducibility.
* The offline stage is used to estimate the decomposition parameters from the training data.
* The online stage uses the estimated parameters for adaptive decomposition of subsequent HD-sEMG signals.
* Additional examples and documentation may be added in future updates.

## Related Publication

This repository accompanies the following publication:

> **Real-time Motor Unit Tracking: Fast and Adaptive Decomposition of Non-Stationary Surface Electromyographic Signals**

The manuscript has been **accepted for publication**. Full publication information, including the journal, volume, issue, page numbers, and DOI, will be added to this section once available.

## Citation

If you use Adam-FastICA or the code provided in this repository in your research, please cite the corresponding publication:

```text
Cui, Z., et al.
"Real-time Motor Unit Tracking: Fast and Adaptive Decomposition of
Non-Stationary Surface Electromyographic Signals."
[Forthcoming publication]
```

The complete citation and BibTeX entry will be updated once the final publication information is available.

## References

1. Hyvärinen, A., Karhunen, J., & Oja, E. (2001). *Independent Component Analysis*. John Wiley & Sons.

2. Kingma, D. P., & Ba, J. (2015). Adam: A Method for Stochastic Optimization. *International Conference on Learning Representations (ICLR).*

3. Negro, F., Muceli, S., Castronovo, A. M., Holobar, A., & Farina, D. (2016). Multi-channel intramuscular and surface EMG decomposition by convolution kernel compensation. *IEEE Transactions on Biomedical Engineering*, 63(7), 1493–1504.

## License

This repository is provided for academic and non-commercial research purposes.

## Contact

For questions, bug reports, or suggestions, please open an Issue in this repository or contact the repository owner.
