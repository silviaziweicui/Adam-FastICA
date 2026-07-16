# Adam-FastICA

> **Notice**
> This repository is released to facilitate the peer-review process. It contains the implementation used in the submitted manuscript. The documentation and examples may be further refined after the peer-review process.

## Overview

Adam-FastICA is a high-density surface electromyography (HD-sEMG) decomposition algorithm that integrates the FastICA framework with the Adam optimization strategy to improve the convergence speed and stability of motor unit (MU) decomposition.

The repository provides the reference MATLAB implementation used in the submitted manuscript, including both offline decomposition and online adaptive decomposition.

## Requirements

- MATLAB R2018b or later
- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox

## Input Data Format

The input HD-sEMG signal should be organized as

```
channels × samples
```

For example,

```matlab
data = randn(64, 120000);
```

where

- rows correspond to recording channels;
- columns correspond to time samples.

## Quick Start

### 1. Prepare the data

Load your HD-sEMG recording into MATLAB.

```matlab
load('your_data.mat');
```

The training and testing datasets should be stored as matrices with dimensions

```
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

- Additional examples and documentation may be added in future updates.
- The implementation is intended for academic research and reproducibility.


## Related Manuscript

This repository accompanies the implementation described in the following manuscript, which is currently under peer review.

> Title **Real-time Motor Unit Tracking: Fast and Adaptive Decomposition of Non-Stationary Surface Electromyographic Signals**

The manuscript information will be updated after publication.


## Citation

If you find this repository useful in your research, please consider citing the corresponding publication after it becomes available.

BibTeX information will be added after the manuscript is accepted.


## References

1. Hyvärinen, A., Karhunen, J., & Oja, E. (2001). *Independent Component Analysis*. John Wiley & Sons.

2. Kingma, D. P., & Ba, J. (2015). Adam: A Method for Stochastic Optimization. *International Conference on Learning Representations (ICLR).*

3. Negro, F., Muceli, S., Castronovo, A. M., Holobar, A., & Farina, D. (2016). Multi-channel intramuscular and surface EMG decomposition by convolution kernel compensation. *IEEE Transactions on Biomedical Engineering*, 63(7), 1493–1504.


## License

This repository is provided for academic and non-commercial research purposes.


## Contact

For questions, bug reports, or suggestions, please open an Issue in this repository or contact the repository owner.
