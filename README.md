# Actin Tracking Toolbox

Welcome to the Actin Tracking Toolbox, a MATLAB toolbox for tracking actin filaments in time-lapse microscope datasets. This toolbox is being developed in collaboration with the Leinwand Lab.

## Features
- Segmentation of actin fibers
- Accurate measurement of fiber length
- Data analysis tools, e.g. measuring fiber migration speed

## Installation and Usage

### Prerequisites

This toolbox requires the following MATLAB toolboxes to be installed:
* Image Processing Toolbox
* Curve Fitting Toolbox

### Additional toolboxes

The following toolboxes are also required. Please visit the pages listed below and download the latest release:
* [Cell Tracking Toolbox](https://github.com/Biofrontiers-ALMC/cell-tracking-toolbox)
* [Bioformats for MATLAB](https://github.com/Biofrontiers-ALMC/bioformats-matlab)

### Actin Tracking Toolbox

Finally, [download](https://github.com/Biofrontiers-ALMC/actin-tracking-toolbox/releases) the latest release.

## Usage

This toolbox defines two main classes: `ActinData` and `ActinTracker`.

### Processing movies

To process movie, use the `ActinTracker` class:

```matlab
%Create a new ActinTracker object
AT = ActinTracker;

%Run processing with default settings
process(AT, 'input.nd2', 'outputDir')
```

The processing script will return two files per movie processed: an AVI file showing both segmentation and tracking results, and a MAT file containing tracked data. A txt file containing the processing settings will also be output each time the `process` function is run.

### Analyzing data

To analyze the resulting data, use the `ActinData` class:

```matlab
AD = ActinData;

AD = importdata(AD, 'path\to\BetaWT_50_1_002.mat');
AD = analyze(AD);
```


## Bug reports and feature requests

Please use the [Issue Tracker](https://github.com/Biofrontiers-ALMC/cell-tracking-toolbox/issues) to file a bug report or to request new features.

## Development

The source code can be cloned from this repository
```git
git clone git@github.com:Biofrontiers-ALMC/cell-tracking-toolbox.git
```

## Acknowledgments

- Jian Wei Tay
- Lindsey Lee
- Sarah Lehman
- Leslie Leinwand