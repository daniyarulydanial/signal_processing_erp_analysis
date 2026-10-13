# EEG Signal Processing and ERP Analysis

MATLAB research script for preprocessing EEG recordings and generating event-related potential (ERP) analyses using EEGLAB and associated ERP functionality.

## Objective
Explore a configurable EEG-to-ERP workflow for N-back recordings from the COG-BCI dataset, including artifact handling, epoching, and ERP visualization.

## Pipeline (as configured in the source)
- Load EEGLAB `.set` recordings and event trigger mappings.
- Downsample recordings to 250 Hz; apply a 1 Hz high-pass filter and 50 Hz line-noise cleaning.
- Perform component/artifact screening with an ICLabel probability threshold of 0.90.
- Epoch signals from -200 to 800 ms and baseline-correct using -200 to 0 ms.
- Apply amplitude and moving-window artifact screening.
- Save processed data, ERP outputs, plots, and processing logs.

These are script configuration defaults, **not a statement that every processing stage was independently validated on every subject**.

## Files
- `final_erp_cogbci.m`: main MATLAB script.

## Requirements and setup
MATLAB, EEGLAB, and the plugins/functions called by the script (check the installed EEGLAB / ERPLAB versions). Obtain the COG-BCI dataset according to its access rules. Before running, edit the `USER SETTINGS` section, including `data_root`, selected subjects and sessions, `task_filename`, and `triggerlist_path`. The checked-in configuration uses `twoBACK.set` and an illustrative subset of subjects; do not assume that it processes the full dataset without adjustments.

## Status and limitations
This is a research script/archive, not a packaged or fully independently reproducible pipeline. It requires external EEG recordings and a correctly formatted trigger list. Future improvements include dependency versions, a dataset layout example, execution instructions, and representative ERP figures.

## Author
Danial Daniyaruly · [LinkedIn](https://www.linkedin.com/in/danial-daniyaruly/)
