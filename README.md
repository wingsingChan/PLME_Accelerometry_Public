# Code and Data for "User- and budget-friendly accelerometers highlight vulnerability to poaching and weather for the Critically Endangered big-headed turtle (*Platysternon megacephalum*)"

This repository contains example data collected from an user- and budget-friendly accelerometer that was installed on the Critically Enangered big-headed turtle (*Platysternon megacephalum*). The accelerometer was developed by our team and the code and resources for building it are open-source and are available at [Hackster.io](https://www.hackster.io/brian-k2/accelerometer-and-data-logger-for-small-animal-research-c877c6). 

In addition to the original accelerometry datasets, all codes generated for the analysis of the collected datasets are included in this repository. 

> Chan, W. S., Katona, B., Sung, Y. H., Bonebrake, T.C. (*In prep*). User- and budget-friendly accelerometers highlight vulnerability to poaching and weather for the Critically Endangered big-headed turtle (*Platysternon megacephalum*).

## Data

All data are stored in the `data/` directory. 

### Metadata

Metadata for each round of accelerometer deployment, including biometric measurements of the monitored individual turtles, is recorded in `PLME_accelerometry_meta.csv`.

A description of the variables in the metadata is provided in the table below:

| Variables                 | Description                                                                                                                      |
|---------------------------|----------------------------------------------------------------------------------------------------------------------------------|
| Species                   | Abbreviation of the species name, *[**Pl**]{.underline}atysternon [**me**]{.underline}gacephalum* (i.e., PLME)                   |
| new.recapture             | Indicates whether the individual is newly captured or recaptured                                                                 |
| turtle.ID                 | Individual identity identified through marginal scale notching                                                                   |
| pit.tag                   | Pit tag ID                                                                                                                       |
| date                      | Date when the turtle was captured                                                                                                |
| sex                       | Sex of the individual turtle, as either male (`M`) or female (`F`)                                                               |
| cl                        | Straight-line carapace length measured with a caliper                                                                            |
| pl                        | Straight-line plastron length measured with a caliper                                                                            |
| wt                        | Body weight measured using an electronic balance                                                                                 |
| hw                        | Head width measured across the head of the animal using a caliper                                                                |
| tl                        | Tail length measured from cloaca to tail tip using a caliper                                                                     |
| leech                     | Number of leeches present on the individual's body                                                                               |
| acce.roundNo              | Accelerometer deployment round number                                                                                            |
| acce.startDate            | Date when the accelerometer was activated                                                                                        |
| acce.startTime            | Datetime when the accelerometer was activated                                                                                    |
| acce.delayTime.d          | Days programmed until the accelerometer wakes from sleep mode                                                                    |
| acce.samplRate.Hz         | Sampling rate of the accelerometer in Hertz (Hz)                                                                                 |
| acce.samplInterval.s      | Sampling duration for each recording interval in seconds                                                                         |
| acce.sleepInterval.m      | Time interval between samples in minutes                                                                                         |
| acce.batteryCap.mAh       | Capacity of the battery in milliampere-hours (mAh) used to power the accelerometer unit                                          |
| acce.estLifespan.d        | Estimated lifespan of the accelerometer in days                                                                                  |
| acce.estEndTime           | Estimated end time for the accelerometer's operation                                                                             |
| acce.netWeight            | Net weight of the accelerometer system in grams                                                                                  |
| transmitterFreq           | Frequency of the radiotransmitter deployed with the accelerometer (if applicable)                                                |
| filename                  | Filename for which the respective accelerometer data is stored                                                                   |

### Raw Accelerometry Data

Raw accelerometer data can be found in `data/accelerometer/`, organized by the round of deployment in a format that specifies both the start date and round number of the deployment as `startDate-roundNo/`. Each subdirectory contains spreadsheets with raw data from individual turtles, named by the abbreviation of the species (`PLME`), individual ID, sex (denoted by `M` for Male and `F` for Female), and the start date when the accelerometer was activated (i.e., `PLME_IDSex_acceleration_logfile_data_startDate.csv`).

Each dataset consists of five columns: the first column represents the time from the built-in Real-Time Clock, followed by acceleration along the x-, y- and z-axes in that order. The last column includes the temperature recorded by temperature sensor onboard. 

### Raw Environmental Data

#### Temperature

Ambient air and water temperatures were measured in the study stream respectively at two localities using iButtons. The raw temperature data are stored in `data/thermochrone/` and are organized into deployment round as `data/thermochrone/startDate-roundNo/`. Each subdirectory contains spreadsheets with raw temperature data recorded by individual iButtons. These files are sorted and named according to the iButton serial number, site name, locality of placement, medium (air or water), and the start date (i.e., `serial_site_loc_medium_startDate.csv`).

Each dataset consists of 4 columns: the date (`Date`) and time (`Time`) when the temperature data was logged, the temperature unit (`Unit`), and the measured temperature value in the specified unit (`Value`). 

#### Rainfall

Hourly precipitation data collected from the weather station at the project site is provided by the Hong Kong Observatory and is stored in `rainfall_kfbg_hourly.csv`. 

### Processed Data

All raw accelerometry data were processed and aggregated into `PLME_accelerometry_aggregated.csv` for fitting individual-specific hidden Markov models (HMMs). Outputs from the most parsimonious HMMs were further summarised into `PLME_accelerometry_hmm.csv` for further investigation using Generalized Additive Models (GAM). 

Description of the variables in the processed acclerometery datasets are listed in the table below: 

| Variables      | Description                                                                                                                          |
|----------------|--------------------------------------------------------------------------------------------------------------------------------------|
| ID             | Track ID specific to individual identity, accelerometer round number, and track number, formatted as `turtleID-acce.roundNo-trackNo` |
| ID_old         | Individual and deployment round-specific identity, represented by `turtleID-acce.roundNo`                                            |
| acce.roundNo   | Accelerometer deployment round number                                                                                                |
| turtle.ID      | Individual identity identified through marginal scale notching                                                                       |
| sex            | Sex of the individual turtle, as either male (`M`) or female (`F`)                                                                   |
| cl             | Straight-line carapace length measured with a caliper                                                                                |
| time           | Datetime when the aggregated acceleration data is recorded                                                                           |
| date           | Date when the aggregated acceleration data is recorded                                                                               |
| hour           | Hour corresponding to the recorded time                                                                                              |
| month          | Month corresponding to the recorded date                                                                                             |
| year           | Year corresponding to the recorded date                                                                                              |
| julian         | Julian date corresponding to the recorded date                                                                                       |
| x              | Dynamic body acceleration in the x-axis                                                                                              |
| y              | Dynamic body acceleration in the y-axis                                                                                              |
| z              | Dynamic body acceleration in the z-axis                                                                                              |
| ODBA           | Overall Dynamic Body Acceleration                                                                                                    |
| ODBA.scaled    | Log-transformed Overall Dynamic Body Acceleration                                                                                    |
| airTemp        | Ambient air temperature at the time of data recording                                                                                |
| waterTemp      | Water temperature at the time of data recording                                                                                      |
| rainfall       | Amount of rainfall recorded during the observation period                                                                            |
| roll_sum_24h   | Cumulative rainfall over the past 24 hours                                                                                           |
| heavyRain      | Boolean expression whether the cumulative rainfall in the past 24 hours (`roll_sum_24h`) exceeds 30 mm                               |
| acce.netWeight | Net weight of the accelerometer system in grams                                                                                      |
| state          | Motion state as either stationary (`1`) or mobile (`2`) identified by the hidden Markov model                                        |

## Scripts

All R codes used to analyze the data are available in the directory `script/`.

-   `PLME_AccelerometryAnalysis_01_Preprocessing.R` contains codes to process the original time series data and aggregate it into formats suitable for the subsequent analyses.

-   `PLME_AccelerometryAnalysis_02_HMM.R` contains codes to fit individual-specific hidden Markov models to derive the motion states of each respective animal at each time step.  

-   `PLME_AccelerometryAnalysis_03_GAM.R` contains code to fit Generalized Additive Model (GAM) to investigate the time activity patterns of *P. megacephalum*

## Contact Information

Questions regarding the code or the data should be directed to the project author Wing Sing Chan at [wschan1021\@gmail.com](mailto:wschan1021@gmail.com).

## License

This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative Commons Attribution - Non-Commercial 4.0 International License</a>.<br />
<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commos License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png" /></a>
