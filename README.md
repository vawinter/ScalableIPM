# A Scalable Integrated Population Model for Estimating Abundance in Gamebird Management

## Integrated Population Model (IPM) for Wild Turkey in Pennsylvania

### Overview

This repository contains two Integrated Population Models (IPMs) for estimating the abundance, survival, and recruitment of wild turkeys in Pennsylvania. Both models account for age class and sex differences:

The **Research IPM** operates at the scale of three Wildlife Management Units (WMUs).

The **Operational IPM** leverages informative priors from the Complex IPM posteriors to estimate parameters at the broader **wildlife management region** scale.

### Key Demographic Parameters:

**Survival**

**Abundance**

**Recruitment**

## Research IPM (R-IPM)

The **Research IPM** estimates demographic parameters at the WMU scale for males and females.

### Estimations by:

**Season**: Females (November), Males (May)

**Geographic Scale**: 3 WMUs

**Sex**: Male, Female

**Age Class**: Adult, Juvenile

### Models

**Males**:

**Dead Recovery Model**: Estimates harvest and survival rates at the WMU or regional level.

**Lincoln-Peterson Estimator***: Estimates abundance at four biologically relevant time points throughout the year.

**Females**:

**Known-Fate Model**: Estimates annual survival from telemetered females.

**Recruitment Models:**

**Hen with Brood (HWB) Model**: Estimates the number of hens with poults on September 1.

**Poults:Brood (PPB) Model**: Estimates the ratio of poults to broods.

**Lincoln-Peterson Estimator**: Estimates abundance at four biologically relevant time points throughout the year.

### Data File:

`Research_IPM_run.Rdata` – Contains setup data and IPM output.

## Operational IPM (O-IPM)

The **Operational IPM** estimates demographic parameters at the **regional scale** for males and females, using informative priors derived from the *Complex IPM*.

### Estimations by:

**Season**: Females (November), Males (May)

**Geographic Scale**: 9 regions

**Sex**: Male, Female

**Age Class**: Adult, Juvenile

### Models

**Males**:

**Dead Recovery Model**: Estimates harvest and survival rates at the WMU or regional level.

**Lincoln-Peterson Estimator**: Estimates abundance at four biologically relevant time points throughout the year.

**Females**:

**Recruitment Models:**

**Hen with Brood (HWB) Model**: Estimates the number of hens with poults on September 1.

**Poults:Brood (PPB) Model**: Estimates the ratio of poults to broods.

**Lincoln-Peterson Estimator**: Estimates abundance at four biologically relevant time points throughout the year.

### Data File:

`Operational_IPM_run.Rdata` – Contains setup data and IPM output.

## Evaluating O-IPM with Vague Priors (V-IPM)

To assess the impact of informative priors, we also fitted the Operational IPM with vague `Beta(1,1)` priors on female harvest rates and survival.

### Data File:

`vague_IPM.R` – Contains the V-IPM model and the setup data data is `Operational_IPM_run.Rdata`. The set up data does not differ, only priors in the model.