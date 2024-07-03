# Deep depressions in the coastal Red Sea directly evidence oxygen depleted zones in a coral reef domain but sustain resistant ecosystems


## Abstract

Persistent oxygen-depleted zones are known primarily from enclosed basins in temperate regions or the open-ocean (including oxygen minimum and limiting zones). However, little is known about the potential for such zones in tropical coastal domains, even though warmer temperatures and complex geomorphological structures in some domains increase their likelihood. Here, we show two subsurface oxygen-depleted zones within deep depressions of the Red Sea’s Difaht Farasan—a carbonate platform hosting the world’s third-largest contiguous tropical coral reef system3. One zone maintains suboxic oxygen levels (~11-14 µmol kg⁻¹), while the other sustains oxygen levels below detection (<2 µmol kg⁻¹). The suboxic zone shows no fixed nitrogen loss, while the near-anoxic zone hosts anaerobic microbial populations and shows signs of nitrogen loss. We propose that the Red Sea’s warm and saline environment interacts with the semi-enclosed depressions to restrict vertical mixing, limiting oxygen resupply at depth. Yet, unlike most other oxygen-depleted zones, these zones support resistant aerobic respiring taxa, indicating an unusual capacity to reduce aerobic oxygen demands at high temperatures (>22°C). Targeted exploration of deep tropical coastal environments is crucial to determine if similar zones and the resistant taxa they sustain exist beyond the Difaht Farasan and their potential responses to climate change.

## Publication
<insert link to publication>

## Authors
Shannon G. Klein, Larissa Frühe, Anieka J. Parry, Fabio Marchese, Megan K. B. Nolan, Elisa Laiolo, Kah Kheng Lim, Alexandra Steckbauer, Jessica Breavington, Christopher A. Hempel, Kate von Krusenstiern, Froukje M. van der Zwan, Eleonora Re, Taiba Alamoudi, Jacqueline V. Alva Garcia, Silvia Arossa, Carlos Angulo-Preckler, Mattie Rodrigue, Vincent A. Pieribone, Mohammad A. Qurban, Francesca Benzoni, and Carlos M. Duarte.

## Table of Contents
* [Data storage and availability](#data-avail)
* [Sampling Overview](#sampling-overview)
* [Water chemistry data processing](#water-chem-data)
* [Sequence data processing](#rawread-proc)
* [Data analysis](#stats)
* [Contact](#contact)
<!-- * [License](#license) -->

## Data storage and availability

Raw sequencing data can be accessed in the Sequence Read Archive (SRA) of NCBI under PRJNA899789. 
Bacterial sequences can be downloaded using accession numbers SRR22244354 - SRR22244374
Protist sequences can be downloaded using accession numbers SRR22244609 - SRR22244629
Included in the uploaded bioproject are also metabarcoding reads targeting the COI amplicon, which have not been used in the published paper due to bad quality. 

Supplementary files and non-sequencing data can be found the below linked DRYAD repository. 
https://datadryad.org/stash/share/ndIWIGf28mmJ-nUc5HqPEuqPQggpeTYRjIDaMt5imEQ

All code used to analyse and visualise data is stored in this GitHub repository. 

## Water chemistry data processing & Visual Data analysis
- Code for Pearson correlations analysing the relationship between environmental parameters, Mann-Whitney test analysing fish swimming speeds between site 2 and the open-water reference sites can be found here [Statistics_EnvParam&Visuals](https://github.com/lexscience/Klein-2024/blob/main/Statistics_EnvParam%26Visuals.R)
- Code for CTD data processing can be found here 

## Sequence data processing
- Code for amplicon sequence variant (ASV) inference using DADA2 available in [DADA2_Inference](https://github.com/lexscience/Klein-2024/blob/main/DADA2_Inference.R)
- Code for cleaning and removal of contaminant sequences available in [Decontam&Cleaning](https://github.com/lexscience/Klein-2024/blob/main/Decontam%26Cleaning)

Clean data tables as input for data analysis can be accessed via dryad or here https://github.com/lexscience/Klein-2024/blob/main/clean_V3V4.csv for V3V4 SSU rRNA gene bacterial matrix or here https://github.com/lexscience/Klein-2024/blob/main/clean_V4.csv for V4 SSU rRNA gene protists matrix

## Data analysis
- Code to produce Figure 3 - Non-metrical multidimensional scaling (NMDS) [NMDS_Species](https://github.com/lexscience/Klein-2024/blob/main/NMDS_Species) and relative abundance taxonomy charts [Taxonomy](https://github.com/lexscience/Klein-2024/blob/main/Taxonomy) for eDNA metabarcoding datasets is deposited. The accompanying environmental data used for eDNA data anlysis is available here https://github.com/lexscience/Klein-2024/blob/main/samdf.csv

## Contact
Project lead by Shannon Klein [@DrShannonKlein] - shannon.klein@kaust.edu.sa . 
Repository created and curated by [@lexscience] - larissa.fruehe@oceanx.org . 
Feel free to contact us. 

<!-- Optional -->
<!-- ## License -->
<!-- This project is open source and available under the [... License](). -->

<!-- You don't have to include all sections - just the one's relevant to your project -->
