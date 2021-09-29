# Marine Heatwave-Driven Impacts on Fish Assemblages Uncovered Through Archived DNA

Zachary Gold <sup>1 , 2</sup> , Ryan P. Kelly <sup>3</sup>, Andrew Olaf Shelton <sup>2</sup>, Andrew Thompson <sup>4</sup>, Kelly D. Goodwin <sup>5</sup>, Ramon Gallego <sup>2</sup>, Kim Parsons <sup>2</sup>, Luke R. Thompson <sup>5 , 6</sup>,  Dovi Kacev <sup>7</sup>, Paul H. Barber <sup>8</sup>

<sup>1</sup> Cooperative Institute for Climate, Ocean, & Ecosystem Studies, UW, Seattle, WA <br />
<sup>2</sup> Northwest Fisheries Science Center, NMFS/NOAA, Seattle, WA <br />
<sup>3</sup> School of Marine and Environmental Affairs, UW, Seattle, WA <br />
<sup>4</sup> Southwest Fisheries Science Center, NMFS/NOAA, La Jolla, CA <br />
<sup>5</sup> Ocean Chemistry and Ecosystems Division, Atlantic Oceanographic and Meteorological Laboratory, Miami, FL <br />
<sup>6</sup> Northern Gulf Institute, Mississippi State University, Mississippi State, MS <br />
<sup>7</sup> Scripps Institution of Oceanography, UCSD, La Jolla <br />
<sup>8</sup> Department of Ecology and Evolutionary Biology, UCLA, Los Angeles, CA <br />

## Description
This page is dedicated to hosting data and code generated for the manuscript.

Included on this page is
1. **Scripts used to Conduct Analyses**

  **/analysis**

      1. *CalCOFI_results_short_20210928.Rmd* This script does the main analyses in the paper.
      2. *calcofi_metadata_20210907.Rmd* This script organizes metadata.
      3. *Calcofi_satellite_data_2021090-7_.Rmd* This script obtains SST data for analyses.
      4. *Calcofi_edna_vs_morphology_20210908.Rmd* This script formats data for the STAN model.
      5. *Run_Mifish_Joint_Model.R* This script runs the joint STAN model
      6. *Mifish_Joint_Model.stan* This is the STAN model.

  **/decon**

      1. *20210927_Calcofi_decontam_zjg.R* This script processes the "raw" *Anacapa Toolkit* output and runs decontamination to remove poorly sequenced samples and contaminant ASVs
      2. *20210927_Merge_data.Rmd* This script formats the output from decontamination and sums ASVs by taxonomy.
2. **Data**

  **/anacapa_output_20210507**

    1. *combo_q30* *Anacapa Toolkit* Output of MiFish 12S data using the *CRUX* Global Reference database from [Gold et al. 2021](https://onlinelibrary.wiley.com/doi/epdf/10.1111/1755-0998.13450) **Note** *ASV raw file needed for decon is not uploaded due to size*
    2. *fishcard_q30* *Anacapa Toolkit* Output of MiFish 12S data using the *CRUX* Local Reference database from Gold et al. 2021. [See GitHub](https://github.com/zjgold/FishCARD) & [See Dryad](https://doi.org/10.5068/D1H963) **Note** *ASV raw file needed for decon is not uploaded due to size*

  **/CalCOFI_Database_194903-201907_csv_30Apr2020**

    1. Directory includes metadata files **Note** *Bottle data file needed for metadata is not uploaded due to size*

  **/**

    1. *larval_counts_20210305.csv* Microscopy derived morphological data
    2. *mifish_library_prep.csv* Relevant sequence library preparation information needed for the joint model.
    3. *20210622_species_mapping_file.csv* Species mapping file that links taxonomy obtained from morphological and molecular data.
    4. *habitat_association_to_check_art.csv* Habitat association data for all species
    5. Other files are various metadata and intermediate outputs needed to link data to each other.


**Raw sequence inputs, reference databases, and Anacapa scripts  will eventually be  available on Dryad**
