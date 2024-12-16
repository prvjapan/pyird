---
title: 'PyIRD: A Python-Based Data Reduction Pipeline for Subaru/IRD and REACH'
tags:
  - Python
  - astronomy
  - spectroscopy
authors:
  - name: Yui Kasagi
    orcid: 0000-0002-8607-358X
    affiliation: 1 
  - name: Hajime Kawahara
    orcid: 0000-0003-3309-9134
    affiliation: "1, 2" 
  - name: Ziying Gu
    affiliation: 2
  - name: Teruyuki Hirano
    orcid: 0000-0003-3618-7535
    affiliation: "3, 4, 5"
  - name: Takayuki Kotani
    orcid: 0000-0001-6181-3142
    affiliation: "3, 4, 5"
  - name: Masayuki Kuzuhara
    orcid: 0000-0002-4677-9182
    affiliation: "3, 5"
  - name: Kento Masuda
    affiliation: 6

affiliations:
 - name: Institute of Space and Astronomical Science, Japan Aerospace Exploration Agency, 3-1-1 Yoshinodai, Chuo-ku, Sagamihara, Kanagawa, 252-5210, Japan
   index: 1
 - name: Department of Astronomy, Graduate School of Science, The University of Tokyo, 7-3-1 Hongo, Bunkyo-ku, Tokyo 113-0033, Japan
   index: 2
 - name: Astrobiology Center, 2-21-1 Osawa, Mitaka, Tokyo 181-8588, Japan
   index: 3
 - name: Astronomical Science Program, The Graduate University for Advanced Studies, SOKENDAI, 2-21-1 Osawa, Mitaka, Tokyo 181-8588, Japan
   index: 4
 - name: National Astronomical Observatory of Japan, 2-21-1 Osawa, Mitaka, Tokyo 181-8588, Japan
   index: 5
 - name: Department of Earth and Space Science, Osaka University, Osaka 560-0043, Japan
   index: 6
date: 12 December 2024
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary
<!-- 
comments

-->
`PyIRD` is a Python-based pipeline for reducing spectroscopic data obtained with IRD (InfraRed Doppler; @kotani2018) and REACH (Rigorous Exoplanetary Atmosphere Characterization with High dispersion coronagraphy; @kotani2020) on the Subaru Telescope. 
It is designed to process raw images into one-dimensional spectra in a semi-automatic manner. 
Unlike traditional methods, it does not rely on IRAF [@tody1986; @tody1993], a software traditionally used for astronomical data reduction. This approach simplifies the workflow while maintaining efficiency and accuracy.
Additionally, the pipeline includes an updated method for removing readout noise patterns from the detector, enabling efficient extraction of spectra even for faint targets such as brown dwarfs.

# Statement of need

The reduction of high-dispersion spectroscopic data has traditionally been performed using IRAF, one of the most widely used software tools for astronomical data reduction and analysis.
However, the National Optical Astronomy Observatories (NOAO) officially ceased its development and maintenance in 2013.
As a result, there is a growing demand for a modern, flexible solution.

`PyIRD` addresses this need and has already been utilized in several papers [@kawashima2024; @kawahara2024; @tomoyoshi2024]. 


# Key Features

`PyIRD` is designed to perform data reduction semi-automatically by following a general workflow for high-dispersion spectroscopic data reduction, as illustrated in Figure \autoref{fig:reduc_flow}.

![Flowchart of the reduction process for IRD and REACH data. The reduction process follows from top to bottom of this figure. Texts in the grey boxes represent instance names of each reduction step used in `PyIRD`. \label{fig:reduc_flow}](fig/reduc_flowchart.png)

Users can define a set of FITS-format files using a Python class named `FitsSet`, and functions in the `Stream2D` class are applied to generate the one-dimentional spectrum.
Since all functions in `PyIRD` are written in Python rather than IRAF's subset preprocessor language (SPP), the package is easy to develop and maintain.
This also significantly reduces the time required for the reduction process: users only need to execute a single Python script without complex IRAF configuration.
For example, reducing data with `PyIRD` typically takes a few tens of minutes to produce one-dimensional spectra from raw data obtained during a single observing night, compared to approximately half a day with traditional IRAF methods.

Moreover, `PyIRD` achieves a higher level of readout noise pattern removal on the detector.
This feature is particularly important for processing data from faint objects such as brown dwarfs, where the astronomical signal is often comparable in strength to systematic noise.
The dominant noise source is the readout pattern from the H2RG detector used in IRD.
To address this, `PyIRD` models the noise by calculating a median profile for each readout channel and applying a 2D Gaussian Process using `gpkron` [@gpkron2022].
This innovative method effectively mitigates the readout pattern, as shown in Figure \autoref{fig:pattern}, and improve data quality for faint targets.

![(Left) Raw image; (Middle) Readout pattern model created by `PyIRD`; (Right) Pattern-corrected image \label{fig:pattern}](fig/clean_pattern.png)

# Acknowledgements

Y.K. acknowledges support from JST SPRING, Grant Number JPMJSP2104 and JSPS KAKENHI grant No. 24K22912.

# References