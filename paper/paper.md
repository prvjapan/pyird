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
Unlike traditional methods, it does not rely on `IRAF` [@tody1986; @tody1993], a software used for astronomical data reduction. This approach simplifies the workflow while maintaining efficiency and accuracy.
Additionally, the pipeline includes an updated method for removing readout noise patterns from raw images, enabling efficient extraction of spectra even for faint targets such as brown dwarfs.

# Statement of need

The reduction of high-dispersion spectroscopic data has traditionally been performed using `IRAF`, one of the most widely used software tools for astronomical data reduction and analysis.
Although the National Optical Astronomy Observatory (NOAO) officially ceased development and maintenance of `IRAF` in 2013, community-based maintenance has continued. 
However, the official IRAF community distribution[^iraf_community] and the Space Telescope Science Institute (STScI)[^stsci] have both recommended that researchers transition away from `IRAF` due to its legacy architecture and lack of institutional support.

[^iraf_community]: IRAF Community Distribution website: https://iraf-community.github.io/
[^stsci]: STScI Newsletter (2018), "Removing the Institute’s Dependence on IRAF: You Can Do It Too": https://www.stsci.edu/contents/newsletters/2018-volume-35-issue-03/removing-the-institutes-dependence-on-iraf-you-can-do-it-too

In recent years, several open-source, Python-based pipelines for the reduction of near-infrared echelle spectrographs have been developed.
Some pipelines utilize `PyRAF`, a Python interface to IRAF, such as `WARP` for the WINERED spectrograph [@Hamano2024], while others, including `PLP` for IGRINS [@Sim2014] and `PypeIt` [@pypeit:joss_pub], do not rely on PyRAF or IRAF-based components.
While these pipelines provide general frameworks or instrument-specific solutions, `PyIRD` is designed to offer a simple pipeline optimized for IRD and REACH data reduction.
Furthermore, recent advances combining adaptive optics with these instruments have enabled high-resolution spectroscopic observations of faint companions orbiting bright main-sequence stars.
To support such observations, `PyIRD` implements improved detector noise reduction to extract high-quality spectra from these faint targets.

Together, these developments underscore the need for actively maintained, scalable, and flexible software for high-resolution spectroscopic data reduction.
`PyIRD` addresses this need by providing a modern, Python-based pipeline and has already been utilized in several studies [@kasagi2025; @kawashima2024; @kawahara2024; @tomoyoshi2024]. 


# Key Features

`PyIRD` performs semi-automatic data reduction by following a general workflow for high-dispersion spectroscopic data, as illustrated in \autoref{fig:reduc_flow}.
It primarily utilizes `Astropy` [@astropy:2022], `NumPy` [@harris2020array], `SciPy` [@2020SciPy-NMeth], and `pandas` [@reback2020pandas].

![Flowchart of the reduction process for IRD and REACH data. The reduction process follows from top to bottom of this figure. Texts in the grey boxes represent instance names of each reduction step used in `PyIRD`. \label{fig:reduc_flow}](fig/reduc_flowchart.png)

To simplify the handling of a large number of input FITS-format files, `PyIRD` introduces a Python class called `FitsSet`.
Once initialized with parameters such as the file IDs and the directory containing those files, `FitsSet` automatically organizes and manages the input files and their metadata.
It also allows users to apply reduction functions collectively to a specified list of FITS IDs, enabling efficient and consistent data processing through the `Stream2D` class.

Since all functions in `PyIRD` are written in Python rather than IRAF's subset preprocessor language (SPP), the package is easy to develop and maintain.
This also significantly reduces the time required for the reduction process: users only need to execute a single Python script without complex IRAF configuration.
For example, reducing data with `PyIRD` typically takes a few tens of minutes to produce one-dimensional spectra from raw data obtained during a single observing night, compared to approximately half a day with traditional IRAF methods.

Moreover, `PyIRD` achieves a higher level of readout noise pattern removal on final results.
This feature is particularly important for processing data from faint objects such as brown dwarfs, where the astronomical signal is often comparable in strength to systematic noise.
The dominant noise source is the readout pattern from the H2RG detector used in IRD.
To address this, `PyIRD` models the noise by calculating a median profile for each readout channel and applying a 2D Gaussian Process using `gpkron` [@gpkron2022].
This innovative method effectively mitigates the readout pattern, as shown in \autoref{fig:pattern}, and improves data quality for faint targets.

![(Left) Raw image; (Middle) Readout pattern model created by `PyIRD`; (Right) Pattern-corrected image \label{fig:pattern}](fig/clean_pattern.png)

# Acknowledgements

Y.K. acknowledges support from JST SPRING, Grant Number JPMJSP2104 and JSPS KAKENHI grant No. 24K22912.
Z.G. acknowledges support from Forefront Physics and Mathematics Program to Drive Transformation (FoPM), a World-leading Innovative Graduate Study (WINGS) Program, the University of Tokyo.

# References