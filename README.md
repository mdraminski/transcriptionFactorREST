---
output:
  pdf_document: default
  word_document: default
  html_document: default
---
# transcriptionFactorREST
REST acts as a transcriptional repressor of neuronal genes in non-neuronal cells. In some contexts, REST has been reported as a transcriptional activator. We aimed at defining the role of REST in IDH mutation related phenotype in glioma because REST was 1) described as a factor involved in a blockage of cell differentiation and 2) an important oncogenic factor, and knowing that 3) IDH mutations are oncogenic drivers in glioma that cause significant changes in epigenome, leading to block of differentiation. A pair of REST silenced IDH mutant and IDH WT U87 cell lines were used as a model. Whole genome and transcriptome analyses revealed different patterns of REST binding and its proximal TF motifs in IDH mutated and WT cells, and identified the genes downstream of REST related to ECM organization and cell differentiation. This study shows that the REST role in gliomas is dependent on IDH mutation status. 

This project contains all R scripts used to produce all the results related to the paper:
**'Comprehensive analysis of the transcription factor REST regulatory networks in IDH-mutant and IDH-wild type glioma cells and gliomas'**

## Table of contents
* [Installation](#installation)
* [Usage](#usage)
* [Authors](#authors)
* [License](#license)
* [Acknowledgments](#acknowledgments)

# Installation
To rerun some parts of the analysis please make sure you downloaded input data by running the following bash scripts:

* './data/download.data.sh'
* './RData/RDatadownload.RData.sh'

# Usage
The project cosists of several folders that contain R source code to run some parts of the analysis described in the paper.

## Project Structure

* DiffMeth (Author: Michał Dramiński)
* BW (Author: Bartosz Wojtas)
* SurvivalAnalysis (Author: Adria-Jaume Roura)

# Version
- Version: 1.0.0
- Date: 09.06.2022

# Authors
The paper and related set of scripts has been created and implemented by:

- Małgorzata Perycz [1,2] (author) ()
- Michal J. Dabrowski [2] (author, developer, maintainer)
- Marta Jardanowska [2] (author)
- Adria-Jaume Roura [1] (author, developer, maintainer)
- Bartlomiej Gielniewski [1] (author)
- Karolina Stepniak [1] (author)
- Michał Dramiński [2] (author, developer, maintainer)
- Bartosz Wojtas [1,3] (author, developer, maintainer) 
- Bozena Kaminska [1] (author)

1. Laboratory of Molecular Neurobiology, Nencki Institute of Experimental Biology, Polish Academy of Sciences
2. Computational Biology Group, Institute of Computer Science of the Polish Academy of Sciences, Warsaw, Poland
3. Sequencing Laboratory, Nencki Institute of Experimental Biology, Polish Academy of Sciences 

# License
This project and the accompanying materials are made available under the terms of the GNU Public License v3.0 which accompanies this distribution, and is available at [http://www.gnu.org/licenses/gpl.html](http://www.gnu.org/licenses/gpl.html). For more information please see LICENSE file.

# Acknowledgments

- We would like to thank to everyone who helped and contributed to this work. 
- This work is supported by a grant from the Polish National Science Centre [DEC-2015/16/W/NZ2/00314].
