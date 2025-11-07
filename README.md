### **SigRescueR**
SigRescueR R package provides a pan-system framework for noise correction and mutational signature identification across various mutation classes and mutational channels.

<p align="center"> <img src="https://github.com/user-attachments/assets/72691df9-2d77-4272-b534-ba5487bec70f" alt="SigRescueR logo" width="600"> </p>

### **Introduction**
SigRescueR is a specialized R package designed to improve the accuracy of mutational signature analyses by removing baseline or background mutational signatures from exposure data. It leverages Bayesian inference to isolate true exposure-associated signals, improving downstream biological interpretation and reducing confounding effects in complex cancer genomics datasets.

This tool is particularly valuable for researchers studying environmental mutagens and cancer etiology, enabling cleaner, more interpretable mutational signatures.

**Features**

- Bayesian inference-based mutational signature cleaning
- Support for complex baseline signature removal
- Flexible model configuration including warm-up, iterations, and chains
- Parallelized computation for faster processing
- Integration with multiple **single-base substitutions** (SBS) mutation classification: SBS96, SBS288, SBS1536
- Versatile application for not only SBS mutations, but also **insertion/deletions (indels)**, and **doublet-base substitutions (DBS)**
- Visualization support for cleaned mutational profiles

### **Installation**
```
install.packages("devtools")
devtools::install_github("nguyenptuan/SigRescueR")
```

### **Usage Example**
This example demonstrates analyzing mutational signatures from human intestinal organoids treated with colibactin versus non-mutagenic controls (dataset from Pleguezuelos-Manzano et al., 2020, Nature).

```
library(SigRescueR)

data("treatment")
data("control")

# Load and setup the model to enable faster parallelized runs
SigRescueR::SigRescueSetup(model="COM-Poisson")

# Prepare input data
treatment <- treatment %>%
  column_to_rownames(var = "MutationType") %>%
  round(0) %>%
  rownames_to_column(var = "MutationType")

control <- control %>%
  column_to_rownames(var = "MutationType") %>%
  sweep(2, colSums(.), "/") %>%
  rownames_to_column(var = "MutationType")

# Define output path and filename
output_path <- "/full/path/to/output_directory/"
filename <- "testdata"
objects <- c("control", "treatment", "model", "output_path", "filename")

# Run Bayesian inference
SigRescueR::SigRescueRun(objects = objects)

# Load results and run analysis
load(paste0(output_path, filename, ".rda"))
SigRescueR::SigRescueAnalyze(res = res)

# Visualize cleaned profiles
SigRescueR::SigRescuePlot(samp.cleaned)
```

### **Documentation**

Detailed package documentation is available in the `docs` folder



### **References and Resources**

- [DOI:10.11309/n12078-218-57901-y](url)


 ### **COPYRIGHT**

Copyright (c) 2025, Peter Nguyen [Zhivagui Lab] All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

### **Contact Information**

Please address any queries or bug reports to Peter Nguyen at [Peter.nguyen@unlv.edu](mailto:peter.nguyen@unlv.edu). Additional support can be provided by Maria Zhivagui at [Maria.zhivagui@unlv.edu](mailto:mara.zhivagui@unlv.edu).
