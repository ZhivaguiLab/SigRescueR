## SigRescueR
A Package For Removing Baseline Mutational Signatures From Exposure Mutational Signature

#### Install _SigRescueR_ from GitHub
```
install.packages("devtools")
devtools::install_github("nguyenptuan/SigRescueR")
```

#### How To Use _SigRescueR_: A Practical Example

In this section, we demonstrate the use of _SigRescueR_ with human intestinal organoids exposed to colibactin 
and non-mutagenic control from Pleguezuelos-Manzano et al. 2020 Nature in SBS-96 classification.
```
data("treatment")
data("control")
```
The first step is to load the model using the `SigRescueSetup()` function. This process may take some time, but it 
enables faster parallelized runs by preventing the model from being reloaded each time the function is executed. 
By default, it sets the number of warm-up, iterations and chains, but you can easily modify these settings 
using the `warmup`, `iter`, and `chains` parameters within the function.
```
SigRescueR::SigRescueSetup(model="COM-Poisson")
```

The next step is to verify that the treatment data consist of integer values, since the function does not accept 
floating-point inputs, and that the control data are normalized to 1.
```
treatment <- treatment %>% column_to_rownames(var = "MutationType") %>%
  round(., digits = 0) %>% 
  rownames_to_column(var = "MutationType")
  
control <- control %>% column_to_rownames(var = "MutationType") %>%
  sweep(x = ., MARGIN = 2, STATS = colSums(.), FUN = "/") %>%
  rownames_to_column(var = "MutationType")
```
Then, specify the `filename` and `output_path`. Next, create a vector of the required objects as `objects` and
run the function `SigRescueRun()` to start the Bayesian inference, which executes the main function of _SigRescueR_. 
It is essential that the `MutationType` column matches between the two objects, as inconsistencies will prevent correct 
execution of the function. The approx. runtime for this example is 5 minutes.
```
output_path <- "/full/path/to/output_directory/"
filename <- "testdata"
objects <- c("control", "treatment", "model", "output_path", "filename")
SigRescueR::SigRescueRun(objects = objects)
```
After completion, the data can be loaded using `output_path` and `filename`, and the `SigRescueAnalyze()` function can be 
executed to generate the cleaned profile object `samp.cleaned` as well as the summary object `samp.res`.
```
load(paste0(output_path, filename, ".rda"))
SigRescueR::SigRescueAnalyze(res = res)
```
Finally, you can use `samp.cleaned` to plot SBS-96 profiles with plotmydata(). By default, the `output_path` will be the 
current working directory with the `filename` will be named "clean_res", but you can easily modify these settings 
using the `output_path` and `filename` parameters within the function
```
SigRescueR::SigRescuePlot(clean = samp.cleaned, output_path = "~/", filename = "clean_res")
```
