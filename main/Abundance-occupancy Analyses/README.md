# README

To run these scripts, you will need to install [Stan](https://mc-stan.org/) and some additional *R* packages.

## Install Stan

Follow directions on this website to install Stan: https://mc-stan.org/cmdstanr/articles/cmdstanr.html#installing-cmdstan-1

## Install CmdStanR

**CmdStanR** is the *R* interface to Stan. **posterior** is a package for inference from the fitted model. Use this code to install both:

```
# we recommend running this is a fresh R session or restarting your current session
install.packages(c("cmdstanr", "posterior"), repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
```

## Install other *R* packages

Install additional *R* packages you'll need using this code:

```
install.packages(c("conflicted", "cowplot", "dplyr", "ggdist", "ggplot2", "readr", "stringr", "tidyr"))
```

## Fit model

The code to fit model is in `01_fit-model.R`. You can run using:

```
source("01_fit-model.R")
```

## Make figure

You can also skip fitting the model and use the saved output to make the figure using `02_make-figure.R`.

```
source("02_make-figure.R")
```
