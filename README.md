## Effective Biophysical Model of Transcription and Translation

### Introduction
This repository contains code for effective synthetic circuit models written in the [Julia](https://www.julialang.org) programming language. The model equations are described in the publication:

[Adhikari et al (2020) Effective Biophysical Modeling of Dynamic Cell Free Transcription and Translation Processes
bioRxiv 139774; doi: https://doi.org/10.1101/139774](https://www.biorxiv.org/content/early/2017/12/28/139774)


### Installation and Requirements
[Julia](https://www.julialang.org) v1.3 or greater must be installed on your machine. In addition, the solution and
analysis of the effective TX/TL models require a few additional packages:

Package | Description
--- | ---
[DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl) | A suite for numerically solving differential equations written in Julia and available for use in Julia, Python, and R
[DiffEqSensitivity.jl](https://github.com/JuliaDiffEq/DiffEqSensitivity.jl) | Encodes sensitivity analysis utilities
[JSON.jl](https://github.com/JuliaIO/JSON.jl) | Package for parsing and printing JSON
[Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl) | Univariate and multivariate optimization solvers
[DataFrames.jl](https://github.com/JuliaData/DataFrames.jl) | Tools for working with tabular data in Julia
[CSV.jl](https://github.com/JuliaData/CSV.jl) | Package for working with delimited files
[Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl) | Package implements a variety of interpolation schemes for the Julia language
[Distributions.jl](https://github.com/JuliaStats/Distributions.jl) | A Julia package for probability distributions and associated functions
[PyPlot.jl](https://github.com/JuliaPy/PyPlot.jl) | Julia interface to the [Matplotlib](https://matplotlib.org) plotting library from Python
[POETs.jl](https://github.com/varnerlab/POETs.jl) | A Julia package that implements the Pareto Optimal Ensemble Techniques (POETs) method for multiobjective optimization

These packages will be installed and compiled the first time you execute the model code (which can take a bit).
To get the model codes, you can download the model repository as a zip file, clone or pull it using the command:

	git pull https://github.com/varnerlab/Biophysical-TXTL-Model-Code

or

	git clone https://github.com/varnerlab/Biophysical-TXTL-Model-Code

In the ``src`` directory there are two subdirectories ``P70-deGFP-model`` and ``P28-cIssrA-model`` which contain
the code for the P70-deGFP model and the P28 negative feedback model. 
