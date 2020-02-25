## Effective Biophysical Model of Transcription and Translation

### Introduction
This repository contains code for effective synthetic circuit models written in the [Julia](https://www.julialang.org) programming language. The model equations are described in the publication:

[Adhikari et al (2020) Effective Biophysical Modeling of Dynamic Cell Free Transcription and Translation Processes
bioRxiv 139774; doi: https://doi.org/10.1101/139774](https://www.biorxiv.org/content/early/2017/12/28/139774)


### Installation and Requirements
[Julia](https://www.julialang.org) v1.3 or greater must be installed on your machine. In addition, solution and
analysis of the effective TX/TL models require a few additional packages:

Package | Description
--- | ---
[DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl) | a suite for numerically solving differential equations written in Julia and available for use in Julia, Python, and R
[DiffEqSensitivity.jl](https://github.com/JuliaDiffEq/DiffEqSensitivity.jl) | encodes sensitivity analysis utilities
[JSON.jl](https://github.com/JuliaIO/JSON.jl) | package for parsing and printing JSON


You can download the model repository as a zip file, clone or pull it using the command:

	git pull https://github.com/varnerlab/Biophysical-TXTL-Model-Code

or

	git clone https://github.com/varnerlab/Biophysical-TXTL-Model-Code
