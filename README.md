# Perple_X

[![Build Status](https://github.com/jadconnolly/Perple_X/workflows/CI/badge.svg)](https://github.com/jadconnolly/Perple_X/actions)


This repository contains source code, data files, and executables for [Perple_X](https://www.perplex.ethz.ch), a suite of Fortran 77[^1] programs for calculating phase diagrams, thermobarometry, manipulating thermodynamic data, and modeling phase fractionation and reactive transport. 

Executables for **Windows**, **Linux**, and **macOS** are available on the [releases page](https://github.com/jadconnolly/Perple_X/releases/). See the [Perple_X installation instructions](https://www.perplex.ethz.ch/perple_x_installation.html).
 
Executables for all modern operating systems are also available via the **Julia** package manager. They can be used directly within Julia or, after extraction, as command-line programs. See the [Julia installation instructions](https://github.com/jadconnolly/Perple_X/blob/main/markdown/julia.md). [StatGeochem](https://github.com/brenhinkeller/StatGeochem.jl) provides a basic Julia interface to Perple_X.

[^1]: Perple_X makes limited use of Fortran 90 features.