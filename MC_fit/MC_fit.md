# MC_fit

# Table 1

| Method | Local Eq | Bulk-Eq | Missing Components | Misfit | Bayes |
|-------:|--|---------|--------------------|--------|-------|
| Method | Local Eq | Bulk Eq | Missing Components | Misfit | Bayes |


More information can be found on the [Perple_X](https://www.perplex.ethz.ch) webpage.

We also provide precompiled binaries for all modern operating systems, which can be downloaded through the julia package manager.
The binaries can be used from the command-line, or directly from within julia.  

## 1. Downloading Perple_X binaries
Follow these steps:
1) Download [julia](https://julialang.org). Easiest is if you download the [pre-compiled binaries](https://julialang.org/downloads/).

2) Start julia, and go to the package manager:
```julia
julia> ]
(@v1.7) pkg> 
```
(note that the version numbers may differ if you downloaded a more recent version of julia).

3) Download Perple_X:
```julia
(@v1.7) pkg> add Perple_X_jll
    Updating registry at `~/.julia/registries/General.toml`
   Resolving package versions...
   Installed Perple_X_jll ─ v6.9.1+0
  Downloaded artifact: Perple_X
    Updating `~/.julia/environments/v1.7/Project.toml`
  [9c18f8f7] + Perple_X_jll v6.9.1+0
    Updating `~/.julia/environments/v1.7/Manifest.toml`
  [9c18f8f7] + Perple_X_jll v6.9.1+0
Precompiling project...
  5 dependencies successfully precompiled in 44 seconds (478 already precompiled)
```
This will download the most recent binaries for your system (note that the version may differ on your system, depending on when you read this).

4) If you, at a later stage, want to *upgrade* the binaries to the latest version do this:
```julia
(@v1.7) pkg> update Perple_X_jll
```

5) You can also install a specific version of Perple_X:
```julia
(@v1.7) pkg> add Perple_X_jll@6.9.1
```
6) Return to the REPL (command-window) of julia by using the `<backspace>` button.

## 2. Running the Perple_X binaries from julia
You can call any of the Perple_X binaries directly from within julia. 
For example, you can run `build` with: 
```julia
julia> using Perple_X_jll
julia> run(Perple_X_jll.build())

Perple_X version 6.9.1, source updated October 14, 2022.

Copyright (C) 1986-2022 James A D Connolly <www.perplex.ethz.ch/copyright.html>.


NO is the default (blank) answer to all Y/N prompts


Enter a name for this project (the name will be used as the
root for all output file names) [default = my_project]:
```

All other Perple_X programs can be accessed in the same way, e.g.:
```julia
julia> run(Perple_X_jll.vertex())
julia> run(Perple_X_jll.psvdraw())
```
etc.

Generally you will not be in your working directory if you start julia. An easy way to change your directory is by using the build-in terminal, which you can reach by typing `;` in the REPL:
```julia
julia> ;
shell> 
```
Once you are in there, you can change directories just as you would in the windows or linux shell.

Note that this simply runs the executables, but you won't have the data back into julia (for plotting etc.). If you are interested in that see below.

## 3. Using the downloaded binaries from the terminal
Alternatively, you can also use the downloaded binaries directly from the terminal/shell. 
Yet, where on Earth did it put those files? 
There is a simple way to determine that:
```julia
julia> using Perple_X_jll
julia> Perple_X_jll.PATH_list
1-element Vector{String}:
 "/Users/kausb/.julia/artifacts/b04b94c4088232c1ff10e3536d0edb88eaa54ce4/bin"
julia> ;
shell> ls /Users/kausb/.julia/artifacts/b04b94c4088232c1ff10e3536d0edb88eaa54ce4/bin
actcor  build   convex  ctransf fluids  frendly meemum  pspts   pssect  pstable psvdraw pt2curv vertex  werami
```

Yet, in addition to the binaries itself, you will also have to tell your system where some required dynamic libraries (such as libgfortran) are located.
```julia
julia> Perple_X_jll.LIBPATH_list
2-element Vector{String}:
 "/Applications/Julia-1.7.app/Contents/Resources/julia/lib/julia"
 "/Users/kausb/.julia/artifacts/b04b94c4088232c1ff10e3536d0edb88eaa54ce4/lib"
```
On your system, you will have to set the correct path (note: they will *not* be the same as shown above, as they are generally different on every system).
Unsurprisingly, things are a bit different on windows *vs* linux/mac:

#### *linux/mac*
```
$ export DYLD_LIBRARY_PATH=/Applications/Julia-1.7.app/Contents/Resources/julia/lib/julia:$DYLD_LIBRARY_PATH
$ export PATH=/Users/kausb/.julia/artifacts/b04b94c4088232c1ff10e3536d0edb88eaa54ce4/bin:$PATH
```
Once this is done, you should be able to call Perple_X files:
```
$ build
```

#### *windows*
```
$ set DYLD_LIBRARY_PATH=:$DYLD_LIBRARY_PATH
$ set PATH=/Users/kausb/.julia/artifacts/b04b94c4088232c1ff10e3536d0edb88eaa54ce4/bin;/Applications/Julia-1.7.app/Contents/Resources/julia/lib/julia;%PATH%
```
The ending for executables is a bit different on windows:
```
$ build.exe
```


Note that if you upgrade the version of Perple_X, it will put the files in a different directory so you'll have to redo this step.


## 4. A julia interface for Perple_X
A basic julia interface to Perple_X is provided in the [StatGeochem](https://github.com/brenhinkeller/StatGeochem.jl) package by Brenhin Keller.



