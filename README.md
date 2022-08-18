# SelfEnergyAnalysis

*A Julia package to extract the self-energy from Exciton-polariton photoluminescence data.*

This package allows you to apply the self-energy analysis commonly used in condensed matter ARPES experiments to Exciton-Polariton photoluminescence (PL) data. 

## Instalation Guide

Please make sure you have the [Julia programming language](https://julialang.org/downloads/) and [git](https://git-scm.com/downloads) installed. It is also recomended to use [Visual Studio Code](https://code.visualstudio.com/) with its [Julia extenstion](https://github.com/julia-vscode/julia-vscode#installing-juliavs-codevs-code-julia-extension).

To start, create or select a directory to work in and open the git terminal in this directory. Then clone the GitLab repository into your folder by entering:

```
git clone ist-git@git.uwaterloo.ca:QuINLab/Projects/ExcitonPolariton/exciton-polariton-self-energy-analysis/SelfEnergyAnalysis.jl.git
```
You may need to enter the passphrase for your SSH protocol depending on your access level. If you are having issues you should check out the [GitLab documentation](https://git.uwaterloo.ca/help/user/ssh.md). 

Once the `SelfEnergyAnalysis.jl` folder is in your working directory you should start the Julia REPL in vscode. To do this open the directory in vscode, open the command pallet with `Ctrl`+`Shift`+`P` and type `julia`. Find ``Julia: Start REPL`` in the drop down menu and press ``enter``. 

In the Julia REPL enter `]` to open the package manager. Now enter
```
add SelfEnergyAnalysis.jl
```
!!! note
    If this throws an error try:
    ```
    add ./SelfEnergyAnalysis.jl
    ```
This should complete the installation of the package. You may find it useful to add the following packages as as well 
```
add PyPlot
```
```
add Statistics
```
If there are any issues please refer to the Teams walkthrough recording.