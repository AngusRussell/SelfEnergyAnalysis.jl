# Processing raw data

!!! note

    Make sure you have correctly installed the package from GitLab following the ReadMe.md on the [GitLab page](https://git.uwaterloo.ca/QuINLab/Projects/ExcitonPolariton/exciton-polariton-self-energy-analysis/SelfEnergyAnalysis.jl). 

    It is recomeneded to have the data folder you want to analysise in the working directory.

Assuming you have followed the instalation guide correctly you should now create a new Julia file and load the following packages:
```julia
using SelfEnergyAnalysis
using PyPlot
using Statistics
```
## Loading the data 

The next step is to load the data you want to analyse. To do this use the [`load_all()`](@ref) function.

!!! note
    Look at the format of the data file before you use load. Below is two examples with different file formats. Also Julia cannot parse `\`, one needs to use either `\\`, or `/` in the file path string.

For csv with 3 rows of headers:

```julia
BK30, λ30 = load_all("BK30KM32d3L/BK30KM32d3L/raw_data", "csv", 400, 1340, 3);
```

For .dat file with no headers:

```julia
BK40, λ40 = load_all("20220718 FF SPE-BK40KM32d4T4K-20140304", "dat", 400, 1340, false)
```
## Initial plotting
It can be useful to plot all the raw data and look through to see if there are any artifacts in your data you want to crop out. We can do this using the [`Plot_all()`](@ref) function:
```julia
Plot_all(BK30)
```
Look through the plots and then you can choose a problem dataset for the next stage.

!!! hint
    If the plots do not display then run: 
    ```julia
    pygui(true)
    ``` 

## Cropping the data

Now the data has been loaded we need to crop it to our region of interest (ROI). We have two methods for doing this with the [`crop_data()`](@ref) function. 

```julia
crop_data(BK40, λ, 19);
```

This first method will open an interactive cropping process where you will be able to see a plot of the pump-power file determined by the third input argument, in our example here it is `19`. This method is useful when it is your first time encountering the data and you want to carefully choose your lower and upper bounds.

The second method allows you to skip the interactive process in method one by having an input for values you decided to crop at:

```julia 
BK30_cropped, λ_cropped = crop_data(BK30, λ, (42,294), (350,1170));
```
## Correcting the data

The raw data needs to be corrected for changes in attenuation. As the intensity of light emitted from the sample increases, the light needs to be attenuated in order to protect the detector and avoid saturation. In your raw data folder, look for filenames that indicate the same laser pump-power but in different files e.g. 

```P000.4mW774.6nmG1800t0.10sOD0.0FFPos301000.csv``` 

```P000.4mW774.6nmG1800t1.00sOD0.0FFPos301000.csv```

In the above example it can be seen that both files have the same pump power of $P=0.4~(\mathrm{mW})$ but the exposure time is different, $t=0.1~(\mathrm{s})$ and $t=1.0~(\mathrm{s})$ respectively. The lower exposure time is our unsaturated data set and the longer exposure time is our saturated. We can calculate a correction factor to apply to all the data sets with $t=0.1~(\mathrm{s})$ exposure time using the appropriately named [`attenuation_correction_scalar()`](@ref) function.

However, before we apply this correction it is necessary to remove any background noise and set it equal to zero. We do this using the [`remove_background`](@ref) function. 

```julia 
background_value = median(BK30_raw)
BK30, λ30 = remove_background(BK30_cropped, background_value);
```

This will set all values less than or equal to the `background_value` to 0. After this process is complete we can apply our correction scalars. 

## Correction function 

We can combine all the above steps into one Julia function which will then be availble to use whenever we want to analyse that particular set of data.
Two examples of how to do this is given below for the **BK30** and **BK40** datasets. 

### Correcting BK30KM32d3L
```julia
function correct_BK30(Path)
    BK30_raw, λ_raw = load_all(Path, "csv", 400, 1340, 3) #Relative path so have BK30KM32d3L folder in working directory, or alternatively use full path
    BK30_cropped, λ_cropped = crop_data(BK30_raw, λ_raw, (42,294), (350,1170))

    # Important: calculate scalar correction factors before removing background.
    exp_scalar = attenuation_correction_scalar(BK30_cropped, 4, 3)
    OD1_scalar =  attenuation_correction_scalar(BK30_cropped, 8, 9)
    OD2_scalar = attenuation_correction_scalar(BK30_cropped, 14, 15)
    OD5_scalar = attenuation_correction_scalar(BK30_cropped, 16, 17)
    BK30 = remove_background(BK30_cropped, 658) # Remove background


    # Define array with 3rd index corresponding to number of files without any duplicates due to exposure correction
    BK30_corrected = Array{Float64, 3}(undef, size(BK30,1), size(BK30,2) , 13) 

    BK30_corrected[:,:,1] = BK30[:,:,1]
    BK30_corrected[:,:,2] = BK30[:,:,2]
    BK30_corrected[:,:,3] = BK30[:,:,3] .* exp_scalar
    BK30_corrected[:,:,4] = BK30[:,:,5] .* exp_scalar
    BK30_corrected[:,:,5] = BK30[:,:,6] .* exp_scalar
    BK30_corrected[:,:,6] = BK30[:,:,7] .* exp_scalar
    BK30_corrected[:,:,7] = BK30[:,:,9] .* exp_scalar .* OD1_scalar
    BK30_corrected[:,:,8] = BK30[:,:,10] .* exp_scalar .* OD1_scalar
    BK30_corrected[:,:,9] = BK30[:,:,11] .* exp_scalar .* OD1_scalar
    BK30_corrected[:,:,10] = BK30[:,:,12] .* exp_scalar .* OD1_scalar
    BK30_corrected[:,:,11] = BK30[:,:,13] .* exp_scalar .* OD1_scalar
    BK30_corrected[:,:,12] = BK30[:,:,15] .* exp_scalar .* OD1_scalar .* OD2_scalar
    BK30_corrected[:,:,13] = BK30[:,:,17] .* exp_scalar .* OD1_scalar .* OD2_scalar .* OD5_scalar

    return BK30_corrected, λ_cropped
end

```
### Correcting BK40KM32d4T4K
```julia
 
function correct_BK40(Path)
    BK40_raw, λ_raw = load_all(Path, "dat", 400, 1340, false)
    BK40_cropped, λ_cropped = crop_data(BK40_raw, λ_raw, (1,400), (481,815))
    file_index = [4 3; 12 11; 19 20; 22 23; 25 26; 29 30]
    BK40_corrected = Array{Float64, 3}(undef, size(BK40_cropped,1), size(BK40_cropped,2), size(BK40_cropped,3)-size(file_index,1))
    Scalar_V = Vector{Float64}(undef, size(file_index,1))
    for i = 1:size(file_index,1)
        Scalar_V[i] = attenuation_correction_scalar(BK40_cropped, file_index[i,1], file_index[i,2])
    end
    BK40 = remove_background(BK40_cropped, 1)

    BK40_corrected[:,:,1:2] = BK40[:,:,1:2]
    BK40_corrected[:,:,3] = BK40[:,:,3] .* Scalar_V[1]
    BK40_corrected[:,:,4:9] = BK40[:,:,5:10].* Scalar_V[1]
    BK40_corrected[:,:,10] = BK40[:,:,11].* Scalar_V[1] .* Scalar_V[2]
    BK40_corrected[:,:,11:16] = BK40[:,:,13:18].* Scalar_V[1] .* Scalar_V[2]
    BK40_corrected[:,:,17:18] = BK40[:,:,20:21].* Scalar_V[1] .* Scalar_V[2] .* Scalar_V[3]
    BK40_corrected[:,:,19:20] =  BK40[:,:,23:24].* Scalar_V[1] .* Scalar_V[2] .* Scalar_V[3] .* Scalar_V[4]
    BK40_corrected[:,:,21:23] =  BK40[:,:,26:28].* Scalar_V[1] .* Scalar_V[2] .* Scalar_V[3] .* Scalar_V[4] .* Scalar_V[5]
    BK40_corrected[:,:, 24] = BK40[:,:,30] .* Scalar_V[1] .* Scalar_V[2] .* Scalar_V[3] .* Scalar_V[4] .* Scalar_V[5] .* Scalar_V[6]

    return BK40_corrected,λ_cropped
end
```
## Checking symmetry
If you would like to check the symmetry of your measured data, one can use the [`check_symmetry()`](@ref). It is probably a good idea to use the lowest pump-power to compare:
```julia
check_symmetry(BK30_cropped[:,:,1], λ_cropped)
```
This produces four plots
1) A plot of the `BK30_cropped[:,:,1]` data with a red-line indicating the $k=0$ pixel position.
2) A plot of the left-hand-side (LHS) of `BK30_cropped[:,:,1]` relative to the $k=0$ pixel, mirrored on to the right-hand-side (RHS), with red line indicating $k=0$ pixel position.
3) A plot of the RHS mirrored on to the left like 2.
4) A scatter plot of the [`energy_band_max()`](@ref) process applied to the three plots above. 

This enables one to look at the difference between LHS and RHS.

## Converting Pixels to $k$ and $E$
Almost all of the functions will convert the data into appropriate physical values like energy $E~(\mathrm{eV})$ or $\mathrm{(meV)}$ and momentum $k~(\mathrm{m^{-1}})$. This section is included so that one knows how to do this themselves incase they need it for their own functions.

To convert the wavelength data vector into energy values the package uses the function [`λ_to_E(λ)`](@ref). An example:
```julia
E = λ_to_E(λ30)
``` 
To convert the pixels from the PL data into momentum we use the [`pixel_to_k()`](@ref). This looks at the center of the PL data and finds the energy pixel with the highest intensity in that narrow domain in the center. It then fits a Lorentzian to the momentum distribution curve (MDC) for that energy value. The centroid of the lorentzian is determined as the $k=0$ pixel which is then plugged into the conversion equation for pixel to momentum:
$$ k= (\mathrm{pixels} - \mathrm{pixel}(k_0))\cdot 4.22\times10^4 $$

This method works well most of the time and has so far proved fairly robust against artifacts in the data . However, if one plots the energy dispersion produced by the EDC protocol and feels as though the k=0 pixel is off it can be corrected using [`shift_E()`](@ref) function as described in the [EDC](@ref) section.