"""
    load_all(File_path, File_type, k_length, λ_length, header)

Load a folder of raw PL data into an `Array`. Returns 3D `Array` of PL values and `Vector` of wavelength values.

Path of folder specified in `File_path` as a `String` with "`/`" between folders and corresponding pixel dimensions `k_length` and `λ_length`.

`header` argument of function skipps the header positions in the CSV. If there are no headers, enter `false`.

# Example
```julia
julia> BK30, λ30 = load_all("BK30KM32d3L/BK30KM32d3L/raw_data", "csv", 400, 1340, 3);

julia> BK40, λ40 = load_all("20220718 FF SPE-BK40KM32d4T4K-20140304", "dat", 400, 1340, false)
```
*Hint*: `;` supresses output
"""
function load_all(File_path::String, File_type::String, k_length::Int, λ_length::Int, header::Int)
    data_only = "*."*File_type
    File_names = glob(data_only, File_path)
    number_of_files = length(File_names)
    Data3D = Array{Float64, 3}(undef, k_length, λ_length, length(File_names))
    λ = Array{Float64, 1}(undef, λ_length)

    for i = 1:number_of_files
        # This code skips text lines in CSV and extracts data, I assume same for all Folders but if not change values bellow accordingly
        DF = CSV.read(File_names[i], drop=[1], header=header, footerskip=4, Tables.matrix);
        Data3D[:,:,i] = DF[2:k_length+1,:];
        λ[:] = vec(DF[1,:]);
    end
    return Data3D, λ
end
function load_all(File_path::String, File_type::String, k_length::Int, λ_length::Int, header::Bool) # this function adds a method for BK40 with no headers and inverted data structure
    data_only = "*."*File_type
    File_names = glob(data_only, File_path)
    number_of_files = length(File_names)
    Data3D = Array{Float64, 3}(undef, k_length, λ_length, length(File_names))
    λ = Array{Float64, 1}(undef, λ_length)

    for i = 1:number_of_files
        # This code skips text lines in CSV and extracts data, I assume same for all Folders but if not change values bellow accordingly
        DF = CSV.read(File_names[i], header=header, Tables.matrix);
        Data3D[:,:,i] = DF[:,2:k_length+1]' # Applys transpose so that the k_index is same as BK30
        λ[:] = vec(DF[:,1]);
    end
    return Data3D, λ
end

"""
    crop_data(Data3D::Array, λ::Vector, filenumber::Integer)

Interactive process for taking the output of `load_all` and cropping to a user designated region of interest.

Takes an input of a 3D PL data `Array`, 1D wavelength `Vector` and then a `filenumber` `Integer`. 
The `filenumber` value selects the pump-power that is plotted to provide a visual aid when selecting pixel values to crop to.
The output returns the cropped PL `Array` and wavelength `Vector` and also prints the chosen cropping pixels so one can use these values in the other `crop_data` method and skip the interactie process.

# Example
```
julia> crop_data(BK40, λ, 19);
Input lower bound for E pixel (between 1:1340):
580 
Input upper bound for E pixel (between 580:1340):
815
Are you happy with the Energy crop? (Answers: yes/y | no/n | exit )
y
Input lower bound for k pixel (between 1:400):
110
Input upper bound for k pixel (between 110:400):
290
Happy with crop? (Answers: yes/y | no/n )
y
Congrats! You have cropped the data!
Take a note of your cropping coordinates to avoid the GUI process next time!
k tuple: (110, 290)
E tuple: (580, 815)
```
See also [`load_all`](@ref)
"""
function crop_data(Data3D::Array, λ::Vector, filenumber::Integer)
    
    not_cropped = true
    E_lowerbound, E_upperbound, k_lowerbound, k_upperbound = 0, 0, 0, 0;
    skip_factor = 0
    E_max = size(Data3D,2)
    k_max = size(Data3D,1)
    while not_cropped
        fig(10,10)
        plt.contourf(Data3D[:,:,filenumber])
        
        if skip_factor == 0
            plt.title("Input lower bound for E pixel (between 1:$E_max):", weight = "bold")
            display(gcf())
            println("Input lower bound for E pixel (between 1:$E_max):")
            E_lowerbound = try parse(Int64,readline())
                catch e
                    @error "Enter an integer value (between 1:$E_max)"
                    continue
                end
            if E_lowerbound < 1 || E_lowerbound > E_max
                @error "Enter an integer value between 1:$E_max"
                continue
            end
            
            skip_factor = 1
        end

        if skip_factor == 1
            plt.plot(fill(E_lowerbound, k_max-1), 1:k_max-1, color = "red")
            plt.title("Input upper bound for E pixel (between $E_lowerbound:$E_max):", weight = "bold")
            display(gcf())
            println("Input upper bound for E pixel (between $E_lowerbound:$E_max):")
            E_upperbound = try parse(Int64,readline())
            catch e
                @error "Enter an integer value (between $E_lowerbound:$E_max)"
                continue
            end
            if E_upperbound < E_lowerbound || E_upperbound > E_max
                @error "Bad input: Enter an integer value between $E_lowerbound:$E_max"
                continue
            end
            plt.plot(fill(E_upperbound, k_max-1), 1:k_max-1, color = "red")
            display(gcf())
            # Chance to redo E crop
            println("Are you happy with the Energy crop? (Answers: yes/y | no/n | exit )")
            answer = readline()
            if answer == "y" || answer == "yes"
                skip_factor = 2
            end
            if answer == "n" || answer == "no"
                skip_factor = 0
                continue
            end
            if answer == "exit"
                break
            end
        end


        if skip_factor == 2
            plt.plot(fill(E_lowerbound, k_max-1), 1:k_max-1, color = "red")
            plt.plot(fill(E_upperbound, k_max-1), 1:k_max-1, color = "red")

            plt.title("Input lower bound for k pixel (between 1:$k_max):", weight = "bold")
            display(gcf())
            println("Input lower bound for k pixel (between 1:$k_max):")
            k_lowerbound = try parse(Int64,readline())
            catch e
                @error "Enter an integer value (between 1:$k_max)"
                continue
            end
            if k_lowerbound < 1 || k_lowerbound > k_max
                @error "Bad input: Enter an integer value between 1:$k_max"
                continue
            end
            skip_factor = 3
        end

        if skip_factor == 3
            plt.plot(fill(E_lowerbound, k_max-1), 1:k_max-1, color = "red")
            plt.plot(fill(E_upperbound, k_max-1), 1:k_max-1, color = "red")
            plt.plot(1:E_max-1, fill(k_lowerbound, E_max-1), color="red")
            plt.title("Input lower bound for k pixel (between $k_lowerbound:$k_max):", weight = "bold")
            display(gcf())
            println("Input upper bound for k pixel (between $k_lowerbound:$k_max):")
            k_upperbound = try parse(Int64,readline())
            catch e
                @error "Enter an integer value (between $k_lowerbound:$k_max)"
                continue
            end
            if k_lowerbound < 1 || k_lowerbound > k_max
                @error "Bad input: Enter an integer value between  $k_lowerbound:$k_max"
                continue
            end
            plt.plot(1:E_max-1, fill(k_upperbound, E_max-1), color="red")
        end

        plt.xticks([E_lowerbound, E_upperbound])
        plt.yticks([k_lowerbound, k_upperbound])
        display(gcf())
        println("Happy with crop? (Answers: yes/y | no/n )")
        answer = readline()
        if answer == "y" || answer == "yes"
            plt.close("all")
            break
        else
            skip_factor = 0
        end
    end
    println("Congrats! You have cropped the data! ")
    println("Take a note of your cropping coordinates to avoid the GUI process next time!")
    println("k tuple: ($k_lowerbound, $k_upperbound)")
    println("E tuple: ($E_lowerbound, $E_upperbound)")

    k_new_length = k_upperbound-k_lowerbound   # Upperbound - lowerbound gives length 252

    E_new_length = E_upperbound-E_lowerbound # Upperbound - lowerbound gives length 820

    # Julia indexing starts from 1
    k_lowerbound -= 1
    k_upperbound -= 1
    E_lowerbound -= 1
    E_upperbound -= 1

    # Define arrays with dimension size you want after the cropping
    Data3D_cropped = Array{Float64, 3}(undef, k_new_length, E_new_length, size(Data3D,3))
    λ_cropped = Array{Float64,1}(undef, E_new_length)

    for k = 1: size(Data3D,3) # iterate through number of 2D arrays
        for i = 1:k_new_length
            for j = 1:E_new_length     # iterate over E pixels values
                Data3D_cropped[i, j, k] = Data3D[k_lowerbound+i, E_lowerbound+j, k] # iterate from lower bound
                λ_cropped[j] = λ[E_lowerbound+j]  # iterate from lower bound
            end
        end
    end
  
    return Data3D_cropped, λ_cropped
end

"""
    crop_data(Data3D, λ, (k_bounds)::Tuple, (E_bounds)::Tuple)

Crop raw `Array` of PL data, `Data3D`, to region of interest defined by the `k_bounds` and `E_bounds` tuples and then reduces wavelength `Vector`, `λ`, accordingly.

# Example
```
julia> BK30, λ = load_all("BK30KM32d3L/BK30KM32d3L/raw_data", 400, 1340, 3);

julia> BK30_cropped, λ_cropped = crop_data(BK30, λ, (42,294), (350,1170))
([607.0 607.0 … 607.0 607.0; 605.0 605.0 … 607.0 607.0; … ; 605.0 602.0 … 603.0 603.0; 605.0 603.0 … 603.0 603.0]

[606.0 604.0 … 605.0 604.0; 606.0 606.0 … 605.0 605.0; … ; 605.0 602.0 … 604.0 604.0; 605.0 605.0 … 604.0 604.0]

[599.0 599.0 … 602.0 602.0; 599.0 598.0 … 602.0 602.0; … ; 599.0 599.0 … 598.0 598.0; 599.0 598.0 … 598.0 598.0]

...

[602.0 602.0 … 604.0 601.0; 601.0 602.0 … 604.0 604.0; … ; 605.0 603.0 … 605.0 606.0; 603.0 601.0 … 605.0 607.0]

[607.0 606.0 … 604.0 605.0; 601.0 601.0 … 601.0 603.0; … ; 599.0 599.0 … 604.0 606.0; 599.0 598.0 … 606.0 606.0]

[602.0 603.0 … 598.0 600.0; 602.0 602.0 … 598.0 597.0; … ; 598.0 600.0 … 604.0 602.0; 599.0 599.0 … 604.0 602.0], [779.0253992871895, 779.034742820494, 779.0440860556898, 779.0534289927556, 779.0627716316701, 779.072113972412, 779.08145601496, 779.090797759293, 779.1001392053894, 779.109480353228  …  786.4941406902949, 786.5032359412629, 786.512330877488, 786.5214254989506, 786.5305198056308, 786.5396137975085, 786.5487074745646, 786.5578008367789, 786.566893884132, 786.5759866166042])
```
See also [`load_all`](@ref)
"""
function crop_data(Data3D::Array, λ::Vector, k_bounds::Tuple, E_bounds::Tuple)
    k_lowerbound = k_bounds[1] -1       # 42 but indexing starts from 1 in julia
    E_lowerbound = E_bounds[1] -1     # 350 but indexing starts from 1 in julia

    k_new_length = k_bounds[2]-k_bounds[1]
    E_new_length = E_bounds[2]-E_bounds[1]

    # Define arrays with dimension size you want after the cropping
    Data3D_cropped = Array{Float64, 3}(undef, k_new_length, E_new_length, size(Data3D,3))
    λ_cropped = Array{Float64,1}(undef, E_new_length)

    for k = 1: size(Data3D,3) # iterate through number of 2D arrays
        for i = 1:k_new_length
            for j = 1:E_new_length     # iterate over E pixels values
                Data3D_cropped[i, j, k] = Data3D[k_lowerbound+i, E_lowerbound+j, k] # iterate from lower bound
                λ_cropped[j] = λ[E_lowerbound+j]  # iterate from lower bound
            end
        end
    end
    return Data3D_cropped, λ_cropped
end

"""
    attenuation_correction_scalar(Data3D::Array, saturated_filenumber::Int, unsaturated_filenumber::Int)

Calculate the scalar correction factor by comparing two PL 2D `Arrays` with the same pump-power but different attenuation. 

The attenuation can come from either a change in exposure time or a change in the optical density (OD). Shorter exposure time would be the unstaurated, attenuated dataset. Higher OD number would be the unsaturated, attenuated dataset.

Using the information in the filename the user can determine the saturated and unsaturated data sets.

# Example
Looking at the BK30 data set:

"BK30KM32d3L\\BK30KM32d3L\\raw_data\\P000.1mW774.6nmG1800t1.00sOD0.0FFPos301000.csv"

"BK30KM32d3L\\BK30KM32d3L\\raw_data\\P000.2mW774.6nmG1800t1.00sOD0.0FFPos301000.csv"

"BK30KM32d3L\\BK30KM32d3L\\raw_data\\P000.4mW774.6nmG1800t0.10sOD0.0FFPos301000.csv" 

"BK30KM32d3L\\BK30KM32d3L\\raw_data\\P000.4mW774.6nmG1800t1.00sOD0.0FFPos301000.csv"  

"BK30KM32d3L\\BK30KM32d3L\\raw_data\\P000.8mW774.6nmG1800t0.10sOD0.0FFPos301000.csv"

⋮

"BK30KM32d3L\\BK30KM32d3L\\raw_data\\P007.2mW774.6nmG1800t0.10sOD5.0FFPos301000.csv" 

One can see that the third and fourth file have the same pump power but different exposure times. Therefore to calculate the correction factor:

```
julia> BK30, λ = load_all("BK30KM32d3L\\BK30KM32d3L\\raw_data", 400, 1340, 3);

julia> BK30, λ = crop_data(BK30, λ, (42,294), (350, 1170));

julia> exposure_correction = attenuation_correction_scalar(BK30, 4, 3)
8.198551113496109
```
"""
function attenuation_correction_scalar(Data3D::Array, saturated_filenumber::Int, unsaturated_filenumber::Int)
    unsaturated_I_ratio = Array{Union{Missing, Float64}, 2}(missing , size(Data3D,1), size(Data3D,2))

    saturated = Data3D[:,:,saturated_filenumber]
    unsaturated = Data3D[:,:,unsaturated_filenumber]
    maximum_sat_I = maximum(saturated)

    for i = 1:size(Data3D,1)
        for j = 1:size(Data3D,2)
            if 0.25*maximum_sat_I < saturated[i,j] < 0.75*maximum_sat_I
                unsaturated_I_ratio[i,j] = saturated[i,j] / unsaturated[i,j]
            end
        end
    end
    scalar = median(skipmissing(unsaturated_I_ratio))
    return scalar
end

"""
    remove_background(Data3D, Background)

Takes `Array` of cropped PL data and sets all intensity values bellow user defined `Background` value to zero.

# Example
```
BK30 = remove_background(BK30_cropped, 658)
252×820×17 Array{Float64, 3}: ...
```
"""
function remove_background(Data3D::Array, Background::Number)
    Data3D = Data3D .- Background
    for k = 1:size(Data3D,3)
        for i = 1:size(Data3D,1)
            for j = 1:size(Data3D,2)
                if Data3D[i,j, k] < 1
                    Data3D[i,j,k] = 0
                end
            end
        end
    end
    return Data3D
end

"""
    correct_BK30()

Function specific to the BK30KM32d3L data set taken recently at the University of Waterloo. 

This function loads, crops and corrects the BK30 data provided the folder is in the same path as current working directory.
User can look to this function as a guide of how to combine functions `load_all`, `crop_data` and `attenuation_correction_scalar` in order to create their own function for processing PL datasets.

# Example
```
julia> BK30, λ30 = correct_BK30("<path_to_raw_data_here>");
```
"""
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

"""
    correct_BK40(Path)

Function specific to the SPE-BK40KM32d4T4K-20140304 data set taken recently at the Stanford University. 

This function loads, crops and corrects the BK40 data provided the folder is in the same path as current working directory.
User can look to this function as a guide of how to combine functions `load_all`, `crop_data` and `attenuation_correction_scalar` in order to create their own function for processing PL datasets.

# Example
```
julia> BK40, λ40 = correct_BK40("<path_to_raw_data_here>");
```
"""
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

"""
    pixel_to_k(Data2D)

Takes 2D `Array` of corrected PL data and finds k=0 pixel. Returns a pixel to k conversion function and k0 pixel value.

This function is also written in a way specific to the BK30KM32d3L data set with regards to the pixel range to search in. May need ot be changed when applying to different data sets.

# Example
```
julia> k, k0 = pixel_to_k(BK30[:,:,2])
(k, 127.3497705660506)

julia> k
k (generic function with 1 method)

julia> k(1:252)
-5.331960317887335e6:42200.0:5.260239682112665e6
```
"""
function pixel_to_k(Data2D::Array)
    L = size(Data2D,1)
    x = collect(1:L)
    LB = int(L*0.5 -4)
    UB = int(L*0.5 +4)

    max_x_index = argmax(Data2D[LB:UB,:])[1] + LB-1 # These numbers define te pixel window within which we look for 
    maxIdata = maximum(Data2D[LB:UB,:])
    max_y = argmax(Data2D[LB:UB,:])[2] 
    y = Data2D[:, max_y]

    
    param = Array{Float64, 1}(undef, 3)
    FWHM_guess = 8
    Lorentzian(x, p) =  (1/pi .* (0.5.*p[1])./((x.-p[2]).^2 .+ (0.5.*p[1]).^2)).*p[3] 
    while true
        p0 = [FWHM_guess, max_x_index, maxIdata/2]
        fit = curve_fit(Lorentzian, x, y, p0)
        param = fit.param
        error_percent = stderror(fit) ./ abs.(param) .* 100
        if error_percent[2] < 5 || FWHM> 50
            break
        end
        FWHM_guess +=0.5
    end
    k0 = param[2]

    k(x) = (x.-k0).* 4.22e4 # pixel to m^-1 conversion factor
    return k, k0
end


"""
    λ_to_E(wavelength)

Converts `Vector` containing wavelength values to `Vector` containing energy values in units of eV.

# Example
```
julia> E = λ_to_E(λ)
820-element Vector{Float64}:
 1.5926313998361783
 1.5926122982432607
 1.592593197717957
 1.5925740982602512
 1.592554999870126
 1.5925359025475647
 1.5925168062925505
 1.5924977111050667
 1.5924786169850962
 ⋮
 1.5774708971065752
 1.5774526566617348
 1.577434417270037
 1.577416178931463
 1.5773979416459927
 1.577379705413609
 1.5773614702342915
 1.5773432361080215
```
"""
function λ_to_E(wavelength)
    E = Vector{Float64}(undef,length(wavelength))

    h = 6.62607015e-34 # J Hz
    e = 1.60217663e-19 # C
    c = 3e8
    i=1
    for λ in wavelength
        E[i] = ((h*c)/(λ*1e-9)/e)
        i+=1
    end
    return E
end

"""     Plot_all(Data3D)

Plots all of the loaded raw data. This is useful for identifying the file
to use in the interactive cropping process, i.e. the third argument in [`crop_data`](@ref).

# Example 
```julia
julia> Plot_all(BK31_raw)
```
"""
function Plot_all(Data3D)
    for i = 1:size(Data3D, 3)        
        fig()
        plt.contourf(Data3D[:,:,i])
    end
end

"""
    MirrorData(Data2D)

Takes a 2D array of PL data and mirrors around the ``k=0`` point determined by [`pixel_to_k`](@ref). 
Prints the ``k=0`` pixel and returns the mirrored LHS of the original data and the mirrored RHS of the original data.

# Example
```julia
julia> BK31, λ31 = correct_BK31();

julia> MirrorData(BK31[:,:,1]);
203
```
"""
function MirrorData(Data2D)
    k, k0 = pixel_to_k(Data2D)
    centre = int(k0)
    println(centre)
    fig()
    plt.contourf(Data2D)
    plt.plot(1:size(Data2D,2)-1, fill(centre, size(Data2D,2)-1), linestyle="--", color="red")
    plt.yticks([centre])

    mirrored_left = Array{Float64, 2}(undef, centre *2, size(Data2D,2))
    mirrored_right = Array{Float64, 2}(undef, (size(Data2D,1)-centre) * 2, size(Data2D,2))
    a = (size(Data2D,1)-centre)

    # Mirror Left
    mirrored_left[1:centre,:] = Data2D[1:centre,:]
    mirrored_left[end:-1:centre+1,:] = Data2D[1:centre,:]

    mirrored_right[1:a+1, :] = Data2D[end:-1:centre,:]
    mirrored_right[a:end,:] = Data2D[centre:end,:]
    
    fig()
    plt.contourf(mirrored_left)
    plt.plot(1:size(mirrored_left,2)-1, fill(centre, size(mirrored_left,2)-1), linestyle="--", color="red")


    fig()
    plt.contourf(mirrored_right)
    plt.plot(1:size(mirrored_right,2)-1, fill(a, size(mirrored_right,2)-1), linestyle="--", color="red")
        
    return mirrored_left, mirrored_right
end

"""
    check_symmetry(Data2D, wav)

Takes a 2D array of PL data, `Data2D`, and a corresponding wavelength `Vector`, `wav`. It then uses the [`MirrorData`](@ref) function to mirror the data around the ``k=0`` pixel
and then applys the [`energy_band_max`](@ref) to the original, LHS mirrored and RHS mirrored data sets. 

Returns a plot of the three energy bands plotted together. This allows for a visual assessment of the symmmetry around the ``k=0`` pixel.

# Example 
```julia
julia> check_symmetry(BK30[:,:,1], λ30)
203
PyObject <matplotlib.legend.Legend object at 0x000000008C7642B0>
```
"""
function check_symmetry(Data2D, wav)
    mirrored_left, mirrored_right = MirrorData(Data2D)

    Em, kpoints = energy_band_max(Data2D, wav, false)

    E_mleft, k_mleft = energy_band_max(mirrored_left, wav, false)

    E_mright, k_mright = energy_band_max(mirrored_right, wav, false)

    fig()
    plt.scatter(k_mright, Em(k_mright), s=10)
    plt.scatter(k_mright, E_mleft(k_mright), s=10)
    plt.scatter(k_mright, E_mright(k_mright), s=10)
    plt.xlabel(L"k~(\mathrm{m^{-1}})")
    plt.ylabel(L"E~(\mathrm{eV})")
    plt.legend(["EDC of Data", "LHS mirrored", "RHS mirorred"])

end