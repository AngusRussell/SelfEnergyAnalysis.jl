"""
    add_tails(a, n)

Adds `n` many copies of the first value to start of the `Vector`; `a`, and `n` many copies of the last value to the end of the `Vector`; `a`.

# Example
```jldoctest
julia> a = [1, 2, 3]
3-element Vector{Int64}:
 1
 2
 3

julia> add_tails!(a, 3)
9-element Vector{Int64}:
 1
 1
 1
 1
 2
 3
 3
 3
 3
```
"""
function add_tails!(a::Vector, n::Int)
    for i = 1:n
        push!(a, a[end])
        pushfirst!(a, a[1])
    end
    return a
end

"""
    remove_tails!(a::Vector, n::Int)

Removes `n` many copies of the first value at the start of the `Vector`, and `n` many copies of the last value at the end of  the `Vector`.

# Example
```jldoctest
julia> a
9-element Vector{Int64}:
 1
 1
 1
 1
 2
 3
 3
 3
 3

julia> remove_tails!(a,3)
3-element Vector{Int64}:
 1
 2
 3
```
"""
function remove_tails!(a::Vector, n::Int)
    for i = 1:n
        popat!(a, length(a))
        popfirst!(a)
    end
    return a
end

"""
    closest_index(x, val)

Searches `x` for the value `val` and returns the number index of the closest element to `val`.

# Example
```jldoctest
julia> closest_index(BK30, 3.442828e7)
2659737

julia> BK30[2659737]
3.4428284217012994e7
```    
"""
function closest_index(x, val)
    ibest = 1
    dxbest = abs(x[ibest]-val)
    for i in eachindex(x)
        dx = abs(x[i]-val)
        if dx < dxbest
            dxbest = dx
            ibest = i
            end
    end
    return ibest
end

"""
    remove_NaN!(a, b)

Takes two `Vectors` `a` and `b`. If `a` has ``NaN`` value it is removed along with the value at the same index in `b`.

# Example
```jldoctest
julia> a = [1, 2, 3, 4, NaN, 6, 7]
7-element Vector{Float64}:
   1.0
   2.0
   3.0
   4.0
 NaN
   6.0
   7.0

julia> b = [1, 2, 3, 4, 5, 6, 7]
7-element Vector{Int64}:
 1
 2
 3
 4
 5
 6
 7

julia> remove_NaN!(a, b)

julia> a
6-element Vector{Float64}:
 1.0
 2.0
 3.0
 4.0
 6.0
 7.0

julia> b
6-element Vector{Int64}:
 1
 2
 3
 4
 6
 7
```
"""
function remove_NaN!(a, b)
    i = 1
    while any(isnan.(a))
        if isnan(a[i])
            popat!(a, i)
            popat!(b, i) 
            i = 1
        
        else
            i += 1
        end
    end
end


"""
remove_NaN!(a)

Takes a `Vector` `a`. If `a` has ``NaN`` value it is removed.

# Example
```
julia> a = [1, 2, 3, 4, NaN, 6, 7]
7-element Vector{Float64}:
   1.0
   2.0
   3.0
   4.0
 NaN
   6.0
   7.0

julia> remove_NaN!(a)

julia> a
6-element Vector{Float64}:
 1.0
 2.0
 3.0
 4.0
 6.0
 7.0
```
"""
function remove_NaN!(a)
    i = 1
    while any(isnan.(a))
        if isnan(a[i])
            popat!(a, i)

            i = 1
        
        else
            i += 1
        end
    end
end

"""
    int(x::Number)

Returns the nearest integral of x as type `Int` that is less than or equal to x.

# Example
```
julia> int(6.66)
6
```
"""
int(x) = floor(Int, x::Number) 

"""
    kCropGUI(k, E)

Plots results from MDC protocol and interactively crops the momentum bounds to remove bad values.

Returns the chosen upper and lower bound values for negative and positive k.
"""
function kCropGUI(k::Vector, E::Vector, Eb)
    not_cropped = true
    skip_factor = 0
    kL_lowerbound = 0
    kL_upperbound = 0
    while not_cropped
        fig()
        plt.plot(k, E)
        plt.plot(k, Eb(k), color = "orange", linestyle = ":")
        plt.legend(["Data", "Bare-band"])
        plt.xlim([k[1], 0])
        if skip_factor == 0

            plt.title("Input lower bound for negative k:", weight = "bold")
            println("Input lower bound for negative k:")
            kL_lowerbound = try parse(Float64,readline())
                catch e
                    @error "Enter a numerical value"
                    continue
                end
            if kL_lowerbound > 0
                kL_lowerbound *= -1
            end
            plt.plot(fill(kL_lowerbound, 1:length(E)), E, color="red", linestyle = "-.")

            plt.title("Input upper bound for negative k:", weight = "bold")
            println("Input upper bound for negative k:")
            kL_upperbound = try parse(Float64,readline())
            catch e
                @error "Enter a numerical value"
                continue
            end
            if kL_upperbound > 0
                kL_upperbound *= -1
            end
            plt.plot(fill(kL_upperbound, 1:length(E)), E, color="red", linestyle = "--")


            println("Happy with negative k crop? (y/n) or (skip) for mirrored")
            answerL = readline()
            if answerL == "y" || answerL =="yes" || answerL =="Yes"


            elseif answerL == "skip" || answerL =="s" || answerL =="Skip"
                println("$kL_lowerbound, $kL_upperbound, $(-kL_lowerbound), $(-kL_upperbound)")
                return kL_lowerbound, kL_upperbound, -kL_upperbound, -kL_lowerbound
                break
            else
                continue
            end
        end
        plt.xlim([0, k[end]])
        plt.title("Input lower bound for positive k:", weight = "bold")
        println("Input lower bound for positive k:")
        kR_lowerbound = try abs(parse(Float64,readline()))
        catch e
            @error "Enter a numerical value"
            skip_factor =1
            continue
        end
        plt.plot(fill(kR_lowerbound, 1:length(E)), E, color="red", linestyle="-.")

        plt.title("Input upper bound for positive k:", weight = "bold")
        println("Input upper bound for positve k:")
        kR_upperbound = try abs(parse(Float64,readline()))
        catch e
            @error "Enter a numerical value"
            skip_factor =1
            continue
        end
        plt.plot(fill(kR_upperbound, 1:length(E)), E, color="red", linestyle="--")
        println("Happy with positive k crop? (y/n)")
        answerR = readline()
        if answerR == "y" || answerR =="yes" || answerR =="Yes"
            if skip_factor == 1
                plt.plot(fill(kL_upperbound, 1:length(E)), E, color="red", linestyle = "--")
                plt.plot(fill(kL_lowerbound, 1:length(E)), E, color="red", linestyle = "-.")
            end
            plt.xlim([k[1], k[end]])
            plt.xticks([kL_lowerbound, kL_upperbound, kR_lowerbound, kR_upperbound])
            println("Happy with crop? (y/n)")
            answer = readline()
            if answer == "y" || answer =="yes" || answer =="Yes"
                plt.close()
                return kL_lowerbound, kL_upperbound, kR_lowerbound, kR_upperbound
                break
            else
                skip_factor = 0
                continue
            end
        elseif answerR  =="n" || answerR =="no" || answerR =="No"
            skip_factor = 1
            continue
        else
            continue
        end
    end
end