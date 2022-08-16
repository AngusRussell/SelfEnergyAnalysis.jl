"""
    fig()

Creates [Matplotlib figure](https://matplotlib.org/) object with a 6x6 size and custom formating. Quick way of keeping figure layout consistent

# Example
```
julia> fig()

julia> plt.plot(BK30[:,400,1])
1-element Vector{PyCall.PyObject}:
 PyObject <matplotlib.lines.Line2D object at 0x000000006C26F550>

```
# Example
```
julia> fig()

julia> plt.contourf(BK30[:,:,1])
PyObject <matplotlib.contour.QuadContourSet object at 0x000000006C309310>

```
"""
function fig()
    plt.figure(figsize = (6, 6))
    rc("xtick.major", size = 4)
    rc("xtick.major", width = 3)
    rc("ytick.major", size = 4)
    rc("ytick.major", width = 3)
    rc("xtick",labelsize=17)
    rc("ytick",labelsize=17)
    rc("font", weight="bold")
    rc("axes", labelsize=17, linewidth=3)
    plt.subplots_adjust(bottom=0.2)
    plt.subplots_adjust(left=0.25)
end

"""
    fig(width::Float64, height::Float64)

Creates [Matplotlib figure](https://matplotlib.org/) object with a 6x6 size and custom formating. Quick way of keeping figure layout consistent

# Example
```
julia> fig(4, 4)

julia> plt.plot(BK30[:,400,1])
1-element Vector{PyCall.PyObject}:
 PyObject <matplotlib.lines.Line2D object at 0x000000006C26F550>

```
# Example
```
julia> fig(14, 8)

julia> plt.contourf(BK30[:,:,1])
PyObject <matplotlib.contour.QuadContourSet object at 0x000000006C309310>

```
"""
function fig(width::Number, height::Number)
    plt.figure(figsize = (width, height))
    rc("xtick.major", size = 4)
    rc("xtick.major", width = 3)
    rc("ytick.major", size = 4)
    rc("ytick.major", width = 3)
    rc("xtick",labelsize=17)
    rc("ytick",labelsize=17)
    rc("font", weight="bold")
    rc("axes", labelsize=17, linewidth=3)
    plt.subplots_adjust(bottom=0.2)
    plt.subplots_adjust(left=0.25)
end