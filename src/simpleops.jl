
using Statistics
using CurveFit

export demean, detrend, samptimes, levelcrossings
"""
`demean(x::AbstractVector)`

Remove the mean component.
"""
function demean(x::AbstractVector)
    xm = mean(x)
    return x .- xm
end

"""
`detrend(x::AbstractVector)`

Remove a linear trend from the data

"""
function detrend(x::AbstractVector)
    # Least squares fit of x = a₁ + a₂⋅t
    t = eachindex(x)
    a₁, a₂ = linear_fit(t, x)
    return x .- (a₁ .+ a₂ .* t)
end

samptimes(::Type{T}, fs::T, N; t0=0) where {T} = range(T(t0), step=T(1)/T(fs), length=N)
samptimes(fs, N; t0=0.0) = range(t0, step=1/fs, length=N)
samptimes(fs, x::AbstractVector{T}; t0=0) where {T} = samptimes(T, fs, length(x); t0=t0)



    
function levelcrossings(x::AbstractVector, xlev=0; up=true)

    # We also want to search always up! So we invert the signal if the search
    # should be down
    if up
        s = 1
    else
        s = -1
    end

    # indices of segments where a crossing was found
    segments = Int[]
    xprev = s * (x[begin] - xlev)
    for i in firstindex(x)+1:lastindex(x)
        xnext = s * (x[i] - xlev)
        δx = xnext - xprev
        if δx > 0
            if xnext*xprev < 0
                push!(segments, i-1)
            elseif xprev==0
                push!(segments, i-1)
            end
        end
        xprev = xnext
    end

    return segments
end

