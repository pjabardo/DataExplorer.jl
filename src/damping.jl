export logdecrement, vibration_damping

import Polynomials: fit, Polynomial, coeffs

function logdecrement(x::AbstractVector{T}) where {T}

    # Find the maximum value of y
    y = abs.(x)
    idx_max = argmax(y)

    idx = levelcrossings_index(x, 0, bothdirs=true)

    istart = findfirst(>(idx_max), idx)
    nw = 1
    ampl = T[]
    for i in eachindex(idx[istart:end-1])
        k = idx[i]
        kn = idx[i+1]
        yi = y[k:kn]
        imax = argmax(yi)-1 + k
        limlft = max(imax-nw, firstindex(y))
        limrgt = min(imax+nw, lastindex(y))
        ii = limlft:limrgt
        yfit= y[ii]
        p = fit(Polynomial, ii, yfit, 2)
        a₁, a₂, a₃ = coeffs(p)
        tmax = -a₂/2a₃
        xmax = p(tmax)
        push!(ampl, xmax)
    end
    return ampl

end

    
function vibration_damping(x::AbstractVector{T}) where {T}
    x = detrend(x)
    y = abs.(x)
    idx_max = argmax(y)

    idx = levelcrossings_index(x, 0, bothdirs=true)

    N = size(x,1)
    
    t₀ = Float64[]
    
    istart = findfirst(>(idx_max), idx)
    nw = 1
    ampl = T[]
    for i in eachindex(idx[istart:end-1])
        
        k = idx[i]
        kn = idx[i+1]

        push!(t₀, linear_interp((k, k+1), (x[k], x[k+1]), 0))
        
        
        yi = y[k:kn]
        imax = argmax(yi)-1 + k
        limlft = max(imax-nw, firstindex(y))
        limrgt = min(imax+nw, lastindex(y))
        ii = limlft:limrgt
        yfit= y[ii]
        p = fit(Polynomial, ii, yfit, 2)
        a₁, a₂, a₃ = coeffs(p)
        tmax = -a₂/2a₃
        xmax = p(tmax)
        push!(ampl, xmax)
    end
    i = idx[end]
    
    push!(t₀, linear_interp( (i-1, i), (x[i-1], x[i]), 0))
    
    return ampl, t₀

end

