export logdecrement

import Polynomials: fit, Polynomial, coeffs

function logdecrement(x::AbstractVector{T}) where {T}

    # Find the maximum value of y
    y = abs.(x)
    idx_max = argmax(y)

    idx = levelcrossings_index(x, 0, bothdirs=true)

    istart = findfirst(>(idx_max), idx)
    nw = 1
    ampl = Point2{T}[]
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
        println(imax, "  ", tmax)
        push!(ampl, Point2{T}(tmax, xmax))
    end

    return ampl

end

    
#function dampingratio(
