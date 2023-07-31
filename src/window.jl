using GLMakie

export WindowExplorer

mutable struct WindowExplorer{Fig,T,TVec,X,TLim}
    fig::Fig
    fs::T
    t::TVec
    x::X
    tmin::TLim
    tmax::TLim
    names::Vector{String}
end

    
function WindowExplorer(fig, fs, x; names=nothing, t0=0, tmin=nothing, tmax=nothing, tlabel="", xlabel="", title="", legend=true) 

    t = range(t0, step=1/fs, length=size(x,1))

    if isnothing(tmin)
        tmin = Observable(minimum(t))
    else
        tmin = Observable(tmin)
    end

    if isnothing(tmax)
        tmax = Observable(maximum(t))
    else
        tmax = Observable(tmax)
    end

    ND = ndims(x)

    if ND==1
        x = hcat(x)
    end
    
        
    if isnothing(names)
        names = string.(1:size(x,2))
    elseif ND==1 && isa(names, Union{AbstractString,Symbol})
        names = [string(names)]
    else
        @assert length(names) == size(x,2) "Number of columns incompatible with number of names!"
    end
    
    ax = Axis(fig[1,1], xlabel=tlabel, ylabel=xlabel, title="")
    for i in 1:size(x,2)
        lines!(ax, t, x[:,i], label=names[i])
    end
    
    vlines!(ax, tmin, color=:blue, label="tmin = $(round(tmin[], digits=3))")
    vlines!(ax, tmax, color=:red, label="tmax = $(round(tmax[], digits=3))")

    interval = IntervalSlider(fig[2,1], range=t, startvalues=(tmin[],tmax[]))

    lift(interval.interval) do ival
        tmin[] = ival[1]
        tmax[] = ival[2]
    end
    if legend
        axislegend(ax)
    end
   
    WindowExplorer(fig, fs, t, x, tmin, tmax, names)
end

WindowExplorer(fig, x::TSeries; tmin=nothing, tmax=nothing, tlabel="", xlabel="",
               title="", legend=true, names=nothing) =
                   WindowExplorer(fig, samplerate(x), x[];
                                  t0=starttime(x), tmin=tmin, tmax=tmax,
                                  tlabel=tlabel, xlabel=xlabel,
                                  title=title, legend=legend, names=names)

Base.display(win::WindowExplorer) = display(win.fig)

function Base.getindex(d::WindowExplorer)
    tmin = d.tmin[]
    tmax = d.tmax[]
    idx = tmin .< d.t .< tmax
    return TSeries(d.fs, d.x[idx,:], tmin)
end

Base.minimum(win::WindowExplorer) = win.tmin[]
Base.maximum(win::WindowExplorer) = win.tmax[]
Base.extrema(win::WindowExplorer) = (win.tmin[], win.tmax[])

