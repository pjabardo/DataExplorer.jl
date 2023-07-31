export filter_explorer

using DSP
#=
struct FilterExplorer
    fig::Fig
    fs::T
    filt::Filt
    freqs::Tuple{T,T}
    order::Int
end
=#
        

function makefilt(fs, ftype, freq; order=2)

    if ftype == :bandpass
        return digitalfilter(Bandpass(freq[1], freq[2], fs=fs), Butterworth(order))
    elseif ftype==:bandstop
        return digitalfilter(Bandstop(freq[1], freq[2], fs=fs), Butterworth(order))
    elseif ftype==:lowpass
        return digitalfilter(Lowpass(freq[1], fs=fs), Butterworth(order))
    elseif ftype==:highpass
        return digitalfilter(Highpass(freq[1], fs=fs), Butterworth(order))
    end
    error("Unknow filter type $ftype!")
end

        
function filter_explorer(fig, fs, x;
                         n=div(length(x),8), window=nothing,
                         order=2, ftype=:bandpass)
    
    figts = fig[1,1] = GridLayout()
    figfreq = fig[2,1] = GridLayout()
    
    win = WindowExplorer(figts, fs, x, legend=false)
    
    xts = TSeries(fs, x)
    xpts = Observable(Point2.(sampletimes(xts), xts.x))
    
    S = welch_pgram(x, n, fs=fs, window=window)
    S2 = welch_pgram(x, n, fs=fs, window=window)
    
    slines = Observable(Point2.(freq(S)[2:end], power(S)[2:end]))
    slines2 = Observable(Point2.(freq(S2)[2:end], power(S2)[2:end]))
    
    freqmin = Observable(freq(S)[begin+1])
    freqmax = Observable(freq(S)[end])

    if ftype == :bandpass
        freqmin[] = freq(S)[begin+1]
        freqmax[] = freq(S)[end-1]
    elseif ftype == :bandstop
        freqmin[] = freq(S)[end-1]
        freqmax[] = freq(S)[end-1]
    elseif ftype == :lowpass
        freqmin[] = freq(S)[end-1]
        freqmax[] = freq(S)[begin+1]
    elseif ftype == :higpass
        freqmin[] = freq(S)[begin+1]
        freqmax[] = freq(S)[end-1]
    end
    myfilt = Observable(makefilt(fs, ftype, (freqmin[], freqmax[]); order=order))
    reload = Button(figfreq[1,1], label="Load data section", tellwidth=false)
    
    ax = Axis(figfreq[2,1])
    lines!(ax, slines)
    
    interval = IntervalSlider(figfreq[3,1], range=freq(S), startvalues=(freqmin[],freqmax[]))

    on(reload.clicks) do xxx
        x1 = win[]
        S = welch_pgram(x1.x[:,1], n; fs=fs, window=window)
        slines[] = Point2.(freq(S)[2:end], power(S)[2:end])
    end

    labmin = Observable("$(round(freqmin[], digits=3))")
    labmax = Observable("$(round(freqmax[], digits=3))")
    
    vlines!(ax, freqmin, color=:blue)
    if ftype != :lowpass || ftype != highpass
        vlines!(ax, freqmax, color=:red)
    end
    

    lift(interval.interval) do ival
        freqmin[] = ival[1]
        freqmax[] = ival[2]
        labmin[] = "Freq. min. $(round(freqmin[], digits=3)) (Hz)"
        labmax[] = "Freq. max. $(round(freqmax[], digits=3)) (Hz)"
    end


    figfilt = fig[1:2,2] = GridLayout()
    info = figfilt[1,1] = GridLayout(tellwidth=false)
    loadfreqs = Button(info[1,1], label="Load frequencies")
    Label(info[1,2], labmin)
    Label(info[1,3], labmax)
    
    
    
    axfilt = Axis(figfilt[2,1])
    lines!(axfilt, xpts)

    axspec = Axis(figfilt[3,1], title="Periodogram")
    lines!(axspec, slines2)
    
    on(loadfreqs.clicks) do xxx
        x1 = win[]
        xts1 = TSeries(fs, x)
        myfilt[] = makefilt(fs, ftype, (freqmin[], freqmax[]); order=order)
        x2 = filtfilt(myfilt[], xts1.x)
        xpts[] = Point2.(sampletimes(xts1), x2)
        S2 = welch_pgram(x2, n, fs=fs, window=window) 
        slines2[] = Point2.(freq(S2)[2:end], power(S2)[2:end])
        
    end
    
    
    win,fig
    
end
