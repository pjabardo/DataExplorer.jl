export TSeries, samplerate, starttime, sampletimes

struct TSeries{N,X,XA<:AbstractArray{X,N},T} <: AbstractArray{X,N}
    fs::T
    x::XA
    t0::T
end

TSeries(fs::T, x; t0=0) where {T} = TSeries(fs, x, T(t0))

Base.getindex(x::TSeries) = x.x
Base.getindex(x::TSeries,idx...) = x.x[idx...]
Base.size(x::TSeries) = size(x.x)
Base.size(x::TSeries,idx...) = size(x.x,idx...)

samplerate(x::TSeries) = x.fs
starttime(x::TSeries) = x.t0
nsamples(x::TSeries) = size(x.x,1)

sampletimes(x::TSeries) = range(starttime(x), step=1/samplerate(x), length=nsamples(x))

