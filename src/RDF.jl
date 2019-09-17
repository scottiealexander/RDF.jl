module RDF

export RDFFile, read_channel, spike_times, sample_times, RDFParams, trial_labels

const SAMPLE_RATE = 24000.0f0

# ============================================================================ #
mutable struct RDFParams
    d::Dict{String, Any}
end
# ---------------------------------------------------------------------------- #
RDFParams(ifile::AbstractString) = RDFParams(parse_ini_file(ifile))
# ---------------------------------------------------------------------------- #
function trial_labels(p::RDFParams, x::Vector{UInt8})
    nv = length(p.d["values"])
    lab = Vector{Float64}()
    for k in x
        if k <= nv
            push!(lab, p.d["values"][k])
        end
    end
    return lab
end
# ---------------------------------------------------------------------------- #
@inline Base.getindex(p::RDFParams, k::AbstractString) = p.d[k]
# ============================================================================ #
mutable struct RDFFile
    path::String
    ids::Vector{UInt8}
    types::Vector{UInt8}
    type_sizes::Dict{UInt8,Int}
    npacket::Int
    first_packet::Int
    first_packet_time::Float32
end
# ---------------------------------------------------------------------------- #
function RDFFile(path::AbstractString)
    r = open(path, "r") do io

        nchan = read(io, UInt8)
        ids = Vector{UInt8}(undef, nchan)
        read!(io, ids)
        types = Vector{UInt8}(undef, nchan)
        read!(io, types)

        np = read(io, Int32)

        offset = position(io)

        seek(io, offset + sizeof(UInt8))
        first_packet_time = read(io, Float32)

        type_sizes = Dict{UInt8,Int}(x => 0 for x in ids)

        return RDFFile(path, ids, types, type_sizes, np, offset,
            first_packet_time)
    end

    for k in r.ids
        r.type_sizes[k] = sizeof(get_eltype(r, k))
    end

    return r
end
# ---------------------------------------------------------------------------- #
function read_channel(r::RDFFile, id::Integer)
    0 < id < 256 || error("Invalid channel: must be 0 < id < 256")

    elt = get_eltype(r, id)
    data = Vector{elt}(undef, 0)
    ts = Vector{Float32}(undef, 0)
    idx = Vector{Int}(undef, 0)
    open(r.path, "r") do io
        seekend(io)
        len = position(io)
        seek(io, r.first_packet)

        for k in 1:r.npacket
            idp = read(io, UInt8)

            # correct packet arrival times to be relative to the first packet
            # of the file, this way continuous channels and event channels
            # will be synchronized "automatically"
            time = read(io, Float32) - r.first_packet_time
            nel = read(io, Int32)

            if idp == id
                push!(idx, length(data))
                tmp = Vector{elt}(undef, nel)
                read!(io, tmp)
                append!(data, tmp)
                push!(ts, time)
            else
                pos = position(io) + nel * r.type_sizes[idp]
                if 0 < pos <= len
                    seek(io, pos)
                end
            end
        end
    end
    return data, ts, idx
end
# ---------------------------------------------------------------------------- #
function get_eltype(r::RDFFile, id::Integer)
    kchan = findfirst(isequal(id), r.ids)

    kchan == nothing && error("Invalid id: $(id)")

    typ = r.types[kchan]

    if typ & 0x01 == 0x00
        if typ == 0x02
            elt = UInt8
        elseif typ == 0x04
            elt = UInt16
        elseif typ == 0x06
            elt = UInt32
        else
            error("Invalid type: $(typ)")
        end
    else
        if typ == 0x03
            elt = Int8
        elseif typ == 0x05
            elt = Int16
        elseif typ == 0x07
            elt = Int32
        elseif typ == 0x0f
            elt = Float32
        else
            error("Invalid type: $(typ)")
        end
    end

    return elt
end
# ---------------------------------------------------------------------------- #
function spike_times(data::AbstractVector{<:Real})
    win = 64
    thr = 0.97f0
    ts = Vector{Float32}()
    k = 1
    si = (1.0f0 / SAMPLE_RATE)
    while k < length(data)
        if data[k] > thr
            ke = min(k+win-1, length(data))
            _, kmx = findmax(data[k:ke])
            if kmx != nothing
                kmx += k - 1
                push!(ts, (kmx - 1) * si)
                k += win
            else
                k = length(data) + 1
            end
        else
            k += 1
        end
    end
    return ts
end
# ---------------------------------------------------------------------------- #
@inline sample_times(data::AbstractVector{<:Real}) =
    range(0, length(data) * (1.0 / SAMPLE_RATE), length=length(data))
# ---------------------------------------------------------------------------- #
function parse_ini_file(pth::AbstractString)
    isfile(pth) || error("Cannot ind file \"$(pth)\"")
    d = Dict{String, Any}()
    for line in eachline(pth, keep=false)
        m = match(r"([^\=\s]+)\s*\=\s*([^\n]+)", line)
        if m != nothing
            if m[2][1] == '['
                ke = findfirst(isequal(']'), m[2])
                ke = ke == nothing ? length(m[2]) : ke - 1
                val = map(x -> parse(Float64, x), split(m[2][2:ke], r"[,\s]+"))
                d[m[1]] = val
            else
                d[m[1]] = m[2]
            end
        else
            @warn("Invalid statement: \"$(line)\"")
        end
    end
    return d
end
# ============================================================================ #
end
