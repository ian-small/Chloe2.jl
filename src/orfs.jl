#from corresponding NCBI translation table
const ncbi_start_codons = Dict(1 => biore"(ATG)|(GTG)"d)
const ncbi_stop_codons = Dict(1 => biore"(TAG)|(TAA)|(TGA)"d)

const start_codon_penalty = Dict(LongSequence{DNAAlphabet{4}}("GTG") => -3.3692338096657193, LongSequence{DNAAlphabet{4}}("ATG") => 0.0)

const START_CODON_PENALTY_WEIGHTING = 30.0
const INTRUSION_PENALTY_WEIGHTING = 100.0
const EXPANSION_PENALTY_WEIGHTING = 10.0

Codon = LongSubSeq{DNAAlphabet{4}}

function codonmatches(seq::CircularSequence, pattern)::Vector{Vector{Int32}}
    frames = [Int32[] for f in 1:3]
    for m in eachmatch(pattern, seq.sequence[1:seq.length+2])
        n_certain(matched(m)) < 3 && continue
        i::Int32 = m.captured[1]
        push!(frames[mod1(i, 3)], i)
    end
    return frames
end

function getcodons(seq::CircularSequence, pattern)
    positions = [Int32[] for f in 1:3]
    codons = Dict{Int32,Codon}()
    for m in eachmatch(pattern, seq.sequence[1:seq.length+2])
        n_certain(matched(m)) < 3 && continue
        i::Int32 = m.captured[1]
        push!(positions[mod1(i, 3)], i)
        codons[i] = matched(m)
    end
    return positions, codons
end

function writeorfs(writer::FASTA.Writer, id::AbstractString, genome::CircularSequence, strand::Char, stops::Vector{Vector{Int32}}, minORF::Int)
    glength = length(genome)
    for (f, frame) in enumerate(stops)
        nextstop = 0
        for (s, stop) in enumerate(frame)
            nextstop = s < length(frame) ? frame[s + 1] : first(stops[wrapframe(f, glength)]) #frame of next stop when wrapping depends on genome length
            circulardistance(stop+3, nextstop, glength) < minORF && continue
            translation = BioSequences.translate(genome.sequence[range(stop + 3, length=circulardistance(stop+3, nextstop, glength))])
            write(writer, FASTA.Record(id * "*" * strand * "*" * string(f) * "*" * string(stop + 3) * "-" * string(nextstop-1), translation))
        end
    end
end

function orfsearch(tempfile::TempFile, id::AbstractString, genome::CircularSequence, fstops::Vector{Vector{Int32}}, rstops::Vector{Vector{Int32}}, minORF::Int; sensitivity = false)
    out = tempfilename(tempfile, "$id.orfs.fa")
    open(FASTA.Writer, out) do writer
        writeorfs(writer, id, genome, '+', fstops, minORF)
        writeorfs(writer, id, reverse_complement(genome), '-', rstops, minORF)
    end
    hmmsearch = which("hmmsearch")
    results = PipeBuffer()
    hmmpath = joinpath(chloe2models, "cds", "all_cds.hmm")
    cmd = sensitivity ? `$hmmsearch --max -o /dev/null --domtblout /dev/stdout $hmmpath $out` : `$hmmsearch -o /dev/null --domtblout /dev/stdout $hmmpath $out`
    run(cmd, devnull, results, stderr)
    return results
end

function parse_domt(results::IOBuffer, glength::Integer)
    matches = FeatureMatch[]
    for line in readlines(results)
        startswith(line, "#") && continue
        bits = split(line, " ", keepempty=false)
        evalue = parse(Float64, bits[13]) # i-Evalue (indepent Evalue)
        evalue > 1 && continue # filter out poor matches
        orf = bits[1]
        orfbits = split(orf, "*")
        strand = orfbits[2][1]
        seqstart = parse(Int, first(split(last(orfbits), "-")))
        model_from = parse(Int32, bits[16])
        model_to = parse(Int32, bits[17])
        orfalifrom = parse(Int32, bits[20])
        orfalito = parse(Int32, bits[21])
        ali_from = seqstart + 3 * (orfalifrom - 1)
        ali_length = 3 * (orfalito - orfalifrom + 1)
        #ali_from > glength && continue # don't retain matches starting in extension
        # note that model coordinates are converted to nucleotide coordinates
        push!(matches, FeatureMatch(orf, [bits[4]], strand, "CDS", 3 * model_from - 2, 3 * model_to, ali_from, ali_length, evalue))
    end
    close(results)
    #println(filter(x -> startswith(x.query, "ycf2"), matches))
    rationalise_matches!(matches, glength)
end

function fix_start_codon!(gene_model, (starts, startcodons, stops), glength)

    function is_possible_start(start, model_start, leftwindow, rightwindow, glength)
        d = circulardistance(start, model_start, glength)
        if d > glength / 2
            d -= glength
        end
        #must be in frame
        mod(d, 3) ≠ 0 && return false
        d > leftwindow && return false
        d < -rightwindow && return false
        return true
    end

    firstexon = first(gene_model)
    frame = mod1(firstexon.target_from, 3)
    inframe_stops = sort(stops[frame])
    local stop_idx::Int
    while(true)
        stop_idx = searchsortedfirst(inframe_stops, firstexon.target_from)
        downstream_stop = inframe_stops[mod1(stop_idx, length(inframe_stops))]
        if circulardistance(firstexon.target_from, downstream_stop, glength) < firstexon.target_length/2 #stop is early in predicted exon
            firstexon.target_length = circulardistance(downstream_stop + 3, firstexon.target_from + firstexon.target_length - 1, glength)
            firstexon.target_from = downstream_stop + 3
        elseif firstexon.target_length/2 < circulardistance(firstexon.target_from, downstream_stop, glength) <  firstexon.target_length #stop is late in predicted exon
            firstexon.target_length = circulardistance(firstexon.target_from, downstream_stop, glength)
        else
            break
        end
    end
    upstream_stop = inframe_stops[mod1(stop_idx - 1, length(inframe_stops))]
    distance_to_upstream_stop = circulardistance(upstream_stop, firstexon.target_from, glength)
    @debug "$(firstexon.target_from), $distance_to_upstream_stop"
   # if firstexon.query == "ndhK_1"; println("$(firstexon.target_from)\t$distance_to_upstream_stop"); end

    possible_starts = Int32[]
    for ps in starts[frame]
        if is_possible_start(ps, firstexon.target_from, min(200, distance_to_upstream_stop), 100, glength)
            push!(possible_starts, ps)
        end
    end
    #if firstexon.query == "ndhK_1"; println(possible_starts); end
    isempty(possible_starts) && return
    @debug possible_starts

    #apply penalties
    penalties = zeros(Float64, length(possible_starts))
    # start codon penalty
    for (i, ps) in enumerate(possible_starts)
        penalties[i] += START_CODON_PENALTY_WEIGHTING * start_codon_penalty[startcodons[ps]]
    end
    # intrusion/expansion penalty
    for (i, ps) in enumerate(possible_starts)
        relative_to_hmm = circulardistance(firstexon.target_from, ps, glength)
        if relative_to_hmm > glength/2; relative_to_hmm -= glength; end
        if relative_to_hmm > 0
            penalties[i] -= INTRUSION_PENALTY_WEIGHTING * relative_to_hmm
        elseif relative_to_hmm < 0
            penalties[i] += EXPANSION_PENALTY_WEIGHTING * relative_to_hmm
        end
    end
    #if firstexon.query == "ndhK_1"; println(penalties); end
    @debug penalties
    beststart = possible_starts[argmax(penalties)]
    
    #modify FeatureMatch
    firstexon.target_length += firstexon.target_from - beststart
    firstexon.target_from = beststart
end

function fix_stop_codon!(gene_model, stop_codon, stops, glength)
    lastexon = last(gene_model)
    #pick first stop
    orfbits = split(lastexon.target_id, "*")
    frame = parse(Int, orfbits[3])
    stop_idx = searchsortedfirst(stops[frame], lastexon.target_from) #index of first in-frame stop following start of last exon
    if stop_idx == length(stops[frame]) + 1 #no stop before end of genome
        stop_idx = searchsortedfirst(stops[wrapframe(frame, glength)], 1) #restart search from start of genome in the correct frame
    end
    next_stop = stops[frame][stop_idx]
    @debug "first stop: $next_stop"

    #modify FeatureMatch
    if stop_codon ∉ [dna"TAA", dna"TGA", dna"TAG", dna"CAA", dna"CGA", dna"CAG"] || circulardistance(lastexon.target_from, next_stop, glength) < lastexon.target_length
        lastexon.target_length = next_stop - lastexon.target_from
    end
    lastexon.target_length += 3 # to include stop in last exon
end

function fill_missing_exon!(part, gene_model)
    #gene-specific cases due to joint CDS/intron models
    if only(first(gene_model).queryparts) ∈ ["petB_2a", "petD_2a", "rpl16_2a"]
        intron = first(gene_model)
        pushfirst!(gene_model, FeatureMatch("$(intron.target_id)*$(intron.strand)*$(mod1(intron.target_from, 3))*$(string(intron.target_from))-$(string(intron.target_from + model_lengths[part] - 1))",
             [part], intron.strand, "CDS", 1, model_lengths[part], intron.target_from, model_lengths[part], intron.evalue))
        intron.target_from += model_lengths[part]
        intron.target_length -= model_lengths[part]
    elseif part == "rps12_5"
        if only(partorder(last(gene_model))) == 4   #rps12 has second intron
            intron = last(gene_model)
            exonstart = intron.target_from + intron.target_length
            frame = mod1(exonstart + 2, 3) #rps12_6 is always phase 2
            push!(gene_model, FeatureMatch("$(intron.target_id)*$(intron.strand)*$(string(frame))*$(string(exonstart))-$(string(exonstart + 25))", ["rps12_5"], intron.strand, "CDS", 1, 26, exonstart, 26, intron.evalue))
        elseif only(partorder(last(gene_model))) == 3 #rps12 lacks second intron
            exon = last(gene_model)
            push!(gene_model, FeatureMatch(exon.target_id, ["rps12_5"], exon.strand, "CDS", 1, 26, exon.target_from + exon.target_length, 26, exon.evalue))
        end
    end
end

function ntcoords2aacoords(x::Real)
    ceil(x / 3)
end

#calculate new frame when frame wraps the end of the genome
function wrapframe(frame, glength)
    mod1(frame - mod1(glength, 3), 3)
end

#relies on target_id having the form id*strand*frame*start*stop
function frame(f::FeatureMatch)
    @assert f.type == "CDS"
    parse(Int, split(f.target_id, "*")[3])
end

#Check if two CDS Features are in-frame
function in_frame(f1::FeatureMatch, f2::FeatureMatch)
    @assert f1.type == "CDS" && f2.type == "CDS"
    frame(f1) == frame(f2)
end

function splice(gm::Vector{FeatureMatch}, genome, rev_genome)
    mrna = LongDNA{4}()
    for feature in gm
        if feature.type == "CDS"
            start = feature.target_from
            stop = feature.target_from + feature.target_length - 1
            dnaseq = feature.strand == '+' ? genome : rev_genome
            append!(mrna, dnaseq[range(feature.target_from; length = feature.target_length)])
        end
    end
    mrna
end



