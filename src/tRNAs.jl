#safer to construct this dynamically fom the model file...
const model2anticodonpos = Dict("trnA-UGC_1" => 34, "trnC-GCA_1" => 33, "trnD-GUC_1" => 35, "trnE-UUC_1" => 35, "trnF-GAA_1" => 34, "trnG-GCC_1" => 33, "trnG-UCC_2" => 10, "trnH-GUG_1" => 36,
     "trnI-CAU_1" => 34, "trnI-GAU_1" => 33, "trnK-UUU_1" => 33, "trnL-CAA_1" => 34, "trnL-UAG_1" => 35, "trnM-CAU_1" => 34, "trnN-GUU_1" => 33, "trnP-UGG_1" => 35,
      "trnQ-UUG_1" => 33, "trnR-ACG_1" => 35, "trnR-UCU_1" => 34, "trnS-GCU_1" => 35, "trnS-GGA_1" => 35, "trnS-UGA_1" => 35, "trnT-GGU_1" => 33, "trnT-UGU_1" => 34, "trnV-GAC_1" => 33,
       "trnV-UAC_1" => 35, "trnW-CCA_1" => 35, "trnY-GUA_1" => 35, "trnfM-CAU_1" => 35)

function anticodon(qfrom::Int, qseq::AbstractString, tseq::AbstractString, pos::Int)
    @assert length(qseq) == length(tseq)
    pos < qfrom && return missing
    pointer = qfrom - 1
    trunc = false
    for (i, c) in enumerate(qseq)
        c == '.' && continue
        c == '*' && continue
        c == ' ' && continue
        isdigit(c) && continue
        if c == ']'
            trunc = false
            continue
        end
        if trunc == true
            continue
        elseif c == '['
            pointer += parse(Int, qseq[i+1:(findnext("]", qseq, i + 2)[1]-1)])
        elseif c == '<'
            trunc = true
            continue
        else
            pointer += 1
        end
        pos < pointer && return missing
        pointer == pos && i+2 <= length(tseq) && return tseq[i:i+2]
    end
    return missing
end

function parse_trn_alignments(results::IOBuffer, glength::Integer)
    trns = FeatureMatch[]
    while !eof(results)
        line = readline(results)
        !startswith(line, ">>") && continue
        readline(results) #header line
        readline(results) #dashes
        bits = split(readline(results), " ", keepempty=false)
        target_from = parse(Int, bits[10])
        tto = parse(Int, bits[11])
        qfrom = parse(Int, bits[7])
        readline(results) #blank
        readline(results) #NC
        readline(results) #CS
        qseqline = strip(readline(results))
        query = qseqline[1:(findfirst(" ", qseqline)[1]-1)]
        qseqline = lstrip(qseqline[length(query)+1:end])
        qseq = qseqline[(findfirst(" ", qseqline)[1]+1):(findlast(" ", qseqline)[1]-1)]
        readline(results) #matches
        tseqline = strip(readline(results))
        target = tseqline[1:(findfirst(" ", tseqline)[1]-1)]
        tseqline = lstrip(tseqline[length(target)+1:end])
        tseq = tseqline[(findfirst(" ", tseqline)[1]+1):(findlast(" ", tseqline)[1]-1)]

        #correct for shattering
        seqrange = match(r"([\+|-])([0-9]+)-([0-9]+)", target)
        tstrand = seqrange.captures[1][1]
        seqstart = parse(Int, seqrange.captures[2])
        target_from += seqstart - 1
        tto += seqstart - 1 

        #anticodon
        if haskey(model2anticodonpos, query)
            expected = only(match(r"^trnf?[A-Z]-([A-Z][A-Z][A-Z])", query).captures)
            anticod = anticodon(qfrom, qseq, tseq, model2anticodonpos[query])
            ismissing(anticod) || anticod ≠ expected && continue
        end
        push!(trns, FeatureMatch(target, [query], tstrand, "tRNA", qfrom, parse(Int, bits[8]), target_from, tto - target_from + 1, parse(Float64, bits[3])))
    end
    close(results)
    rationalise_matches!(trns, glength)
end

function cmsearch_trns(target::String; sensitivity = false)
    cmpath = joinpath(chloe2models, "trns", "all_trns.cm")
    cmscan = which("cmscan")
    results = PipeBuffer()
    cmd = sensitivity ? `$cmscan --max --toponly $cmpath $target` : `$cmscan --toponly $cmpath $target`
    run(cmd, devnull, results, stderr)
    return results
end

function search_shattered_genome(tempfile::TempFile, id::AbstractString, genome::CircularSequence; sensitivity = false)
    shattered_file = tempfilename(tempfile,"$id.shattered.fa")
    revgenome = reverse_complement(genome)
    open(FASTA.Writer, shattered_file) do writer
        overlap = 100
        fragment_size = 2000
        num_fragments = ceil(Int, (length(genome) + overlap)/(fragment_size - overlap))
        for i in 1:num_fragments
            start = (i - 1) * (fragment_size - overlap) + 1
            stop = start + fragment_size - 1
            write(writer, FASTA.Record(id * "+" * string(start) * "-" * string(stop), genome.sequence[start:stop]))
            write(writer, FASTA.Record(id * "-" * string(start) * "-" * string(stop), revgenome.sequence[start:stop]))
        end
    end
    cmsearch_trns(shattered_file; sensitivity = sensitivity)
end