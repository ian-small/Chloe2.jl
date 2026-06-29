function cmsearch_intron(target::String, intron::String; sensitivity = false)  
    cmsearch = which("cmsearch")
    cmpath = joinpath(chloe2models, "introns", intron * ".cm")
    results = PipeBuffer()
    cmd = sensitivity ? `$cmsearch --max -o /dev/null --tblout /dev/stdout --toponly --notrunc --noali -Z 0.3 $cmpath $target` : `$cmsearch -o /dev/null --tblout /dev/stdout --toponly --hmmonly --notrunc --noali -Z 0.3 $cmpath $target`
    run(cmd, devnull, results, stderr)
    results
end

function intronsearch(id::AbstractString, genome::CircularSequence, part, gene_model, tempfile::TempFile; sensitivity = false)
    try
        leeway = 20 # arbitrary grey zone to account for slop in exon placement
        #defaults
        intron_strand = '+'
        donor_site = 1
        acceptor_site = length(genome)
        partordernumber = parse(Int, part.order[1])
        if endswith(part.order, "a") #first half of intron
            donor_idx = findlast(x -> only(partorder(x)) == partordernumber-1, gene_model)
            if ~isnothing(donor_idx)
                donor_match = gene_model[donor_idx]
                donor_site = donor_match.target_from + donor_match.target_length - leeway
                acceptor_site = donor_site + 300 + leeway
                intron_strand = donor_match.strand
            else #no preceding exon, e.g. case for petB, petD, rpl16
                acceptor_idx = findfirst(x -> only(partorder(x)) == partordernumber+1, gene_model)
                if ~isnothing(acceptor_idx)
                    acceptor_match = gene_model[acceptor_idx]
                    acceptor_site = acceptor_match.target_from
                    donor_site = acceptor_site - 2000 - leeway
                    intron_strand = acceptor_match.strand
                end
            end
        else #second half of intron
            acceptor_idx = findfirst(x -> only(partorder(x)) == partordernumber+1, gene_model)
            if ~isnothing(acceptor_idx)
                acceptor_match = gene_model[acceptor_idx]
                acceptor_site = acceptor_match.target_from + leeway
                donor_site = acceptor_site - 300 - leeway
                intron_strand = acceptor_match.strand
            else #no following exon, e.g. case for rps12_4
                donor_idx = findlast(x -> only(partorder(x)) == partordernumber-1, gene_model)
                if ~isnothing(donor_idx)
                    donor_match = gene_model[donor_idx]
                    donor_site = donor_match.target_from + donor_match.target_length
                    acceptor_site = donor_site + 2000 + leeway
                    intron_strand = donor_match.strand
                end
            end
        end

        #don't search if search space is entire genome
        if donor_site == 1 && acceptor_site == length(genome); return missing; end

        search_seq = intron_strand == '+' ? genome[donor_site:acceptor_site] : reverse_complement(genome)[donor_site:acceptor_site]
        intron_name = "$(part.gene)_$(string(part.order))"
        fname = tempfilename(tempfile, "$id.$intron_name.$intron_strand.$(string(donor_site))-$(string(acceptor_site)).fa")
        #println("searching for $intron_name in $intron_strand.$(string(donor_site))-$(string(acceptor_site))")
        open(FASTA.Writer, fname) do writer
            write(writer, FASTA.Record("$id.$intron_name.$intron_strand.$(string(donor_site))-$(string(acceptor_site))", search_seq))
        end
        return cmsearch_intron(fname, "$intron_name"; sensitivity = sensitivity)
    catch
        println("failed intron search: $id $(part.gene)")
        return missing
    end
end

function parse_intron_tbl(results::Union{Missing, IOBuffer}, glength::Int)
    ismissing(results) && return missing
    intron_matches = FeatureMatch[]
    hits = filter!(x -> ~startswith(x, "#"), readlines(results))
    isempty(hits) && return missing
    for hit in hits
        bits = split(hit, " ", keepempty=false)
        query = bits[3]
        m = match(r"(^[^.]+)(?:\.[0-9]+)?\.([^.]+)\.([+|-])\.(-?[0-9]+)-[0-9]+", bits[1]) # (?:\.[0-9]+)? is optional match to .1 version number on accession
        if isnothing(m)
            println("intronseqname: ", bits[1])
            println(hits)
        end
        id = m.captures[1]
        strand = first(m.captures[3])
        seqstart = parse(Int, m.captures[4])
        evalue = parse(Float64, bits[16])
        model_from = parse(Int, bits[6])
        model_to = parse(Int, bits[7])
        seq_from = parse(Int, bits[8])
        seq_to = parse(Int, bits[9])
        if query ∈ ["clpP1_4b", "rps12_4b"] # account for inset models
            model_to -= 20
            seq_to -= 20
        end
        target_from = seqstart + seq_from - 1
        target_to = seqstart + seq_to - 1
        target_length = target_to - target_from + 1
        # note that target coordinates are in stranded nucleotide coordinates
        push!(intron_matches, FeatureMatch(id, [query], strand, "intron", model_from, model_to, target_from, target_length, evalue))
    end
    close(results)
    @debug intronmatch
    first(rationalise_matches!(intron_matches, glength))
end
