const minORF = 39
const transspliced = ["rps12"]
const fusedparts = ["petB_1", "petD_1", "rpl16_1", "rps12_5"]

MayBeString = Union{Nothing,String}

function build_transspliced_genes!(transsplicedparts, transsplicedgms)
    gmstoadd = Set{Vector{Chloe2.FeatureMatch}}()
    extended = false
    for tgm in transsplicedgms
        o = 1
        while any(x -> only(partorder(x)) == o, tgm)
            o += 1;
        end
        newparts = []
        for gm in transsplicedparts
            if any(x -> only(partorder(x)) == o, gm)
                push!(newparts, gm)
            end
        end
        isempty(newparts) && break
        extended = true
        for part in newparts[2:end]
            push!(gmstoadd, append!(copy(tgm), part)) # copy tgm before adding new part
        end
        append!(tgm, first(newparts))
    end
    if ~extended
        return transsplicedgms;
    end
    union!(transsplicedgms, gmstoadd)
    build_transspliced_genes!(transsplicedparts, transsplicedgms)
end

function fix_exon_borders!(gm::Vector{FeatureMatch}, genome, fstarts, fstartcodons, fstops, rev_genome, rstarts, rstartcodons, rstops)
    isempty(gm) && return gm
    glength = length(genome)
    firstexon = first(gm)
    #fix start & stop codons
    if ~ismissing(firstexon) && firstexon.type == "CDS"
        if firstexon.strand == '+'
            hmm_start_codon = genome[firstexon.target_from:(firstexon.target_from+2)]
            (starts, startcodons, stops) = (fstarts, fstartcodons, fstops)
        else
            hmm_start_codon = rev_genome[firstexon.target_from:(firstexon.target_from+2)]
            (starts, startcodons, stops) = (rstarts, rstartcodons, rstops)
        end
        if hmm_start_codon ∉ [dna"ATG", dna"GTG", dna"ACG"] # assume that if predicted start is ACG it is edited to AUG
            fix_start_codon!(gm, (starts, startcodons, stops), glength)
        end
    end
    #if gene(gm) == "rps4"; println(gm); end
    lastexon = last(gm)
    if ~ismissing(lastexon) && lastexon.type == "CDS"
        if lastexon.strand == '+'
            hmm_stop_codon = genome[(lastexon.target_from+lastexon.target_length):(lastexon.target_from+lastexon.target_length+2)]
            stops = fstops
        else
            hmm_stop_codon = rev_genome[(lastexon.target_from+lastexon.target_length):(lastexon.target_from+lastexon.target_length+2)]
            stops = rstops
        end
        fix_stop_codon!(gm, hmm_stop_codon, stops, glength)
    end
    #if gene(gm) == "rps4"; println(gm); end
    fix_splice_junctions!(gm, glength)
    gm
end

function maybeio(f::Function, io::Union{String,IO}, mode::AbstractString)
    if isa(io, IO)
        f(io)
    else
        open(io, mode) do out
            f(out)
        end
    end
end

function chloeone(tempfile::TempFile, infile::Union{String,IO}, edits::MayBeString; sensitivity = false)
    id, fwd_seq = maybeio(infile, "r") do infa
        record = first(FASTA.Reader(infa))
        identifier(record), FASTA.sequence(LongDNA{4}, record)
    end

    rev_seq = BioSequences.reverse_complement(fwd_seq)
    if ~isnothing(edits)
        glength = length(rev_seq)
        edits_record = readgff(edits)[1]
        for site in edits_record.genes
            nt = first(eachposition(locus(site)))
            if fwd_seq[nt] == DNA_C
                fwd_seq[nt] = DNA_T
            elseif rev_seq[glength-nt+1] == DNA_C
                rev_seq[glength-nt+1] = DNA_T
            end
        end
    end
    chloeone(tempfile, id, fwd_seq, rev_seq; sensitivity = sensitivity)
end

function chloeone(tempfile::TempFile, id::AbstractString, fwd_target::LongDNA{4}, rev_target::LongDNA{4}; sensitivity = false)
    t0 = time()
    glength = length(fwd_target)
    genome = CircularSequence(fwd_target)
    rev_genome = CircularSequence(rev_target)
    @info "$(id): $glength bp"

    #extend genome
    extended_genome = genome[1:(glength+4000)]

    extended_file = tempfilename(tempfile, "$id.extended.fa")
    open(FASTA.Writer, extended_file) do writer
        write(writer, FASTA.Record(id, extended_genome))
    end
    t1 = time()
    @info "$(id): time taken to prepare genome: $(t1 - t0)"

    #find tRNAs
    trn_matches = parse_trn_alignments(search_shattered_genome(tempfile, id, genome; sensitivity = sensitivity), glength)
    filter!(x -> x.target_from <= glength, trn_matches)
    #rationalise_trn_alignments(trn_matches)
    @debug trn_matches
    @info "$(id): found $(length(trn_matches)) tRNA exons"
    t2 = time()
    @info "$(id): time taken to find tRNAs: $(t2 - t1)"
    #println(filter(x -> startswith(only(x.queryparts), "trnI-GAU"), trn_matches))

    #find rRNAs
    #search for rrns using hmmsearch
    rrn_matches = parse_tbl(rrnsearch(extended_file; sensitivity = sensitivity), glength)
    filter!(x -> x.evalue < 1e-10, rrn_matches)
    @debug rrn_matches
    #fix_rrn_ends!(v, ftRNAs, rtRNAs, glength)
    @info "$(id): found $(length(rrn_matches)) rRNA exons"
    t3 = time()
    @info "$(id): time taken to find rRNAs: $(t3 - t2)"

    #find CDSs
    startcodon = ncbi_start_codons[1]
    stopcodon = ncbi_stop_codons[1]
    fstarts, fstartcodons = getcodons(genome, startcodon)
    fstops = codonmatches(genome, stopcodon)
    @debug fstops

    rstarts, rstartcodons = getcodons(rev_genome, startcodon)
    rstops = codonmatches(rev_genome, stopcodon)

    cds_matches = parse_domt(orfsearch(tempfile, id, genome, fstops, rstops, minORF; sensitivity = sensitivity), glength)
    @debug cds_matches
    @info "$(id): found $(length(cds_matches)) CDS exons"
    t4 = time()
    @info "$(id): time taken to find CDSs: $(t4 - t3)"
    #println(filter(x -> gene(x) == "rps4", cds_matches))
    record = GenomicAnnotations.Record{Gene}()
    record.name = id
    record.sequence = genome[1:glength]
    record.circular = true

    templates = CSV.File("$chloe2models/templates.tsv") |> DataFrame
    grouped_templates = groupby(templates, :gene)
    intron_search_time = 0
    for (key, parts) in pairs(grouped_templates)
        primary_model = Vector{FeatureMatch}()
        gene_models = [primary_model]
        #match exons
        for part in eachrow(DataFrames.sort!(parts, :order))
            if part.feature == "CDS"
                addmatches2genemodels!(part, gene_models, cds_matches, glength)
            elseif part.feature == "rRNA"
                addmatches2genemodels!(part, gene_models, rrn_matches, glength)
            elseif part.feature == "tRNA"
                addmatches2genemodels!(part, gene_models, trn_matches, glength)
            end
        end
        #if last(key) == "rps12"; println(gene_models); end
        #match introns
        for part in eachrow(parts)
            if part.feature == "intron"
                t6 = time()
                for gm in gene_models
                    intron = parse_intron_tbl(intronsearch(id, genome, part, gm, tempfile; sensitivity = sensitivity), glength)
                    if ~ismissing(intron)
                        push!(gm, intron)
                    end
                end
                t7 = time()
                intron_search_time += t7 - t6
            end
        end
        #if last(key) == "rps4"; println(gene_models); end
        #finalise gene models
        transsplicedparts = Set{Vector{FeatureMatch}}()
        for gm in gene_models
            isempty(gm) && continue
            sort!(gm; by = x -> (only(partorder(x)), x.target_from))
            for part in eachrow(parts)
                if "$(part.gene)_$(part.order)" ∈ fusedparts
                    fill_missing_exon!("$(part.gene)_$(part.order)", gm)
                    #if gene(gm) == "petB"; println(gm); end
                end
            end
            if last(key) ∈ transspliced
                push!(transsplicedparts, gm)
            else
                fix_exon_borders!(gm, genome, fstarts, fstartcodons, fstops, rev_genome, rstarts, rstartcodons, rstops)
                addgene2record!(record, genome, rev_genome, parts, gm)
            end
        end
        if ~isempty(transsplicedparts)
            transsplicedgms = Set{Vector{FeatureMatch}}()
            push!(transsplicedgms, Vector{FeatureMatch}(undef, 0))
            tgenes = build_transspliced_genes!(transsplicedparts, transsplicedgms)
            for tgene in tgenes
                fix_exon_borders!(tgene, genome, fstarts, fstartcodons, fstops, rev_genome, rstarts, rstartcodons, rstops)
                addgene2record!(record, genome, rev_genome, parts, tgene)
            end
        end
    end
    flag_duplicates!(record)
    t5 = time()
    @info "$(id): time taken to build and verify gene models: $(t5 - t4)"
    @info "$(id): of which searching for introns took: $intron_search_time"
    return record
end

function chloe(
    tempfile::TempFile,
    infile::Union{String,IO};
    edits::MayBeString = nothing,
    outfile_gff::Union{String,IO,Nothing} = nothing,
    sensitivity = false,
    reportpseudos = false
)
    record = chloeone(tempfile, infile, edits; sensitivity = sensitivity)
    t6 = time()

    if ~isnothing(outfile_gff)
        writeGFF(record, outfile_gff; reportpseudos = reportpseudos)
    end

    t7 = time()
    @info "$(record.name): time taken to prepare and write outputs: $(t7 - t6)"
end

function chloe(
    infile::Union{String,IO};
    edits::MayBeString = nothing,
    outfile_gff::Union{String,IO,Nothing} =  nothing=nothing ,
    tempdir::MayBeString = nothing,
    reportpseudos = false,
    sensitivity = false
)
    if tempdir === nothing
        tempdir = "."
    end
    tempfile = TempFile(tempdir)
    try
        chloe(tempfile, infile; edits = edits, outfile_gff = outfile_gff, sensitivity = sensitivity, reportpseudos = reportpseudos)
    finally
        cleanfiles(tempfile)
    end
end
