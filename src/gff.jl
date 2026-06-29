function printgff(io::IO, chr::GenomicAnnotations.Record, blacklist::Vector{UInt}; reportpseudos = false)
    iobuffer = IOBuffer()
    ### Header
    if occursin(r"^##gff-version 3", chr.header)
        print(iobuffer, chr.header)
    else
        println(iobuffer, "##gff-version 3")
    end
    ### Body
    for gene in chr.genes
        if index(gene) ∈ blacklist
            feature(gene) ≠ :pseudo && continue # never report non-gene pseudo features
            ~reportpseudos && continue # only report gene pseudo features if reportpseudos is requested
        end
        print(iobuffer, GenomicAnnotations.GFF.gffstring(gene))
    end
    ### Footer
    if ~isempty(chr.sequence)
        println(iobuffer, "##FASTA")
        println(iobuffer, ">", chr.name)
        for s in Iterators.partition(chr.sequence, 80)
            println(iobuffer, join(s))
        end
    end
    print(io, String(take!(iobuffer)))
end

function writeGFF(record::GenomicAnnotations.Record, outfile::String; reportpseudos = false)
    gffrecord = copy(record)
    gffrecord.sequence = dna""
    addgene!(gffrecord, :region, ClosedSpan(1:length(record.sequence));
            Name = record.name,
            ID = record.name,
            Is_circular = record.circular)
    #catalogue pseudogenes
    pseudogenes = filter(x -> feature(x) == :pseudo, gffrecord.genes)
    pseudoids = get.(pseudogenes, :ID, "")
    for g in gffrecord.genes
        if get(g, :Parent, "") ∈ pseudoids
            push!(pseudoids, get(g, :ID, ""))
            push!(pseudogenes, g)
        end
    end
    blacklist = sort!(unique!(index.(pseudogenes)))
    sort!(gffrecord.genes)
    open(GFF.Writer, outfile) do out
        printgff(out.output, gffrecord, blacklist; reportpseudos = reportpseudos)
    end
end
