
const LACKS_ESSENTIAL_FEATURE = "lacks an essential part of the gene"
const LACKS_START_CODON = "lacks a start codon"
const PREMATURE_STOP_CODON = "has a premature stop codon"
const CDS_NOT_DIVISIBLE_BY_3 = "CDS is not divisible by 3"
const INFERIOR_COPY = "probable pseudogene as better copy exists in the genome"
const OVERLAPPING_FEATURE = "better-scoring feature overlaps with this one"

#group features into transcription units, i.e. split trans-spliced genes into separate parts
function features2transcriptionunits(features::Vector{FeatureMatch}, glength::Int)
    tus = Vector{Vector{FeatureMatch}}()
    current_tu = FeatureMatch[]
    for f in features
        if isempty(current_tu) || genome_adjacent(last(current_tu), f, glength)
            push!(current_tu, f)
        else
            push!(tus, current_tu)
            current_tu = FeatureMatch[f]
        end
    end
    push!(tus, current_tu)
    tus
end

function tu2genelocus(features::Vector{FeatureMatch}, glength::Int)
    ff = first(features)
    lf = last(features)
    start = ff.strand == '+' ? ff.target_from : reverse_complement(lf.target_from + lf.target_length - 1, glength)
    stop = ff.strand == '+' ?  lf.target_from + lf.target_length - 1 : reverse_complement(ff.target_from, glength)
    locus = ff.strand == '+' ? ClosedSpan(start:stop) : Complement(ClosedSpan(start:stop))
end

function feature2locus(f::FeatureMatch, glength::Int)::AbstractLocus
    start = f.strand == '+' ? f.target_from : reverse_complement(f.target_from + f.target_length - 1, glength)
    stop = f.strand == '+' ?  f.target_from + f.target_length - 1 : reverse_complement(f.target_from, glength)
    ClosedSpan(start:stop) #warning: locus still needs to made Complement if feature is on reverse strand!!
end

function features2exonloci(features::Vector{FeatureMatch}, glength::Int)
    exons = filter(x -> x.type ≠ "intron", features)
    loci::Vector{AbstractLocus} = feature2locus.(exons, glength)
    if length(loci) == 1
        return first(exons).strand == '+' ? first(loci) : Complement(first(loci))
    elseif all(getproperty.(exons, :strand) .== '+')
        return Join(loci)
    elseif all(getproperty.(exons, :strand) .== '-')
        return Complement(Join(loci))
    else
        for (i, f) in enumerate(exons)
            if f.strand == '-'
                loci[i] = Complement(loci[i])
            end
        end
        return Join(loci)
    end
end

function addgene2record!(uid::UUID, record::GenomicAnnotations.Record, genome, rev_genome, parts, features::Vector{FeatureMatch})
    isempty(features) && return
    glength = length(record.sequence)

    #merge adjacent features of same type (e.g. because intron has been lost in this gene)
    merged_features = FeatureMatch[]
    previous_feature = first(features)
    for f in features[2:end]
        if f.strand == previous_feature.strand && f.type == previous_feature.type && (circulardistance(previous_feature.target_from + previous_feature.target_length, f.target_from, glength) < 500) # arbitrary limit below which we think the feature matches are part of the same feature
            previous_feature = merge_matches(previous_feature, f, glength)
        else
            push!(merged_features, previous_feature)
            previous_feature = f
        end
    end
    push!(merged_features, previous_feature)
    features = merged_features
    #if gene(first(features)) == "ycf2"; println(features); end

    #validate gene before adding it to the record
    #println(features)
    problems = String[]
    # check it has all its essential parts
    essential_parts = filter(x -> x.essential == 1, parts)
    essential_part_strings = essential_parts.gene .* "_" .* essential_parts.order
    if any(essential_part_strings .∉ Ref(reduce(vcat, getproperty.(features, :queryparts))))
        #if gene(first(features)) == "clpP1"; println(filter(x -> x.essential == 1, parts).order, "\t", partorder.(features)); end
        push!(problems, LACKS_ESSENTIAL_FEATURE)
    end
    # if CDS, check it is translatable
    if first(features).type == "CDS"
        splicedseq = splice(features, genome, rev_genome)
        if mod(length(splicedseq),3) ≠ 0
            push!(problems, CDS_NOT_DIVISIBLE_BY_3)
        else
            protein = translate(splicedseq)
            if isempty(protein) || first(protein) ∉ [AA_M, AA_V, AA_T] #approximate; intended to allow ATG, GTG or ACG, but is more forgiving...
                push!(problems, LACKS_START_CODON)
            end
            if ~isempty(protein) && last(protein) ∉ [AA_Term, AA_R, AA_Q] #approximate; intended to allow TAA, TAG, TGA or CAA, CAG, CGA but is more forgiving...
                #println(protein)
                push!(problems, PREMATURE_STOP_CODON)
            end
        end
    end
    #if gene(first(features)) == "rps4";println(features); println(problems);end
    #~isempty(problems) && return
    
    #gene
    tus = features2transcriptionunits(features, glength)
    #if startswith(first(features).query, "trnK-UUU"); println(tus); end
    genename = string(gene(first(tus)))
    gene_id = uuid5(uid, genename * string(first(features).target_from))
    gene_locus = length(tus) == 1 ? tu2genelocus(tus[1], glength) : Order(tu2genelocus.(tus, glength))
    featuretype = isempty(problems) ? :gene : :pseudo
    addgene!(record, featuretype, gene_locus; ID = string(gene_id), Name = genename, source = "Chloe2", score = @sprintf("%.2E", minimum(getproperty.(features, :evalue))))

    featuretype == :pseudo && return

    #RNA
    rna_locus = features2exonloci(features, glength)
    rna_id = uuid5(gene_id, genename)
    if first(features).type == "rRNA"
        addgene!(record, :rRNA, rna_locus; ID = string(rna_id), Parent = string(gene_id), Name = genename, source = "Chloe2", score = @sprintf("%.2E", minimum(getproperty.(features, :evalue))))
    elseif first(features).type == "tRNA"
        addgene!(record, :tRNA, rna_locus; ID = string(rna_id), Parent = string(gene_id), Name = genename, source = "Chloe2", score = @sprintf("%.2E", minimum(getproperty.(features, :evalue))))
    end

    #CDS
    if first(features).type == "CDS"
        addgene!(record, :mRNA, rna_locus; ID = string(rna_id), Parent = string(gene_id), Name = genename, source = "Chloe2", score = @sprintf("%.2E", minimum(getproperty.(features, :evalue))))
        cds_id = uuid5(rna_id, genename)
        addgene!(record, :CDS, rna_locus; ID = string(cds_id), Parent = string(rna_id), Name = genename, source = "Chloe2", score = @sprintf("%.2E", minimum(getproperty.(features, :evalue))), phase = 0)
    end
end

function flag_duplicates!(record::GenomicAnnotations.Record)
    #group genes by id
    for gn in unique(record.genedata.Name)
        copies = @genes(record, gene, :Name == gn)
        length(copies) == 1 && continue
        sort!(copies; by = x -> parse(Float64, get(x, :score, 100.0)))
        for c in copies[2:end]
            if parse(Float64, get(c, :score, 100.0)) > 1.05 * parse(Float64, get(copies[1], :score, 100.0))
                feature!(c, :pseudo)
            end
        end
    end
    record
    #all but min evalue set to pseudogene
end