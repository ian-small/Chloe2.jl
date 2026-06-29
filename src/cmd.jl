import Logging

function get_args()
    args = @dictarguments begin
        @helpusage "Chloe2/src/command.jl [options] <FASTA_files or directories>"
        @helpdescription """
            If there is more than one fasta file to annotate
            then if the options (--gff etc.) are *not* directories
            they will be used as suffixes for the output filenames and
            they will be put alongside the input fasta files.
            """
        @argumentflag sensitivity "--max"
        @arghelp "uses --max setting for cmscan and hmmsearch; slow but sensitive"
        @argumentflag reportpseudos "--pseudo"
        @arghelp "reports incomplete or otherwise problematic features as pseudogenes"
        @argumentoptional String GFF_in "--edits"
        @arghelp "file/dir for gff input containing edit site information"
        @argumentoptional String FA_out "--fa"
        @arghelp "file/dir for fasta output"
        @argumentoptional String GFF_out "--gff"
        @arghelp "file/dir for gff output"
        @argumentoptional String GB_out "--tbl"
        @arghelp "file/dir for .tbl output (for GenBank submissions)"
        @argumentdefault String "info" loglevel "--loglevel"
        @arghelp "loglevel (info,warn,error,debug)"
        @argumentoptional String tmpdir "--tempdir"
        @arghelp "directory to write temporary files into (defaults to /tmp or similar...)"
        @argumentflag failearly "--fail-early"
        @arghelp "if Chloe fails on multiple FASTA inputs then fail immediately"
        # @positionalrequired String FASTA_file
        @positionalleftover String FASTA_files "fastafiles"
        # @arghelp "files/directories for fasta input"
    end
    args
end
const LOGLEVELS = Dict("info" => Logging.Info, "debug" => Logging.Debug, "warn" => Logging.Warn,
    "error" => Logging.Error)

function main()
    args = get_args()
    llevel = get(LOGLEVELS, lowercase(args[:loglevel]), Logging.Warn)
    global_logger(ConsoleLogger(stderr, llevel, meta_formatter=Logging.default_metafmt))

    function getout(accession, out, ext)
        function de(ext)
            if !startswith(ext, ".")
                return ".$(ext)"
            end
            ext
        end
        if out === nothing
            return nothing
        end
        if isdir(out)
            return joinpath(out, basename(accession) * ext)
        end
        return accession * de(out)
    end
    function getout1(accession, out, ext)
        if out === nothing
            return nothing
        end
        if isdir(out)
            return joinpath(out, basename(accession) * ext)
        end
        return out
    end

    tmpdir = args[:tmpdir]
    if tmpdir === nothing
        tmpdir = tempdir()
    end

    if all([isnothing(a) for a in [args[:GFF_out], args[:FA_out], args[:GB_out]]])
        println(stderr, "no output specified! type --help")
        return
    end
    function readfiles(d, ext)
        if isdir(d)
            return filter(x -> endswith(x, ext), readdir(d, join=true))
        end
        [d]
    end
    fastafiles = [fa for d in args[:FASTA_files] for fa in readfiles(d, r"\.(fa|fna|fasta)")]
    if length(fastafiles) != 1
        ofunc = getout
    else
        ofunc = getout1
    end
    gfffiles = args[:GFF_in]
    if ~isnothing(gfffiles)
        gfffiles = readfiles(gfffiles, ".gff")
    else
        gfffiles = fill(nothing, length(fastafiles))
    end
    @assert length(gfffiles) == length(fastafiles)
    function doone(fasta, edits; sensitivity = false, reportpseudos = false)
        ncid = Ref{String}("")
        #try
            accession = first(splitext(basename(fasta)))
            if ~isnothing(edits)
                @assert startswith(basename(edits), accession)
            end
            ncid[] = basename(accession)
            outfile_gff = ofunc(accession, args[:GFF_out], ".gff")
            outfile_fa = ofunc(accession, args[:FA_out], ".fa")
            outfile_gb = ofunc(accession, args[:GB_out], ".tbl")
            if isfile(outfile_gff); return; end
            @info "$fasta"
            chloe(fasta; edits=edits, outfile_gff=outfile_gff, outfile_gb=outfile_gb, outfile_fa=outfile_fa, tempdir=tmpdir, sensitivity = sensitivity, reportpseudos = reportpseudos)
        #= catch e
            if e isa InterruptException
                @info "Abort!"
                exit(0)
            end
            @error "\"$(ncid[])\" failed! $(e)"
            if args[:failearly]
                rethrow()
            end
        end =#
    end
    sensitivity = args[:sensitivity]
    reportpseudos = args[:reportpseudos]
    #read model lengths from .hmm and .cm files
    get_model_lengths()
    Base.exit_on_sigint(false)
    if Threads.nthreads() == 1
        for (fasta, edits) in zip(fastafiles, gfffiles)
            doone(fasta, edits; sensitivity = sensitivity, reportpseudos = reportpseudos)
        end
    else
        asyncmap(x -> doone(x[1], x[2]; sensitivity = sensitivity, reportpseudos = reportpseudos), collect(zip(fastafiles, gfffiles)); ntasks = Threads.nthreads())
    end
end


