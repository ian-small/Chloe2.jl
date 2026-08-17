import Logging

function parse_commandline(args)
    version = pkgversion(Chloe2)
    s = ArgParseSettings(;
        prog = "Chloe2",
        description = "annotates angiosperm plastid genomes",
        usage = "Chloe2/src/command.jl [options] <FASTA_files or directories>",
        epilog = "If there is more than one fasta file to annotate then if the options (--gff etc.) are *not* directories
              they will be used as suffixes for the output filenames and they will be put alongside the input fasta files.",
        version = string(version),
        add_version = true
    )
    #! format: off
    @add_arg_table! s begin
        "--edits"
            help = "file/dir for GFF input containing edit site information"
            arg_type = String
        "--gff"
            help = "file/dir for gff output"
            arg_type = String
            required = true
        "--loglevel"
            help = "loglevel (info,warn,error,debug)"
            arg_type = String
            default = "info"
        "--tempdir"
            help = "directory to write temporary files into (defaults to /tmp or similar...)"
            arg_type = String
        "--overwrite"
            help = "overwrite existing ouput files"
            action = :store_true
        "--fail-early"
            help = "if Chloe fails on multiple inputs then fail immediately"
            action = :store_true
        "--max"
            help = "uses --max setting for cmscan and hmmsearch; slow but sensitive"
            action = :store_true
        "--pseudo"
            help = "reports incomplete or otherwise problematic features as pseudogenes"
            action = :store_true
        "infiles"
            help = "files/directories for fasta input"
            nargs = '+'
    end
    return parse_args(args, s)
end
const LOGLEVELS = Dict("info" => Logging.Info, "debug" => Logging.Debug, "warn" => Logging.Warn,
    "error" => Logging.Error)

function chloe_main(args=ARGS)
    args = parse_commandline(args)
    llevel = get(LOGLEVELS, lowercase(args["loglevel"]), Logging.Warn)
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

    tmpdir = args["tempdir"]
    if tmpdir === nothing
        tmpdir = tempdir()
    end

    if all([isnothing(a) for a in [args["gff"]]])
        @warn "no output specified! type --help"
        return -1
    end
    function readfiles(d, ext)
        if isdir(d)
            return filter(x -> endswith(x, ext), readdir(d, join=true))
        end
        [d]
    end
    fastafiles = [fa for d in args["infiles"] for fa in readfiles(d, r"\.(fa|fna|fasta)")]
    if length(fastafiles) != 1
        ofunc = getout
        if ~isnothing(args["gff"]) && ~isdir(args["gff"])
            @warn("if multiple fasta files are given then --gff must be a directory")
            return -1
        end
    else
        ofunc = getout1
    end
    gfffiles = args["edits"]
    if ~isnothing(gfffiles)
        gfffiles = readfiles(gfffiles, ".gff")
    else
        gfffiles = fill(nothing, length(fastafiles))
    end
    @assert length(gfffiles) == length(fastafiles)
    function doone(fasta, edits; overwrite = false, sensitivity = false, reportpseudos = false)
        ncid = Ref{String}("")
        #try
            accession = first(splitext(basename(fasta)))
            if ~isnothing(edits)
                @assert startswith(basename(edits), accession)
            end
            ncid[] = basename(accession)
            outfile_gff = ofunc(accession, args["gff"], ".gff")
            if ~overwrite && isfile(outfile_gff)
                @warn "$outfile_gff exists and overwrite is false"
                return -1
            end
            @info "$fasta"
            chloe(fasta; edits=edits, outfile_gff=outfile_gff, tempdir=tmpdir, sensitivity = sensitivity, reportpseudos = reportpseudos)
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
    overwrite = args["overwrite"]
    sensitivity = args["max"]
    reportpseudos = args["pseudo"]
    #read model lengths from .hmm and .cm files
    get_model_lengths()
    Base.exit_on_sigint(false)
    if Threads.nthreads() == 1 || length(fastafiles) == 1
        for (fasta, edits) in zip(fastafiles, gfffiles)
            doone(fasta, edits; overwrite = overwrite, sensitivity = sensitivity, reportpseudos = reportpseudos)
        end
    else
        asyncmap(x -> doone(x[1], x[2]; overwrite = overwrite, sensitivity = sensitivity, reportpseudos = reportpseudos), collect(zip(fastafiles, gfffiles)); ntasks = Threads.nthreads())
    end
    0 # return 0 for success
end
