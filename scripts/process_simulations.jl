#================================================================
Evaluate RDFs and exposure distributions from simulations,
then move the simulation files to external hard drive.

Exposure pdfs are evaluated per individual.
RDFs are evaluated by pooling together all the
timepoints and individuals in the simulation.
================================================================#

using DrWatson
@quickactivate "ResidenceTimes"
using ResidenceTimes
using CSV
using DataFrames
using StatsBase
using KernelDensity

function archive!(filename)
    mv(
        datadir("abm", filename),
        joinpath("/media/Elements/ResidenceTimes/data/abm/", filename);
        force=true
    )
    open(datadir("abm_archived.txt"), "a") do io
        println(io, filename)
    end
end

filenames = readdir(datadir("abm"))
L = 1e3 # L is always 1mm
npoints = 100 # points for the rdf sampling
r = range(1, L/2; length=npoints)
for filename in filenames
    prefix, config, suffix = parse_savename(filename)
    df = CSV.read(datadir("abm", filename), DataFrame)
    # distribution of individual exposures
    df_exposure = combine(groupby(df, :id), :c => sum)
    CSV.write(
        datadir("exposure", savename("exposure", config, "csv")),
        df_exposure
    )
    # radial distribution functions
    # each Cs value is compared to corresponding Cs=0 simulation
    iszero(config["Cs"]) && continue
    config_random = copy(config)
    config_random["Cs"] = 0.0
    filename_random = savename(prefix, config_random, suffix)
    !isfile(datadir("abm", filename_random)) && continue
    df_random = CSV.read(datadir("abm", filename_random), DataFrame)
    Pr = kde(df.r; bandwidth=25)
    Pr0 = kde(df_random.r; bandwidth=25)
    k = pdf(Pr, r)
    k0 = pdf(Pr0, r)
    g = k ./ k0
    df_rdf = DataFrame(; r, g)
    CSV.write(
        datadir("rdf", savename("rdf", config, "csv")),
        df_rdf
    )
    # move abm data to hard drive
    archive!(filename)
end
# the Cs=0 files have been kept, remove them now
for filename in filenames
    !isfile(datadir("abm", filename)) && continue
    archive!(filename)
end
