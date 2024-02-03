include("load.jl")
include("./parameters_phi.jl")

for i in length(enumerate(ps))
    mycommand = `sbatch run_snowden.sh $i`
    run(mycommand)
end
