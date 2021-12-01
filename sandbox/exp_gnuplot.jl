using DelimitedFiles
using NLS_Models # for spectrum

const RegisteredData_UUID = typeof(hash(1))

# convert data id to gnuplot id
#
to_gnuplot_uuid(uuid::RegisteredData_UUID) = "\$G"*string(uuid)


mutable struct GnuplotScript
    _registered_data::Dict{RegisteredData_UUID,Any}
    _script::String
end 

GnuplotScript() = GnuplotScript(Dict{RegisteredData_UUID,AbstractArray}(),String(""))

# check if data has already been registered
#
function is_registered(gp::GnuplotScript,uuid::RegisteredData_UUID)
    haskey(gp._registered_data,uuid)
end

# Register data and return associated data uuid.
#
function register_data!(gp::GnuplotScript,data::AbstractVecOrMat;
                       copy_data::Bool=true)::RegisteredData_UUID

    # already registered
    uuid = hash(data)

    if !is_registered(gp,uuid)
        if copy_data
            data = copy(data)
        end
        gp._registered_data[uuid]=data
    end

    uuid
end

function register_data!(gp::GnuplotScript,data::Spectrum;
                        copy_data::Bool=true)::RegisteredData_UUID
    register_data!(gp,hcat(data.X,data.Y),copy_data=copy_data)
end


function plot!(gp::GnuplotScript,uuid::RegisteredData_UUID,plot_arg::String)
    @assert is_registered(gp,uuid)

    gp._script = gp._script * "plot $(to_gnuplot_uuid(uuid)) " * plot_arg * "\n"

    gp
end

function replot!(gp::GnuplotScript,uuid::RegisteredData_UUID,plot_arg::String)

    gp._script = gp._script * "re"

    plot!(gp,uuid,plot_arg)
end

# ================

function write_data(io::IO,gp::GnuplotScript)
    for (k,d) in gp._registered_data
        println(io,"$(to_gnuplot_uuid(k)) << EOD")
        writedlm(io,d)
        println(io,"EOD")
    end 
end
function write(script_file::String,gp::GnuplotScript)
    io = open(script_file, "w");
    write_data(io,gp)
    print(io,gp._script)
    close(io)
end
    
# ****************************************************************
# DEMO 
# ****************************************************************
gp = GnuplotScript()

id_1 = register_data!(gp,rand(5))
id_2 = register_data!(gp,rand(5,2))

gp = plot!(gp,id_1,"u 1 w l")
gp = replot!(gp,id_2,"u 1:2 w l")

write("demo.gp",gp)
