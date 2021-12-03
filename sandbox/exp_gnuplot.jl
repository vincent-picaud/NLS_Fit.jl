using DelimitedFiles
using NLS_Models # for spectrum

const RegisteredData_UUID = typeof(hash(1))

# convert data id to gnuplot id
#
to_gnuplot_uuid(uuid::RegisteredData_UUID) = "\$G"*string(uuid)


mutable struct GnuplotScript
    _registered_data::Dict{RegisteredData_UUID,Any}
    _script::String
    _any_plot::Bool
end 

GnuplotScript() = GnuplotScript(Dict{RegisteredData_UUID,AbstractArray}(),String(""),false)

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

function free_form(gp::GnuplotScript,gp_line::AbstractString)
    gp._script *= gp_line * "\n"

    gp
end

function plot!(gp::GnuplotScript,uuid::RegisteredData_UUID,plot_arg::String)
    @assert is_registered(gp,uuid)

    gp._script = gp._script * "plot $(to_gnuplot_uuid(uuid)) " * plot_arg * "\n"
    gp._any_plot = true
    gp
end

function replot!(gp::GnuplotScript,uuid::RegisteredData_UUID,plot_arg::String)

    # prevent from using replot for the first plot
    if gp._any_plot
        gp._script = gp._script * "re"
    end
    
    plot!(gp,uuid,plot_arg)
end

# vertical ----------------
#
function add_vertical_line!(gp::GnuplotScript,position::Float64;name::Union{AbstractString,Nothing})
    if name != Nothing
        gp._script  *= "set label at $position, 0.0 '$name' rotate by 90 front left offset -1,1,0 tc ls 1\n"
    end
    gp._script *= "set arrow from $position, graph 0 to $position, graph 1 nohead front\n"

    gp
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
    println(io,gp._script)
    # add a final replot to be sure that everything is plotted
    println(io,"replot")
    close(io)
end
    
# ****************************************************************
# DEMO 
# ****************************************************************
gp = GnuplotScript()

id_1 = register_data!(gp,10*rand(5))
id_2 = register_data!(gp,10*rand(5,2))

gp = plot!(gp,id_1,"u 1 w l")
gp = replot!(gp,id_2,"u 1:2 w l")
gp = add_vertical_line!(gp,5.0,name="toto")
gp = add_vertical_line!(gp,2.0,name="titititito")

write("demo.gp",gp)
