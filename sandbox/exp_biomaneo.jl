using Revise
using NLS_Fit
using NLS_Solver
using NLS_Models
using DelimitedFiles
include("exp_gnuplot.jl")

# ================================================================

import NLS_Models: create_model

# Used to tag isotopic motif within group
#
struct EmbeddedData_IsotopicMotif_Model
    _group_idx::Int
    _local_idx::Int
end

# Used to tag a complete group model
#
# This model is the sum of Isotopic Model contained in the considered group
#
struct EmbeddedData_Group_Model
    _group_idx::Int
end

function NLS_Models.create_model(grouped::GroupedBySupport{IsotopicMotif},idx_group::Int)

    model = Model2Fit_Empty()

    n_isotopicmotif = objects_in_group_size(grouped,idx_group)
    for idx_isotopicmotif in 1:n_isotopicmotif
        isotopic_motif = get_object(grouped,idx_group,idx_isotopicmotif)

        
        model += Model2Fit_TaggedModel(Peak_Motif(Gaussian_Peak(),get_profile_matrix(isotopic_motif)),
                                       EmbeddedData_IsotopicMotif_Model(idx_group,idx_isotopicmotif))
    end

    model = Model2Fit_TaggedModel(model,EmbeddedData_Group_Model(idx_group))
   
    model
end

# Create a vector of models. Each model associated to one group
#
function create_vector_of_models(grouped::GroupedBySupport{IsotopicMotif})
    n_group = group_size(grouped)
    map(idx_group -> create_model(grouped,idx_group), 1:n_group)
end 

# ================================================================

# for X int vcat(collect.(ROI_range)...)

#
# Create a new spectrum restricted to ROI.
#
# In the same time return the length of each ROIs
#
# Attention: use functions from NLS_Models
#
function extract_ROI_spectrum_range_pair(grouped::GroupedBySupport{IsotopicMotif},sp::Spectrum)
    ROI_interval = group_support_as_interval(grouped)
    ROI_range = create_range_from_interval(ROI_interval,spectrum.X)

    Spectrum(extract_ROIs(ROI_range,sp.X),
             extract_ROIs(ROI_range,sp.Y)),
    length.(ROI_range)
end 

# ================================================================

# create stacked model and associated ROI spectrum
# 
function create_stacked_model_ROI_spectrum_pair(grouped::GroupedBySupport{IsotopicMotif},sp_before_ROI::Spectrum)
    vector_of_models =  create_vector_of_models(grouped)
    roi_spectrum, roi_lengths = extract_ROI_spectrum_range_pair(grouped,sp_before_ROI)

    Model2Fit_Stacked(vector_of_models,roi_lengths),roi_spectrum
end

# ****************************************************************

# Inputs
# ================

# not ok
#raw_spectrum = read_spectrum_Biomaneo("/home/picaud/Data/Spectres_Biomaneo/txt-minos NMOG/0000100787.txt")
#raw_spectrum = read_spectrum_Biomaneo("/home/picaud/Data/Spectres_Biomaneo/txt-minos NMOG/0000128803(d5).txt")


# raw_spectrum = read_spectrum_Biomaneo("/home/picaud/Data/Spectres_Biomaneo/Spectres_Biomaneo_MF/Heterozygote HbE/0000000036_digt_MF.txt")
 raw_spectrum = read_spectrum_Biomaneo("/home/picaud/Data/Spectres_Biomaneo/January_2020_normalized/Heterozygote HbE B Thal/0000000017_digt_0001_J4_(Manual)_19-12-20_14-19_0001.txt")
# raw_spectrum = read_spectrum_Biomaneo("/home/picaud/GitHub/NLS_Models.jl/data/0000000095.txt")
# raw_spectrum = read_spectrum_Biomaneo("/home/picaud/GitHub/NLS_Models.jl/data/0000000001.txt")
#raw_spectrum = read_spectrum_Biomaneo("/home/picaud/GitHub/NLS_Models.jl/data/spectrum.txt")

raw_spectrum.Y ./= maximum(raw_spectrum.Y)
vect_of_isotopicmotif = hardcoded_IsotopicMotifVect()

# Remove baseline
# ================
Y_baseline = compute_baseline_snip(raw_spectrum,
                                   snip_halfwindow = 50,
                                   smoothing_halfwindow = 20)

spectrum = raw_spectrum - Y_baseline

begin
    gp = GnuplotScript()
    id = register_data!(gp,hcat(raw_spectrum.X,raw_spectrum.Y,Y_baseline.Y))
    plot!(gp,id,"u 1:2 w l t 'raw'")
    replot!(gp,id,"u 1:3 w l t 'baseline'")
    write("baseline.gp",gp)
end

# Create ROI model & spectrum
# ================
grouped = groupbysupport(vect_of_isotopicmotif,by=get_ROI_interval)
stacked_models,ROI_spectrum = create_stacked_model_ROI_spectrum_pair(grouped,spectrum)

# Add a σ law
# ================

# The created model is a succession of isotopic model, θ=(h,σ)
# now we put an affine model for σ. For that we need:
# 1/ σ_index in θ
# 2/ position m/z of each isopoic model
#
σ_index = collect(2:2:NLS_Fit.parameter_size(stacked_models))
isotopicmotif_centers = map(get_position,grouped.objects)

map_mz_to_σ = Map_Affine_Monotonic(1000.0  => 2.0, 3000.0 => 6.0) # the σ map: m/z -> σ(m/z)

stacked_models_σ_law = Model2Fit_Mapped_Parameters(stacked_models,map_mz_to_σ,σ_index,isotopicmotif_centers)

# Add affine calibration
# ================
recalibration_map = Map_Affine_Monotonic(ROI_spectrum.X[1],ROI_spectrum.X[end])
stacked_models_σ_law_recalibration = Model2Fit_Recalibration(stacked_models_σ_law,recalibration_map)

# Initialize θ
# ================

# create θ : h1,...h21, slaw1,slaw2, cal1,cal2
#
n_θ = NLS_Fit.parameter_size(stacked_models_σ_law_recalibration)
θ_init = ones(n_θ)
θ_lb = zeros(n_θ)
θ_ub = zeros(n_θ) .+ 6
θ_init = solve_linear_parameters(stacked_models_σ_law_recalibration,
                                 ROI_spectrum.X,
                                 ROI_spectrum.Y,
                                 θ_init,
                                 [1:21;])

Y_fit_init = eval_y(stacked_models_σ_law_recalibration,ROI_spectrum.X,θ_init) # for gnuplot

# Solve the problem
# ================
bc = BoundConstraints(θ_lb,θ_ub)
nls = NLS_ForwardDiff_From_Model2Fit(stacked_models_σ_law_recalibration,ROI_spectrum.X,ROI_spectrum.Y)
conf = Levenberg_Marquardt_BC_Conf()

result = NLS_Solver.solve(nls,θ_init,bc,conf)

# Plot solution
# ================
Y_fit = eval_y(stacked_models_σ_law_recalibration,ROI_spectrum.X,solution(result))

recalibrated_spectrum = Spectrum(eval_calibrated_x(stacked_models_σ_law_recalibration,spectrum.X,solution(result)),spectrum.Y)

writedlm("poub_calibration.txt",hcat(spectrum.X,recalibrated_spectrum.X,spectrum.Y))

# Here prepare for local fittings
# ****************************************************************
# TODO
# Compute the calibrated spectrum and keep the model without the
# recalibration extra-layer
# -> must implement visit before 
#recalibrated_spectrum = Spectrum(eval_calibrated_x,

# save text file, to be used by gnuplot
#
# plot "poub.txt" u 1:2 w l
# replot "poub.txt" u 1:4 w l t "init"
# replot "poub.txt" u 1:3 w l lw 2 t "finited"
#
writedlm("poub.txt",hcat(ROI_spectrum.X,ROI_spectrum.Y,Y_fit,Y_fit_init))

# ****************************************************************
solution(result)


NLS_Fit.visit_debug(stacked_models_σ_law_recalibration,ROI_spectrum.Y,ROI_spectrum.X,solution(result))

# This is struct only make sense when used inside a GlobalFitResult struct
# (by example we need grouped or calibrated spectrum to extract useful information)
#
Base.@kwdef struct LocalFit
    data::EmbeddedData_Group_Model
    model::NLS_Fit.Abstract_Model2Fit
    ROI_calibrated_spectrum::Spectrum # CAVEAT: do not replace by
                                      # range, as this allows us to
                                      # use local recalibration
    θ::AbstractVector
end

# Allows model modification, see:
#     add_calibration_shift()
# for instance.
function update_model(lf::LocalFit,new_model::NLS_Fit.Abstract_Model2Fit,new_θ::AbstractVector)
    @assert NLS_Fit.parameter_size(new_model) == length(new_θ)
    
    LocalFit(data = lf.data,
             model = new_model,
             ROI_calibrated_spectrum = lf.ROI_calibrated_spectrum,
             θ = new_θ)
end 

struct GlobalFitResult
    grouped::GroupedBySupport{IsotopicMotif}
    calibrated_spectrum::Spectrum
    local_fit::Vector{LocalFit}
end

import Base: length
Base.length(gfr::GlobalFitResult) = length(gfr.local_fit)

function get_spectrum(gfr::GlobalFitResult,idx::Int)
    gfr.local_fit[idx].ROI_calibrated_spectrum
end 

function get_model(gfr::GlobalFitResult,idx::Int)
    gfr.local_fit[idx].model 
end 
function get_model_parameter(gfr::GlobalFitResult,idx::Int)
    gfr.local_fit[idx].θ
end 

function get_group_idx(gfr::GlobalFitResult,idx::Int)
    # often idx = group_idx but this two indices are different (by
    # example if there are fewer ROI model than group)
    #
    gfr.local_fit[idx].data._group_idx
end

# Give the number of isotopic motifs in ROI idx
function get_isotopicmotif_count(gfr::GlobalFitResult,idx::Int)

    group_idx = get_group_idx(gfr,idx)

    objects_in_group_size(gfr.grouped,group_idx)
end 

# Return isotopicmotif name
function get_isotopicmotif_name(gfr::GlobalFitResult,idx::Int,local_idx::Int)

    group_idx = get_group_idx(gfr,idx)

    object = get_object(gfr.grouped,group_idx,local_idx)

    get_name(object)
end

# Return isotopicmotif position
function get_isotopicmotif_position(gfr::GlobalFitResult,idx::Int,local_idx::Int)

    group_idx = get_group_idx(gfr,idx)

    object = get_object(gfr.grouped,group_idx,local_idx)

    get_position(object)
end

# Extract fit result, grouping model per group and providing the right X,Y (after calibration)
#
# caveat: vector index has no reason to be the group index. The right
# way to get the group idx is to use data field.
#
# Model associated to each group is detected thanks to
#
#   Model2Fit_TaggedModel{ _ , EmbeddedData_Group_Model}
#
# type checking
#
function extract_fit_result_per_group_helper(grouped::GroupedBySupport{IsotopicMotif}, # not necessary Isotopic... as we only use interval
                                             global_model::NLS_Fit.Abstract_Model2Fit,
                                             ROI_spectrum::Spectrum,
                                             θ_fit::AbstractVector,
                                             calibrated_spectrum::Spectrum)
    
    collected_model_per_group = Vector{LocalFit}(undef,0)
    
    visit(global_model, ROI_spectrum.Y,ROI_spectrum.X, θ_fit) do model,Y,X,θ
        # filter model
        #
        if get_tagged_data_type(model) === EmbeddedData_Group_Model
            # process local model 
            data  = get_tagged_data(model)
            model = get_tagged_model(model)
            
            # extract the right ROIs: compared to our
            # initial ROI_spectrum, the recalibration
            # procedure may have changed ROI range when
            # computed from calibrates_spectrum
            group_idx = data._group_idx
            group_interval = group_support_as_interval(grouped,group_idx)
            group_range =  create_range_from_interval(group_interval,calibrated_spectrum.X)
            ROI_calibrated_spectrum = calibrated_spectrum[group_range] # TODO: implement spectrum view
            
            # store the result
            #
            push!(collected_model_per_group,
                  LocalFit(data = data,
                           model = model,
                           ROI_calibrated_spectrum = ROI_calibrated_spectrum,
                           θ=θ))
            return false
        end
        
        true
    end
    
    collected_model_per_group
end

# A wrapper that perform calibration before calling
#
#    extract_fit_result_per_group_calibrated(...)
#
function extract_fit_result_per_group(grouped::GroupedBySupport{IsotopicMotif}, # not necessary Isotopic... as we only use interval
                                      global_model::Model2Fit_Recalibration,
                                      ROI_spectrum::Spectrum,
                                      θ_fit::AbstractVector,
                                      uncalibrated_spectrum::Spectrum)

    # Perform calibration and remove calbibration model
    #
    calibrated_X = eval_calibrated_x(global_model,uncalibrated_spectrum.X,θ_fit)
    global_model_after_calibration = get_calibrated_model(global_model)
    global_model_after_calibration_θ_fit = get_calibrated_model_θ(global_model,θ_fit)
    
    calibrated_spectrum = Spectrum(calibrated_X,spectrum.Y)

    # Continue to process each group model
    #
    local_fit_vect = extract_fit_result_per_group_helper(grouped,
                                                         global_model_after_calibration,
                                                         ROI_spectrum,
                                                         global_model_after_calibration_θ_fit,
                                                         calibrated_spectrum)

    
    GlobalFitResult(grouped,calibrated_spectrum,local_fit_vect)
end 


global_fit_result = extract_fit_result_per_group(grouped,
                                         stacked_models_σ_law_recalibration,
                                         ROI_spectrum,
                                         solution(result),
                                         spectrum)

# include("exp_gnuplot.jl")

function plot_fit(gp_output_file::String,global_fit_result::GlobalFitResult)
    gp = GnuplotScript()

    # set title
    #
    free_form(gp,"set title '$gp_output_file'")

    # plot shaded ROI
    #
    n_ROI = length(global_fit_result)
    for idx_ROI in 1:n_ROI
        ROI_spectrum_uuid = register_data!(gp,get_spectrum(global_fit_result,idx_ROI))

        if idx_ROI==1
            plot!(gp,ROI_spectrum_uuid,"with filledcurve y1=-0.1 lc rgb 'gray90' notitle")
        else
            replot!(gp,ROI_spectrum_uuid,"with filledcurve y1=-0.1 lc rgb 'gray90' notitle")
        end             
    end

    # Plot global spectrum
    #
    calibrated_spectrum_uuid = register_data!(gp,global_fit_result.calibrated_spectrum)
    replot!(gp,calibrated_spectrum_uuid,"u 1:2 w l t 'Calibrated spectrum'")

    # Plot fit result
    #
    for idx_ROI in 1:n_ROI
        ROI_spectrum = get_spectrum(global_fit_result,idx_ROI)
        ROI_model = get_model(global_fit_result,idx_ROI)
        ROI_model_parameter = get_model_parameter(global_fit_result,idx_ROI)

        # ROI fit (only plot this is there are overlapping isotopicmotif)
        #
        if get_isotopicmotif_count(global_fit_result,idx_ROI) > 0
            ROI_Y_fit = eval_y(ROI_model,ROI_spectrum.X,ROI_model_parameter)
            ROI_id = register_data!(gp,hcat(ROI_spectrum.X,ROI_Y_fit))
            replot!(gp,ROI_id,"u 1:2 w l lw 2 dt 3 lc rgb 'blue' notitle")
        end
        
        # like there are potentially several isotopic motifs per ROI,
        # we use the visit method to perform individual drawings
        #
        visit(ROI_model,ROI_spectrum.Y,ROI_spectrum.X,ROI_model_parameter) do m,Y,X,θ
            if get_tagged_data_type(m) === EmbeddedData_IsotopicMotif_Model
                # get isotopic name & position using tagged data
                #
                embedded_data = get_tagged_data(m)
                name = get_isotopicmotif_name(global_fit_result,embedded_data._group_idx,embedded_data._local_idx)
                position = get_isotopicmotif_position(global_fit_result,embedded_data._group_idx,embedded_data._local_idx)

                # compute and plot model Y value
                #
                Y_fit = eval_y(m,X,θ)
                motif_id = register_data!(gp,hcat(X,Y_fit)) 
                replot!(gp,motif_id,"u 1:2 w l t '$name'")

                # plot the vertical 
                #
                add_vertical_line!(gp,position,name=name)
                
                return false
            end
            true
        end
    end
    
    write(gp_output_file,gp)
end

plot_fit("demo.gp",global_fit_result)

# Add a calibration shift
#
function add_calibration_shift!(global_fit_result::GlobalFitResult;scale=1)
    # Create the calibration map ================
    #
    map = Map_Translate(scale)

    # Create a new vector of LocalFit
    #
    map!(global_fit_result.local_fit,global_fit_result.local_fit) do lf
        calibrable_model = Model2Fit_Recalibration(lf.model,map)
        calibrable_model_θ = vcat(lf.θ,1)
        update_model(lf,calibrable_model,calibrable_model_θ)
    end
end

# Perform local fits and update θ
#
function local_fit!(global_fit_result::GlobalFitResult)
    conf = Levenberg_Marquardt_BC_Conf()

    for local_fit in global_fit_result.local_fit
        # Refit the model ----------------
        #
        # Create new bound with some margin...
        # TODO: maybe add this into Local_Fit <- This is mandatory!!!
        #
        θ = local_fit.θ

        abs_θ = @. max(1,abs(θ))
        bc = BoundConstraints(max.(0.0,θ-0.8* abs_θ),θ+1.2*abs_θ)

        ROI_model = local_fit.model
        ROI_X = local_fit.ROI_calibrated_spectrum.X
        ROI_Y = local_fit.ROI_calibrated_spectrum.Y
        
        nls = NLS_ForwardDiff_From_Model2Fit(ROI_model,ROI_X,ROI_Y)

        result = NLS_Solver.solve(nls,θ,bc,conf)

        # Update parameter value
        #
        if NLS_Solver.converged(result)
            local_fit.θ .= solution(result)
        else
            @warn "Local fit did not converged for ROI $(local_fit.data._group_idx)"
        end
    end

    global_fit_result
end

# ATTENTION fit ok , mais plot buggé 
add_calibration_shift!(global_fit_result,scale = 10)
local_fit!(global_fit_result)
# remove_calibration!(global_fit_result,scale = 10)

plot_fit("demo_local.gp",global_fit_result)

