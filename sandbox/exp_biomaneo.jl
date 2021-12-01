# TODO:
# - get_tagged_model_type & data_type
# - isomotif tag -> get local number
# - add gnuplot bar + name
# - perform local fit
# - replot
#
using Revise
using NLS_Fit
using NLS_Solver
using NLS_Models
using DelimitedFiles

# ================================================================

import NLS_Models: create_model

# Create the model associated to one group
#
# This model is the sum of Isotopic Model contained in the considered group
#
# This model is tagged (for the moment by its group idx)
#
struct Group_Model_EmbeddedData
    _group_idx::Int
end

function NLS_Models.create_model(grouped::GroupedBySupport{IsotopicMotif},idx_group::Int)

    model = Model2Fit_Empty()

    n_isotopicmotif = objects_in_group_size(grouped,idx_group)
    for idx_isotopicmotif in 1:n_isotopicmotif
        isotopic_motif = get_object(grouped,idx_group,idx_isotopicmotif)

        model += Peak_Motif(Gaussian_Peak(),get_profile_matrix(isotopic_motif))
    end

    model = Model2Fit_TaggedModel(model,Group_Model_EmbeddedData(idx_group))
   
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
                                   snip_halfwindow = 60,
                                   smoothing_halfwindow = 20)

spectrum = raw_spectrum - Y_baseline

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
stacked_models_σ_law_recalibration = Recalibration(stacked_models_σ_law,recalibration_map)

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

writedlm("poub_calibration.txt",hcat(spectrum.X,recalibrated_spectrum.X))

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
    data::Group_Model_EmbeddedData
    model::NLS_Fit.Abstract_Model2Fit
    ROI_calibrated_spectrum::Spectrum # CAVEAT: do not replace by
                                      # range, as this allows us to
                                      # use local recalibration
    θ::AbstractVector
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

# Give the number of isotopic motifs in ROI idx
function get_isotopicmotif_count(gfr::GlobalFitResult,idx::Int)
    objects_in_group_size(gfr.grouped,gfr.local_fit[idx].data._group_idx)
end 
    
# Extract fit result, grouping model per group and providing the right X,Y (after calibration)
#
# caveat: vector index has no reason to be the group index. The right
# way to get the group idx is to use data field.
#
# Model associated to each group is detected thanks to
#
#   Model2Fit_TaggedModel{ _ , Group_Model_EmbeddedData}
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
        if get_tagged_data_type(model) === Group_Model_EmbeddedData
            # process local model 
            data = NLS_Fit.get_data(model)
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
                                      global_model::Recalibration,
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

include("exp_gnuplot.jl")

function plot_fit(global_fit_result::GlobalFitResult)
    gp = GnuplotScript()


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
        if get_isotopicmotif_count(global_fit_result,idx_ROI) > 1
            ROI_Y_fit = eval_y(ROI_model,ROI_spectrum.X,ROI_model_parameter)
            ROI_id = register_data!(gp,hcat(ROI_spectrum.X,ROI_Y_fit))
            replot!(gp,ROI_id,"u 1:2 w l lw 2 dt 3 lc 'blue' notitle")
        end
        
        # like there are potentially several isotopic motifs per ROI,
        # we use the visit method to perform individual drawings
        #
        visit(ROI_model,ROI_spectrum.Y,ROI_spectrum.X,ROI_model_parameter) do m,Y,X,θ
            if m isa Peak_Motif
                Y_fit = eval_y(m,X,θ)
                motif_id = register_data!(gp,hcat(X,Y_fit))
                replot!(gp,motif_id,"u 1:2 w l notitle")

                return false
            end
            true
        end
    end
    
    write("demo.gp",gp)
end

plot_fit(global_fit_result)

# function perform_local_fit(for_local_fit)
#     conf = Levenberg_Marquardt_BC_Conf()

#     for ((m,Y,X,θ)) in for_local_fit
#         abs_θ = @. max(1,abs(θ))
#         bc = BoundConstraints(θ-0.8* abs_θ,θ+1.2*abs_θ)
#         nls = NLS_ForwardDiff_From_Model2Fit(m,X,Y)
#         result = NLS_Solver.solve(nls,θ,bc,conf)

#         println(result)
#     end
# end 

#nperform_local_fit(for_local_fit)
