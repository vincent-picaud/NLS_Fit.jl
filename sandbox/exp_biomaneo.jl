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
    _name::String
    _position::Float64
end

# Used to tag a complete group model
#
# This model is the sum of Isotopic Model contained in the considered group
#
struct EmbeddedData_Group_Model
    # maybe one can add number of iso motif
    _group_interval::Interval 
end

function NLS_Models.create_model(grouped::GroupedBySupport{IsotopicMotif},idx_group::Int)

    model = Model2Fit_Empty()

    n_isotopicmotif = objects_in_group_size(grouped,idx_group)
    for idx_isotopicmotif in 1:n_isotopicmotif
        isotopic_motif = get_object(grouped,idx_group,idx_isotopicmotif)

        profile_matrix = get_profile_matrix(isotopic_motif)

        metadata = EmbeddedData_IsotopicMotif_Model(get_name(isotopic_motif),
                                                    get_position(isotopic_motif))
        
        model += Model2Fit_TaggedModel(Peak_Motif(Gaussian_Peak(),profile_matrix),
                                       metadata)
    end

    group_interval = group_support_as_interval(grouped,idx_group)
    model = Model2Fit_TaggedModel(model,EmbeddedData_Group_Model(group_interval))
   
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
# raw_spectrum = read_spectrum_Biomaneo("/home/picaud/Data/Spectres_Biomaneo/txt-minos NMOG/0000128803(d5).txt")
# raw_spectrum = read_spectrum_Biomaneo("/home/picaud/Data/Spectres_Biomaneo/Reunion_25_Oct/Data_Input/Validation_04/129913.txt")
# raw_spectrum = read_spectrum_Biomaneo("/home/picaud/Data/Spectres_Biomaneo/Spectres_NewBorn/2221206155(d2).txt")

# raw_spectrum = read_spectrum_Biomaneo("/home/picaud/Data/Spectres_Biomaneo/Spectres_Biomaneo_MF/Heterozygote HbE/0000000036_digt_MF.txt")
# raw_spectrum = read_spectrum_Biomaneo("/home/picaud/Data/Spectres_Biomaneo/January_2020_normalized/Heterozygote HbE B Thal/0000000017_digt_0001_J4_(Manual)_19-12-20_14-19_0001.txt")
raw_spectrum = read_spectrum_Biomaneo("/home/picaud/GitHub/NLS_Models.jl/data/0000000095.txt")
#raw_spectrum = read_spectrum_Biomaneo("/home/picaud/GitHub/NLS_Models.jl/data/0000000001.txt")
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
#
recalibration_map = Map_Affine(ROI_spectrum.X[1],ROI_spectrum.X[end])

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

# Solve the problem
# ================
bc = BoundConstraints(θ_lb,θ_ub)
nls = NLS_ForwardDiff_From_Model2Fit(stacked_models_σ_law_recalibration,ROI_spectrum.X,ROI_spectrum.Y)
conf = Levenberg_Marquardt_BC_Conf()

result = NLS_Solver.solve(nls,θ_init,bc,conf)

# # Plot solution
# # ================
# Y_fit = eval_y(stacked_models_σ_law_recalibration,ROI_spectrum.X,solution(result))

# recalibrated_spectrum = Spectrum(eval_calibrated_x(stacked_models_σ_law_recalibration,spectrum.X,solution(result)),spectrum.Y)

# writedlm("poub_calibration.txt",hcat(spectrum.X,recalibrated_spectrum.X,spectrum.Y))


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
            #
            group_interval = data._group_interval
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





# Plot routine ****************************************************************

# A structure to store model value: this struct in only used in plot_fit() function
#
struct ROI_IsotopicModel_ExtractedData
    Y::AbstractVector
    name::String
    position::Float64
end

function plot_fit(gp::GnuplotScript,local_fit_vect::AbstractVector{LocalFit})

    # Transparent
    #
    free_form(gp,"set style fill transparent solid 0.5 noborder")
    
    # Plot local fits result
    #
    for local_fit in local_fit_vect

        ROI_spectrum = local_fit.ROI_calibrated_spectrum
        ROI_model = local_fit.model
        ROI_model_parameter = local_fit.θ

        ROI_spectrum_calibrated = deepcopy(ROI_spectrum)
        ROI_local_model_data = ROI_IsotopicModel_ExtractedData[]
        
        visit(ROI_model,ROI_spectrum.Y,ROI_spectrum.X,ROI_model_parameter) do m,Y,X,θ
            # if a recalibration is performed, record it for future plot
            if m isa NLS_Fit.Model2Fit_Recalibration
                ROI_spectrum_calibrated = Spectrum(eval_calibrated_x(m,X,θ),Y)
                return true
            end
            
            if get_tagged_data_type(m) === EmbeddedData_IsotopicMotif_Model
                # get isotopic name & position using tagged data
                #
                embedded_data = get_tagged_data(m)
                name = embedded_data._name
                position = embedded_data._position

                # compute model value and store it for future use
                # (plotting...)
                #
                Y_fit = eval_y(m,X,θ)
                push!(ROI_local_model_data,
                      ROI_IsotopicModel_ExtractedData(Y_fit,name,position))

                
                return false
            end
            true
        end

        # Concatenate all data vectors X,Y1,Y2... in a single matrix
        #
        # Note by puttingisotopicModel_data we preserve the identification:
        #
        #    "curve number" <-> index of all_ROI_isotopicModel_extractedData
        #
        data = hcat( (isotopicModel_data.Y for isotopicModel_data in ROI_local_model_data) ... ,
                     ROI_spectrum_calibrated.X,
                     ROI_spectrum_calibrated.Y)
        data_m = size(data,2)
        data_m_last = data_m -2 # number of the last model
        data_m_X = data_m-1
        data_m_Y = data_m
        
        data_id = register_data!(gp,data)
        replot!(gp,data_id,"u $data_m_X:$data_m_Y with filledcurve y1=-0.1 lc rgb 'gray90' notitle")

        for j in 1:data_m_last
            name     = ROI_local_model_data[j].name
            position = ROI_local_model_data[j].position

            replot!(gp,data_id,"u $data_m_X:$j w l lw 2 t '$name'")
                            # plot the vertical 
                
            add_vertical_line!(gp,position,name=name)
        end

        # Sum all
        #    reduce((l,r)->l*"+"*r,map(x->"\$$x",1:4))
        # produce "\$1+\$2+\$3+\$4"
        #
        if data_m_last>1
            cmd_str = reduce((l,r)->l*"+"*r,map(x->"\$$x",1:data_m_last))
            cmd_str = "u $data_m_X:("*cmd_str*") w l dt 3 lc rgb 'blue' notitle"
            replot!(gp,data_id,cmd_str)
        end 
    end

    gp
end

function plot_fit(gp_output_file::String,
                  uncalibrated_spectrum::Spectrum,
                  local_fit_vect::AbstractVector{LocalFit})
    gp = GnuplotScript()

        
    # set title
    # 
    free_form(gp,"set title '$gp_output_file' noenhanced")

    sp_id = register_data!(gp,hcat(uncalibrated_spectrum.X,uncalibrated_spectrum.Y))
    replot!(gp,sp_id,"u 1:2 with filledcurve y1=-0.1 lc rgb 'gray50' t 'uncalibrated spectrum'")

    plot_fit(gp,local_fit_vect)
    write(gp_output_file,gp)

    gp
end 

# ****************************************************************

# Plot result with global calibration
#
plot_fit("demo.gp",spectrum, global_fit_result.local_fit)

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

# Perform local fit ================
#
add_calibration_shift!(global_fit_result,scale = 10)
local_fit!(global_fit_result)

# Plot result
#
plot_fit("demo_local.gp",spectrum,global_fit_result.local_fit)
