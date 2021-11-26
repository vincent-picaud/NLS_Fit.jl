# Cet exemple est complex
#
# Il faut empiler les modèles en fonction des ROI Mais comme la
# fonction de calibration dépend d'un X = 1:N entier, il faut
# *anticiper* prévoir le même découpage pour cet X entier

using Revise
using NLS_Fit
using NLS_Solver
using NLS_Models

# ================================================================

import NLS_Models: create_model

# Create the model associated to one group
#
# This model is the sum of Isotopic Model contained in the considered group
#
function NLS_Models.create_model(grouped::GroupedBySupport{IsotopicMotif},idx_group::Int)

    model = Model2Fit_Empty()

    n_isotopicmotif = objects_in_group_size(grouped,idx_group)
    for idx_isotopicmotif in 1:n_isotopicmotif
        isotopic_motif = get_object(grouped,idx_group,idx_isotopicmotif)

        model += Peak_Motif(Gaussian_Peak(),get_profile_matrix(isotopic_motif))
    end

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

# Perform fit
# ================

# Inputs
#
spectrum = read_spectrum_Biomaneo("/home/picaud/GitHub/NLS_Models.jl/data/0000000001.txt")
spectrum.Y ./= maximum(spectrum.Y)
vect_of_isotopicmotif = hardcoded_IsotopicMotifVect()

# Remove baseline
#
Y_baseline = compute_baseline_snip(spectrum,
                                   snip_halfwindow = 60,
                                   smoothing_halfwindow = 20)

spectrum = spectrum - Y_baseline

# Create ROI model & spectrum
#
grouped = groupbysupport(vect_of_isotopicmotif,by=get_ROI_interval)
stacked_models,ROI_spectrum = create_stacked_model_ROI_spectrum_pair(grouped,spectrum)

# create θ : (h,σ) x number of ROIs
#
n_θ = NLS_Fit.parameter_size(stacked_models)
θ_init = zeros(n_θ) .+ 0.5
θ_lb = zeros(n_θ)
θ_ub = zeros(n_θ) .+ 6

# prepare solver
#
bc = BoundConstraints(θ_lb,θ_ub)
nls = NLS_ForwardDiff_From_Model2Fit(stacked_models,ROI_spectrum.X,ROI_spectrum.Y)
conf = Levenberg_Marquardt_BC_Conf()

result = NLS_Solver.solve(nls,θ_init,bc,conf)

Y_fit = eval_y(stacked_models,ROI_spectrum.X,solution(result))

# gnuplot
#
using DelimitedFiles

writedlm("poub.txt",hcat(ROI_spectrum.X,ROI_spectrum.Y,Y_fit))

