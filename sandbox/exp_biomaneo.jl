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
function NLS_Models.create_model(grouped::GroupedBySupport{IsotopicMotif})
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
function extract_ROI_spectrum(grouped::GroupedBySupport{IsotopicMotif},sp)
    ROI_interval = group_support_as_interval(grouped)
    ROI_range = create_range_from_interval(ROI_interval,spectrum.X)

    Spectrum(extract_ROIs(ROI_range,sp.X),
              extract_ROIs(ROI_range,sp.Y)),
    length.(ROI_range)
    
end 



# ****************************************************************

# Perform fit 
spectrum = read_spectrum_Biomaneo("/home/picaud/GitHub/NLS_Models.jl/data/0000000095.txt")
spectrum.Y ./= maximum(spectrum.Y)
vect_of_isotopicmotif = hardcoded_IsotopicMotifVect()
grouped = groupbysupport(vect_of_isotopicmotif,by=get_ROI_interval)

create_model(grouped)

extract_ROI_spectrum(grouped,spectrum)

