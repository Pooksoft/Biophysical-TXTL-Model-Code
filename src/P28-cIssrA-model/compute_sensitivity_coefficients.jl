# include the includes -
include("Include.jl")

# mean center -
function mean_center_array(results_array::Array{Float64,2})::Array{Float64,2}

    # get the size -
    (NR,NC) = size(results_array)
    scaled_array = zeros(NR,NC)

    for col_index = 1:NC
        
        data_col = results_array[:,col_index]
        mu_value = mean(data_col)
        std_value = std(data_col)

        for row_index = 1:NR
            scaled_array[row_index,col_index] = (data_col[row_index] - mu_value)/(std_value)
        end
    end

    return scaled_array
end

# computes the model performance -
function model_performance(parameter_guess_array,index)

    # what is the host_type?
    host_type = :cell_free

    # load the default data_dictionary -
    time_start = 0.0
    time_stop = 16.0
    time_step_size = 0.01

    # what is the host_type?
    host_type = :cell_free

    # path to parameters -
    path_to_biophysical_constants_file = "./CellFree.json"
    #path_to_data_dir = "$(pwd())/data"

    # Load the data dictionary (uses the default biophysical_constants file)
    default_data_dictionary = build_data_dictionary((time_start,time_stop,time_step_size), path_to_biophysical_constants_file, host_type)

    # Set some values for the weights, gene lengths etc -
    model_data_dictionary = customize_data_dictionary(default_data_dictionary,host_type)

    # load the experimental data -
    # exp_data_dictionary = load_experimental_data_dictionary(path_to_data_dir)

    # update the paramaters in the model data dictionary -
    # for now - lets only search over dG's -
    R = model_data_dictionary["R"]
    T_K = model_data_dictionary["T_K"]

    # what is the size of the parameter_guess_array?
    number_of_parameters = length(parameter_guess_array)

    # Phase 1: parameter update =========================================================================== #
    # update the paramaters in the model data dictionary -
    # for now - lets only search over dG's -
    R = model_data_dictionary["R"]
    T_K = model_data_dictionary["T_K"]

    # compute W -
    tmp_W_array = Float64[]
    for index = 1:8
        parameter_guess = parameter_guess_array[index]
        value = exp(-1*parameter_guess/(R*T_K))
        push!(tmp_W_array,value)
    end

    # update the control W's -
    control_parameter_dictionary = model_data_dictionary["control_parameter_dictionary"]
    control_parameter_dictionary["W_clssrA_RNAP"] = tmp_W_array[1]          # 1
	control_parameter_dictionary["W_clssrA_sigma_28"] = tmp_W_array[2]      # 2
	control_parameter_dictionary["W_deGFPssrA_RNAP"] = tmp_W_array[3]       # 3
	control_parameter_dictionary["W_deGFPssrA_sigma_70"] = tmp_W_array[4]   # 4
	control_parameter_dictionary["W_deGFPssrA_clssrA"] = tmp_W_array[5]     # 5
	control_parameter_dictionary["W_sigma_28_RNAP"] = tmp_W_array[6]        # 6
	control_parameter_dictionary["W_sigma_28_sigma_70"] = tmp_W_array[7]    # 7
	control_parameter_dictionary["W_sigma_28_clssrA"] = tmp_W_array[8]      # 8
    model_data_dictionary["control_parameter_dictionary"] = control_parameter_dictionary

    binding_parameter_dictionary = model_data_dictionary["binding_parameter_dictionary"]
    binding_parameter_dictionary["n_clssrA_sigma_28"] = parameter_guess_array[9]
	binding_parameter_dictionary["K_clssrA_sigma_28"] = parameter_guess_array[10]
	binding_parameter_dictionary["n_deGFPssrA_sigma_70"] = parameter_guess_array[11]
	binding_parameter_dictionary["K_deGFPssrA_sigma_70"] = parameter_guess_array[12]
	binding_parameter_dictionary["n_deGFPssrA_clssrA"] = parameter_guess_array[13]
	binding_parameter_dictionary["K_deGFPssrA_clssrA"] = parameter_guess_array[14]
	binding_parameter_dictionary["n_sigma_28_sigma_70"] = parameter_guess_array[15]
	binding_parameter_dictionary["K_sigma_28_sigma_70"] = parameter_guess_array[16]
	binding_parameter_dictionary["n_sigma_28_clssrA"] = parameter_guess_array[17]
    binding_parameter_dictionary["K_sigma_28_clssrA"] = parameter_guess_array[18]
    model_data_dictionary["binding_parameter_dictionary"] = binding_parameter_dictionary

    # time constant modifier -
    time_constant_modifier_array = [
		0.0	                          ;	  # 1	clssrA
		0.0	                          ;	  # 2	deGFPssrA
		0.0	                          ;	  # 3	sigma_28
		0.0                           ;	  # 4	sigma_70
		parameter_guess_array[19]     ;	  # 5	mRNA_clssrA
		parameter_guess_array[20]     ;	  # 6	mRNA_deGFPssrA
		parameter_guess_array[21]     ;	  # 7	mRNA_sigma_28
		1.0 	                      ;	  # 8	mRNA_sigma_70
		parameter_guess_array[22]	  ;	  # 9	protein_clssrA
		parameter_guess_array[23]	  ;	  # 10	protein_deGFPssrA
		parameter_guess_array[24]	  ;	  # 11	protein_sigma_28
		1.0	                          ;	  # 12	protein_sigma_70
	]
    model_data_dictionary["time_constant_modifier_array"] = time_constant_modifier_array

    # setup degradation_modifier_array -
    degradation_modifier_array = [
		0.0	                            ;	# 1	clssrA
		0.0	                            ;	# 2	deGFPssrA
		0.0	                            ;	# 3	sigma_28
		0.0	                            ;	# 4	sigma_70
		parameter_guess_array[25]	    ;	# 5	mRNA_clssrA 0.0085
		parameter_guess_array[26]		;	# 6	mRNA_deGFPssrA 0.05
		parameter_guess_array[27]		;	# 7	mRNA_sigma_28
		0.0                             ;   # 8 mRNA_sigma_70                       
		parameter_guess_array[28]		;	# 9	protein_clssrA
		parameter_guess_array[29]		;	# 10	protein_deGFPssrA
		parameter_guess_array[30]		;	# 11	protein_sigma_28
	    parameter_guess_array[31]		;	# 12	protein_sigma_70
    ]
    
    # update the translation time -
    model_data_dictionary["half_life_translation_capacity"] = parameter_guess_array[32]
    
    # lastly, update KL -
    biophysical_constants_dictionary = model_data_dictionary["biophysical_constants_dictionary"]
    biophysical_constants_dictionary["translation_saturation_constant"] = parameter_guess_array[33]
    model_data_dictionary["biophysical_constants_dictionary"] = biophysical_constants_dictionary

    # grab defaults -
    species_symbol_type_array = model_data_dictionary["species_symbol_type_array"]
    protein_coding_length_array = model_data_dictionary["protein_coding_length_array"]
    gene_coding_length_array = model_data_dictionary["gene_coding_length_array"]
    time_constant_modifier_array = model_data_dictionary["time_constant_modifier_array"]
    initial_condition_array = model_data_dictionary["initial_condition_array"]

    # # get gene IC -
    idx_gene = findall(x->x==:gene,species_symbol_type_array)
    gene_abundance_array = initial_condition_array[idx_gene]

    # Precompute the translation parameters -
	translation_parameter_array = precompute_translation_parameter_array(biophysical_constants_dictionary, protein_coding_length_array, time_constant_modifier_array,host_type)
    model_data_dictionary["translation_parameter_array"] = translation_parameter_array

	# Precompute the kinetic limit of transcription -
	transcription_kinetic_limit_array = precompute_transcription_kinetic_limit_array(biophysical_constants_dictionary, gene_coding_length_array, gene_abundance_array, time_constant_modifier_array, host_type)
    model_data_dictionary["transcription_kinetic_limit_array"] = transcription_kinetic_limit_array

    # Dilution degrdation matrix -
    dilution_degradation_matrix = build_dilution_degradation_matrix(biophysical_constants_dictionary,species_symbol_type_array,degradation_modifier_array)
    model_data_dictionary["dilution_degradation_matrix"] = dilution_degradation_matrix
    # ===================================================================================================== #

    # Phase 2:  solve model equations ===================================================================== #
    # solve the balance equations -
    (TSIM,XSIM) = SolveBalances(time_start,time_stop,time_step_size,model_data_dictionary)
    # ===================================================================================================== #

    # Phase 3: compute the model performance metrics ====================================================== #
    p_GFP_AUC = integrate(TSIM,XSIM[:,index])
    # ===================================================================================================== #

    # return the performance_array -
    return p_GFP_AUC
end

# Function to set some values for the weights, gene lengths etc -
function customize_data_dictionary(default_data_dictionary::Dict{String,Any}, host_type::Symbol)::Dict{String,Any}

    # make a deepcopy -
    customized_data_dictionary = deepcopy(default_data_dictionary)

    # set the gene lengths -
    gene_length_array = customized_data_dictionary["gene_coding_length_array"]
    gene_length_array[1] = 744.0    # clssrA
    gene_length_array[2] = 711.0    # deGFPssrA
    gene_length_array[3] = 720.0    # S28
    gene_length_array[4] = 720.0    # S70 => doesn't matter, gene expression = 0

    # set the protein lengths -
    protein_coding_length_array = customized_data_dictionary["protein_coding_length_array"]
    protein_coding_length_array[1] = 248.0  # clssrA
    protein_coding_length_array[2] = 237.0  # deGFPssrA
    protein_coding_length_array[3] = 240.0  # S28
    protein_coding_length_array[4] = 240.0  # S70 => doesn't matter, translation = 0

    # set ICs -
    initial_condition_array = customized_data_dictionary["initial_condition_array"]
    initial_condition_array[1] = 0.001  # clssrA gene 
    initial_condition_array[2] = 0.008  # deGFP gene
    initial_condition_array[3] = 0.0015 # S28 gene
    initial_condition_array[4] = 0.0    # no S70 gene, just S70 protein in the extract
    #initial_condition_array[11] = 0.020 # sigma28 protein in muM
    initial_condition_array[12] = 0.035 # S70 protein in muM
    
    # setup the W's -
    control_parameter_dictionary = customized_data_dictionary["control_parameter_dictionary"]
    control_parameter_dictionary["W_clssrA_RNAP"] = 0.000014
	control_parameter_dictionary["W_clssrA_sigma_28"] = 1.0
	control_parameter_dictionary["W_deGFPssrA_RNAP"] = 0.000014
	control_parameter_dictionary["W_deGFPssrA_sigma_70"] = 10.0
	control_parameter_dictionary["W_deGFPssrA_clssrA"] = 100.0
	control_parameter_dictionary["W_sigma_28_RNAP"] = 0.00014
	control_parameter_dictionary["W_sigma_28_sigma_70"] = 100.0
	control_parameter_dictionary["W_sigma_28_clssrA"] = 0.1
	control_parameter_dictionary["W_sigma_70_RNAP"] = 0.0

    # setup the binding parameters -
    bP = readdlm("bP-Opt.dat")
    binding_parameter_dictionary = customized_data_dictionary["binding_parameter_dictionary"]
    bp_symbol_name_array = String[]
    push!(bp_symbol_name_array, "n_clssrA_sigma_28")
    push!(bp_symbol_name_array, "K_clssrA_sigma_28")
    push!(bp_symbol_name_array, "n_deGFPssrA_sigma_70")
    push!(bp_symbol_name_array, "K_deGFPssrA_sigma_70")
    push!(bp_symbol_name_array, "n_deGFPssrA_clssrA")
    push!(bp_symbol_name_array, "K_deGFPssrA_clssrA")
    push!(bp_symbol_name_array, "n_sigma_28_sigma_70")
    push!(bp_symbol_name_array, "K_sigma_28_sigma_70")
    push!(bp_symbol_name_array, "n_sigma_28_clssrA")
    push!(bp_symbol_name_array, "K_sigma_28_clssrA")
    for (index,key_text) in enumerate(bp_symbol_name_array)
        
        # get value -
        value = bP[index]

        # cache -
        binding_parameter_dictionary[key_text] = value
    end
    
    # setup degradation_modifier_array -
    degradation_modifier_array = [
		0.0	    ;	# 1	clssrA
		0.0	    ;	# 2	deGFPssrA
		0.0	    ;	# 3	sigma_28
		0.0	    ;	# 4	sigma_70
		0.05    ;	# 5	mRNA_clssrA 0.0085
		0.05	;	# 6	mRNA_deGFPssrA 0.05
		0.05	;	# 7	mRNA_sigma_28
		1.0	    ;	# 8	mRNA_sigma_70
		30.0	;	# 9	protein_clssrA
		28.0	;	# 10	protein_deGFPssrA
		1.5	    ;	# 11	protein_sigma_28
		1.0	    ;	# 12	protein_sigma_70
	]
    customized_data_dictionary["degradation_modifier_array"] = degradation_modifier_array
    biophysical_constants_dictionary = customized_data_dictionary["biophysical_constants_dictionary"]
    species_symbol_type_array = customized_data_dictionary["species_symbol_type_array"]
    protein_coding_length_array = customized_data_dictionary["protein_coding_length_array"]
    gene_coding_length_array = customized_data_dictionary["gene_coding_length_array"]

    # get gene IC -
    idx_gene = findall(x->x==:gene,species_symbol_type_array)
    gene_abundance_array = initial_condition_array[idx_gene]

    # Dilution degrdation matrix -
	dilution_degradation_matrix = build_dilution_degradation_matrix(biophysical_constants_dictionary,species_symbol_type_array,degradation_modifier_array)
    customized_data_dictionary["dilution_degradation_matrix"] = dilution_degradation_matrix

    # time constant modifier -
    time_constant_modifier_array = [
		0.0	                    ;	  # 1	clssrA
		0.0	                    ;	  # 2	deGFPssrA
		0.0	                    ;	  # 3	sigma_28
		0.0                     ;	  # 4	sigma_70
        0.5549555127597618      ;	  # 5	mRNA_clssrA
        0.9961201331095983      ;	  # 6	mRNA_deGFPssrA
		0.9573079006791378      ;	  # 7	mRNA_sigma_28
        2.47595196225463	    ;	  # 8	mRNA_sigma_70
		1.2340314042937286	    ;	  # 9	protein_clssrA
		2.0833904410438735 	    ;	  # 10	protein_deGFPssrA
		4.99999996459193	    ;	  # 11	protein_sigma_28
        2.47595196225463	    ;	  # 12	protein_sigma_70
	]
    customized_data_dictionary["time_constant_modifier_array"] = time_constant_modifier_array

    # Precompute the translation parameters -
	translation_parameter_array = precompute_translation_parameter_array(biophysical_constants_dictionary, protein_coding_length_array, time_constant_modifier_array,host_type)
    customized_data_dictionary["translation_parameter_array"] = translation_parameter_array

	# Precompute the kinetic limit of transcription -
	transcription_kinetic_limit_array = precompute_transcription_kinetic_limit_array(biophysical_constants_dictionary, gene_coding_length_array, gene_abundance_array, time_constant_modifier_array, host_type)
    customized_data_dictionary["transcription_kinetic_limit_array"] = transcription_kinetic_limit_array

    # return -
    return customized_data_dictionary
end


function main(path_to_ensemble_file::String,index)

    # setup the sensitivity function -
    SF(P) = model_performance(P,index)

    # setup ranges -
    sample_bounds_array = Array{Tuple,1}()
    ensemble_array = readdlm(path_to_ensemble_file)
    (number_of_parameters,number_of_trials) = size(ensemble_array)
    for parameter_index = 1:(number_of_parameters-1)

        # get row of parameters -
        parameter_row = ensemble_array[parameter_index,:]
        min_value = minimum(parameter_row)
        max_value = maximum(parameter_row)

        # create the tuple -
        tmp_tuple = (min_value,max_value)

        # cache -
        push!(sample_bounds_array,tmp_tuple)
    end

    @show index

    # do the global sensitivity analysis -
    sensitivity_results = gsa(SF,Morris(total_num_trajectory=10000,num_trajectory=1000),sample_bounds_array)

    # return -
    return sensitivity_results
end

# setup paths -
path_to_ensemble_file = "$(pwd())/Ensemble-c1-restriction-T20.dat"

# compute a sensitivity array for the AUC of each species -
species_index_array = [5 6 7 9 10 11]
number_of_species = length(species_index_array)
number_of_parameters = 33
results_array = zeros(number_of_parameters,1)
for species_index in species_index_array
    
    global results_array
    
    # conduct senstivity analysis -
    sensitivity_results = main(path_to_ensemble_file,species_index)

    # get the μ and σ^2
    mu = sensitivity_results.means
    var = sensitivity_results.variances

    #@show mu, var

    results_array = [results_array transpose(mu) transpose(var)]

end

results_array = results_array[:,2:end]