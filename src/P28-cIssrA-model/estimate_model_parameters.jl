# include packages -
include("Include.jl")

function neighbor_function!(proposed_parameter_array, parameter_array)

    # setup -
    sigma = 0.05
    for i in eachindex(parameter_array)
        proposed_parameter_array[i] = parameter_array[i] + sigma*randn()
    end
end

function cooling_function(temperature)

  # define my new temperature -
  alpha = 0.9
  return alpha*temperature
end


function objective_function(parameter_guess_array,time_start,time_step_size,time_stop,model_data_dictionary,exp_data_dictionary)

    # what is the host_type?
    host_type = :cell_free

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

    # Phase 3:  compute simulation error ================================================================== #
    # compute the error - we need to do a bunch of interpolation -
    tsim_exp = exp_data_dictionary["deGFP_mRNA_data_table"][!,:Time_in_hr]
    
    # GFP mRNA -
    itp_gfp_mRNA =  LinearInterpolation(TSIM, (1000)*XSIM[:,6]);
    mRNA_GFP_sim = itp_gfp_mRNA[tsim_exp]  # convert to muM from nM

    # sigma 28 mRNA -
    itp_sigma28_mRNA =  LinearInterpolation(TSIM, (1000)*XSIM[:,7]);
    mRNA_sigma28_sim = itp_sigma28_mRNA[tsim_exp] # convert to muM from nM

    # C1 mRNA -
    itp_C1_mRNA =  LinearInterpolation(TSIM, (1000)*XSIM[:,5]);
    mRNA_C1_sim = itp_C1_mRNA[tsim_exp] # convert to muM from nM

    # GFP protein -
    itp_gfp_protein =  LinearInterpolation(TSIM, XSIM[:,10]);
    protein_GFP_sim = itp_gfp_protein[tsim_exp]

    # get experimental data -
    mRNA_GFP_exp = exp_data_dictionary["deGFP_mRNA_data_table"][!,:Concentration_in_nM]
    mRNA_GFP_std_exp = exp_data_dictionary["deGFP_mRNA_data_table"][!,:stdev]

    mRNA_sigma28_exp = exp_data_dictionary["sigma28_mRNA_data_table"][!,:Concentration_in_nM]
    mRNA_sigma28_std_exp = exp_data_dictionary["sigma28_mRNA_data_table"][!,:stdev]

    mRNA_C1_exp = exp_data_dictionary["c1_ssrA_mRNA_data_table"][!,:Concentration_in_nM]
    mRNA_C1_std_exp = exp_data_dictionary["c1_ssrA_mRNA_data_table"][!,:stdev]

    protein_GFP_exp = exp_data_dictionary["deGFP_protein_data_table"][!,:Concentration_in_muM]
    protein_GFP_std_exp = exp_data_dictionary["deGFP_protein_data_table"][!,:stdev]

    # compute error terms -
    error_term_array = zeros(4)
    
    # mRNA GFP -
    #WM1 = diagm(1.0 ./((mRNA_GFP_std_exp).^2))
    error_vector_1 = (mRNA_GFP_exp .- mRNA_GFP_sim)
    error_term_array[1] = transpose(error_vector_1)*error_vector_1

    # mRNA sigma28 -
    #WM2 = diagm(1.0 ./((mRNA_sigma28_std_exp).^2))
    error_vector_2 = (mRNA_sigma28_exp .- mRNA_sigma28_sim)
    error_term_array[2] = transpose(error_vector_2)*error_vector_2

    # mRNA C1 -
    #WM3 = diagm( 1.0 ./((mRNA_C1_std_exp).^2))
    error_vector_3 = (mRNA_C1_exp .- mRNA_C1_sim)
    error_term_array[3] = transpose(error_vector_3)*error_vector_3

    # protein deGFP -
    #WM4 = diagm( 1.0 ./((mRNA_C1_std_exp).^2) )
    error_vector_4 = (protein_GFP_exp .- protein_GFP_sim)
    error_term_array[4] = transpose(error_vector_4)*error_vector_4
    # ===================================================================================================== #

    # error total -
    error_total = sum(error_term_array)

    #@show error_total
    # return -
    return error_total
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


function main(path_to_data_dir::String, parameter_array; algorithm::Symbol = :LBFGS)

    # load the default data_dictionary -
    time_start = 0.0
    time_stop = 16.0
    time_step_size = 0.01

    # what is the host_type?
    host_type = :cell_free

    # path to parameters -
    path_to_biophysical_constants_file = "./CellFree.json"

    # Load the data dictionary (uses the default biophysical_constants file)
    default_data_dictionary = build_data_dictionary((time_start,time_stop,time_step_size), path_to_biophysical_constants_file, host_type)

    # Set some values for the weights, gene lengths etc -
    model_data_dictionary = customize_data_dictionary(default_data_dictionary,host_type)

    # load the experimental data -
    exp_data_dictionary = load_experimental_data_dictionary(path_to_data_dir)

    # setup the objective function -
    OF(P) = objective_function(P,time_start,time_step_size,time_stop,model_data_dictionary,exp_data_dictionary)

    # setup paramter bounds -
    pvec_bounds = [

        # dG's -
        -50000.0  -100.0    ;   # 1 W_clssrA_RNAP
        50000.0  110000.0   ;   # 2 W_clssrA_sigma_28
        39995.69565886529   48282.2841495398    ;   # 3 W_deGFPssrA_RNAP
        -32689.208748865043  -24993.98186199351    ;   # 4 W_deGFPssrA_sigma_70
        -50000.0  20000.0   ;   # 5 W_deGFPssrA_clssrA
        12000.0  25000.0    ;   # 6 W_sigma_28_RNAP
        -50000.0  25000.0   ;   # 7 W_sigma_28_sigma_70
        -50000.0  25000.0   ;   # 8 W_sigma_28_clssrA

        # binding parameters -
        0.5 10.0            ;   # 9 n_clssrA_sigma_28
        0.001 100.0         ;   # 10 K_clssrA_sigma_28
        0.5037004212911795  1.7091831951968874       ;   # 11 n_deGFPssrA_sigma_70
        25.539812443045424 63.99724158428853         ;   # 12 K_deGFPssrA_sigma_70
        0.5 10.0            ;   # 13 n_deGFPssrA_clssrA
        0.001 100.0         ;   # 14 K_deGFPssrA_clssrA
        0.5 10.0            ;   # 15 n_sigma_28_sigma_70
        0.001 100.0         ;   # 16 K_sigma_28_sigma_70
        0.5 10.0            ;   # 17 n_sigma_28_clssrA
        0.001 100.0         ;   # 18 K_sigma_28_clssrA

        # time constants -
		0.001 100.0         ;   # 19	mRNA_clssrA
		0.0039273045649573235 7.17667775643900         ;	# 20	mRNA_deGFPssrA
		0.001 100.0         ;	# 21	mRNA_sigma_28
		0.001 100.0         ;	# 22	protein_clssrA
		0.35035678021027006 10.203484322453043         ;	# 23	protein_deGFPssrA
		0.001 100.0	        ;	# 24	protein_sigma_28
        
        # degradation mods -
        0.1 10.0            ;	# 25	mRNA_clssrA 0.0085
		0.0011598053735671 3.6284714410194088 	     ;	# 26	mRNA_deGFPssrA 0.05
		0.1 10.0 	        ;	# 27	mRNA_sigma_28
		1.0 100.0 	        ;	# 28	protein_clssrA
		0.03218782255894215 8.897509480666375 	        ;	# 29	protein_deGFPssrA
		0.13 10.0 	        ;	# 30	protein_sigma_28
        0.132047561240967526 10.137258626014857 	        ;	# 31	protein_sigma_70

        # w -
        6.0 10.0            ;   # 32    translation capacity half-life

        # KL -
        10.81768297340375 276.3396750321431         ;   # 33    KL in muM
    ];


    # pvec_initial = Float64[]
    # dG = readdlm("dG-Opt.dat")
    # R = 8.314
    # T_K = (273.15+29)
    # W_array = Float64[]
    # for value in dG
    #     #tmp = exp(-value/(R*T_K))
    #     push!(pvec_initial, (1+0.10*rand())*value)
    # end
    
    # # binding parameters -
    # # bP = readdlm("bP-Opt.dat")
    # # for value in bP
    # #     push!(pvec_initial, (1+0.10*rand())*value)
    # # end

    # push!(pvec_initial, 4.0)      # 9 n_clssrA_sigma_28
    # push!(pvec_initial, 0.001)    # 10 K_clssrA_sigma_28
    # push!(pvec_initial, 1.0)      # 11 n_deGFPssrA_sigma_70
    # push!(pvec_initial, 0.008)    # 12 K_deGFPssrA_sigma_70
    # push!(pvec_initial, 1.00)     # 13 n_deGFPssrA_clssrA
    # push!(pvec_initial, 10.97)    # 14 K_deGFPssrA_clssrA
    # push!(pvec_initial, 2.83)     # 15 n_sigma_28_sigma_70
    # push!(pvec_initial, 0.0227)   # 16 K_sigma_28_sigma_70
    # push!(pvec_initial, 1.43)     # 17 n_sigma_28_clssrA
    # push!(pvec_initial, 0.69)     # 18 K_sigma_28_clssrA
    
    # # time constants -
    # # tC = readdlm("tC-Opt.dat")
    # # for value in tC
    # #     push!(pvec_initial, (1+0.10*rand())*value)
    # # end
	# push!(pvec_initial, 1.5)      # 19	mRNA_clssrA
	# push!(pvec_initial, 2.47)     # 20	mRNA_deGFPssrA
	# push!(pvec_initial, 0.71)     # 21	mRNA_sigma_28
	# push!(pvec_initial, 0.63)	    # 22	protein_clssrA
	# push!(pvec_initial, 2.01)	    # 23	protein_deGFPssrA
	# push!(pvec_initial, 5.0)	    # 24	protein_sigma_28
    
    # # deg mods -
    # # dM = readdlm("dM-Opt.dat")
    # # for value in dM
    # #     push!(pvec_initial, (1+0.10*rand())*value)
    # # end
    # push!(pvec_initial,0.05)      # 25	mRNA_clssrA 0.0085
	# push!(pvec_initial,0.05)	  # 26	mRNA_deGFPssrA 0.05
	# push!(pvec_initial,0.05)	  # 27	mRNA_sigma_28
	# push!(pvec_initial,30.0)	  # 28	protein_clssrA
	# push!(pvec_initial,34.0)	  # 29	protein_deGFPssrA
	# push!(pvec_initial,5.0)	      # 30	protein_sigma_28
	# push!(pvec_initial,1.0)	      # 31	protein_sigma_70

    # # start w/8.0 for TL capacity -
    # push!(pvec_initial,8.0*(1+0.25*randn()))

    # tmp -
    pvec_initial = parameter_array

    # check bounds -
    number_of_parameters = length(pvec_initial)
    for parameter_index = 1:number_of_parameters
        
        # what is the parameter value?
        p_i = pvec_initial[parameter_index]

        # is p_i outside of the bounds?
        lb_value = pvec_bounds[parameter_index,1]
        ub_value = pvec_bounds[parameter_index,2]
        
        if (p_i<lb_value)
            pvec_initial[parameter_index,1] = lb_value
        end

        if (p_i>ub_value)
            pvec_initial[parameter_index,1] = ub_value
        end
    end

    # search -
    options = Optim.Options(show_trace=true,show_every=100,iterations=1000)
    if algorithm == :SimulatedAnnealing
        result = optimize(OF,pvec_bounds[:,1],pvec_bounds[:,2],pvec_initial,Fminbox(SimulatedAnnealing()), options)
        return result
    elseif algorithm == :LBFGS
        options = Optim.Options(show_trace=true,show_every=20,iterations=200)
        result = optimize(OF,pvec_bounds[:,1],pvec_bounds[:,2],pvec_initial,Fminbox(BFGS()),options)
        return result
    end
end

# setup paths -
path_to_data_dir = "$(pwd())/data"
pV = readdlm("pvec_initial.dat")
pV = pV.*(1 .+ 0.15*randn(length(pV)))

# generate an ensemble of parameters -
number_of_p_sets = 10
number_of_parameters = 33 + 1
ensemble_array = zeros(number_of_parameters,number_of_p_sets)
for pset_index = 1:number_of_p_sets

    global pV

    algorithm_flag = :SimulatedAnnealing
    if mod(pset_index,2) == 0 
        algorithm_flag = :LBFGS
    end

    # call main -
    opt_result = main(path_to_data_dir, pV; algorithm = algorithm_flag)

    # grab the best pset, and then dump to disk -
    pV = Optim.minimizer(opt_result)
    ensemble_array[number_of_parameters,pset_index] = Optim.minimum(opt_result)
    for (p_index,p_value) in enumerate(pV)
        ensemble_array[p_index,pset_index] = p_value
    end

    # update -
    pV = pV.*(1 .+ 0.25*randn(length(pV)))
end

# dump -
writedlm("Ensemble-T9.dat",ensemble_array)


