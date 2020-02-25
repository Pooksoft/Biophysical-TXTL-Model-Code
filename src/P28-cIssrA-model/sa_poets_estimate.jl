include("Include.jl")

# Evaluates the objective function values -
function local_refienment_step(path_to_data_dir, parameter_array; sigma=0.05, iteration_max=100)
  
    # initialize -
    number_of_parameters = length(parameter_array)
    BIG = 1e10

    # load the default data_dictionary -
    time_start = 0.0
    time_stop = 16.0
    time_step_size = 0.01

    # what is the host_type?
    host_type = :cell_free

    # path to parameters -
    path_to_biophysical_constants_file = "./CellFree.json"

    # wght array -
    W = diagm(ones(5))
    W[3,3] = 0.1
    W[4,4] = 1000.0
    W[5,5] = 1.0

    # Load the data dictionary (uses the default biophysical_constants file)
    default_data_dictionary = build_data_dictionary((time_start,time_stop,time_step_size), path_to_biophysical_constants_file, host_type)

    # Set some values for the weights, gene lengths etc -
    model_data_dictionary = customize_data_dictionary(default_data_dictionary,host_type)

    # load the experimental data -
    exp_data_dictionary = load_experimental_data_dictionary(path_to_data_dir)

    # setup the functions -
    OF(P) = objective_function(P,time_start,time_step_size,time_stop,model_data_dictionary,exp_data_dictionary)

    function _compute_error_total(objective_array,W)
        value = transpose(objective_array)*W*objective_array 
        return value[1]
    end

    # calculate the starting error -
    parameter_array_best = parameter_array
    error_array = BIG*ones(4)
    error_array[1] = _compute_error_total(OF(parameter_array_best), W)

    # main refinement loop -
    iteration_counter = 1
    while (iteration_counter<iteration_max)

        # take a step up -
        parameter_up = parameter_array_best.*(1 .+ sigma*rand(number_of_parameters))
        parameter_up = check_parameter_bounds(parameter_up)

        # take a step down -
        parameter_down = parameter_array_best.*(1 .- sigma*rand(number_of_parameters))
        parameter_down = check_parameter_bounds(parameter_down)

        # Evaluate the obj function -
        error_array[2] = _compute_error_total(OF(parameter_up),W)
        error_array[3] = _compute_error_total(OF(parameter_down),W)

        # Calculate a correction factor -
        a = error_array[2]+error_array[3] - 2.0*error_array[1]
        parameter_corrected = parameter_array_best
        if (a>0.0)
            amda = -0.5*(error_array[3] - error_array[2])/a
            parameter_corrected = parameter_array_best .+ amda*rand(number_of_parameters)
            parameter_corrected = check_parameter_bounds(parameter_corrected)
            error_array[4] = _compute_error_total(OF(parameter_corrected), W)
        end

        # Which step has the min error?
        min_index = argmin(error_array)
        if (min_index == 1)
            parameter_array_best = parameter_array_best
        elseif (min_index == 2)
            parameter_array_best = parameter_up
        elseif (min_index == 3)
            parameter_array_best = parameter_down
        elseif (min_index == 4)
            parameter_array_best = parameter_corrected
        end

        # Update the local error
        error_array[1] = error_array[min_index]

        @show iteration_counter,error_array[min_index]

        # update local counter -
        iteration_counter = iteration_counter + 1
    end

    return parameter_array_best
end
  
function check_parameter_bounds(parameter_array)

    # setup paramter bounds -
    pvec_bounds = [

        # dG's -
        39995.69565886529   68282.2841495398    ;   # 1 W_clssrA_RNAP
        -50000.0  -1.0   ;   # 2 W_clssrA_sigma_28
        39995.69565886529   68282.2841495398    ;   # 3 W_deGFPssrA_RNAP
        -35689.208748865043  -24993.98186199351    ;   # 4 W_deGFPssrA_sigma_70
        -50000.0  20000.0   ;   # 5 W_deGFPssrA_clssrA
        39995.69565886529   68282.2841495398    ;   # 6 W_sigma_28_RNAP
        -50000.0  25000.0   ;   # 7 W_sigma_28_sigma_70
        -50000.0  25000.0   ;   # 8 W_sigma_28_clssrA

        # binding parameters -
        0.5 10.0            ;   # 9 n_clssrA_sigma_28
        0.001 100.0         ;   # 10 K_clssrA_sigma_28
        0.5037004212911795  2.2       ;   # 11 n_deGFPssrA_sigma_70
        25.539812443045424 100.0         ;   # 12 K_deGFPssrA_sigma_70
        0.5 10.0            ;   # 13 n_deGFPssrA_clssrA
        0.001 100.0         ;   # 14 K_deGFPssrA_clssrA
        0.5 10.0            ;   # 15 n_sigma_28_sigma_70
        0.001 100.0         ;   # 16 K_sigma_28_sigma_70
        0.5 10.0            ;   # 17 n_sigma_28_clssrA
        0.001 100.0         ;   # 18 K_sigma_28_clssrA

        # time constants -
		0.0001 100.0         ;   # 19	mRNA_clssrA
		0.0039273045649573235 7.17667775643900         ;	# 20	mRNA_deGFPssrA
		0.0001 100.0         ;	# 21	mRNA_sigma_28
		0.0001 100.0         ;	# 22	protein_clssrA
		0.35035678021027006 10.203484322453043         ;	# 23	protein_deGFPssrA
		0.001 100.0	        ;	# 24	protein_sigma_28
        
        # degradation mods -
        0.1 10.0            ;	# 25	mRNA_clssrA 0.0085
		0.0011598053735671 3.6284714410194088 	     ;	# 26	mRNA_deGFPssrA 0.05
		0.1 10.0 	        ;	# 27	mRNA_sigma_28
		5.0 10.0 	        ;	# 28	protein_clssrA
		40.0 60.0 	        ;	# 29	protein_deGFPssrA
		0.13 10.0 	        ;	# 30	protein_sigma_28
        0.132047561240967526 10.137258626014857 	        ;	# 31	protein_sigma_70

        # w -
        4.0 10.0            ;   # 32    translation capacity half-life

        # KL -
        10.81768297340375 396.3396750321431         ;   # 33    KL in muM
    ];

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

    # return -
    return pvec_initial
end

function neighbor_function(parameter_array; sigma=0.05)

    # setup -
    number_of_parameters = length(parameter_array)
    
    # calculate new parameter array -
    new_parameter_array = parameter_array.*(1 .+ sigma*randn(number_of_parameters))

    # check the bounds and return -
    return check_parameter_bounds(new_parameter_array)
end

function cooling_function(temperature)

  # define my new temperature -
  alpha = 0.9
  return alpha*temperature
end

function acceptance_probability_function(rank_array,temperature)
    return (exp(-rank_array[end]/temperature))
end

function objective_function(parameter_guess_array,time_start,time_step_size,time_stop,model_data_dictionary,exp_data_dictionary)

    # what is the host_type?
    host_type = :cell_free

    # some global parameters -
    BIG = 1e10
    SMALL = 1e-6

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
    error_term_array = BIG*ones(5,1)
    
    # mRNA GFP -
    #WM1 = diagm(1.0 ./((mRNA_GFP_std_exp).^2))
    error_vector_1 = (mRNA_GFP_exp .- mRNA_GFP_sim)
    error_term_array[1,1] = 1.0*transpose(error_vector_1)*error_vector_1

    # mRNA sigma28 -
    #WM2 = diagm(1.0 ./((mRNA_sigma28_std_exp).^2))
    error_vector_2 = (mRNA_sigma28_exp .- mRNA_sigma28_sim)
    error_term_array[2,1] = 1.0*transpose(error_vector_2)*error_vector_2

    # mRNA C1 -
    #WM3 = diagm( 1.0 ./((mRNA_C1_std_exp).^2))
    error_vector_3 = (mRNA_C1_exp .- mRNA_C1_sim)
    error_term_array[3,1] = 1.0*transpose(error_vector_3)*error_vector_3

    # protein deGFP -
    #WM4 = diagm( 1.0 ./((mRNA_C1_std_exp).^2) )
    error_vector_4 = (protein_GFP_exp .- protein_GFP_sim)
    error_term_array[4,1] = 1.0*transpose(error_vector_4)*error_vector_4

    # the fifth term is an upper bound on C1 -
    prot_C1 = XSIM[:,9]
    UB = 100.0
    MV_C1 = maximum(prot_C1)
    error_term_C1 = 100000*max(0,(MV_C1 - UB))
    error_term_array[5,1] = error_term_C1

    # # 6th term - limit sigma28 mRNA peak to 200 -
    # mRNA_sigma28_actual = (1000)*XSIM[:,7]
    # UB_sigma28 = 800.0
    # MV_sigma28 = maximum(mRNA_sigma28_actual)
    # error_term_mRNA_sigma28_overshoot = 100000*max(0, (MV_sigma28 - UB_sigma28))
    # error_term_array[6,1] = error_term_mRNA_sigma28_overshoot
    # ===================================================================================================== #

    # error total -
    #error_total = sum(error_term_array)

    # return -
    return error_term_array
end

function main(path_to_data_dir::String, initial_parameter_array::Array{Float64,1}; rank_cutoff::Int64=4, maximum_number_of_iterations::Int64=100)

    # load the default data_dictionary -
    time_start = 0.0
    time_stop = 16.0
    time_step_size = 0.01

    # what is the host_type?
    host_type = :cell_free
    number_of_parameters = length(initial_parameter_array)

    # path to parameters -
    path_to_biophysical_constants_file = "./CellFree.json"

    # Load the data dictionary (uses the default biophysical_constants file)
    default_data_dictionary = build_data_dictionary((time_start,time_stop,time_step_size), path_to_biophysical_constants_file, host_type)

    # Set some values for the weights, gene lengths etc -
    model_data_dictionary = customize_data_dictionary(default_data_dictionary,host_type)

    # load the experimental data -
    exp_data_dictionary = load_experimental_data_dictionary(path_to_data_dir)

    # setup the functions -
    OF(P) = objective_function(P,time_start,time_step_size,time_stop,model_data_dictionary,exp_data_dictionary)
    NF(P) = neighbor_function(P;sigma=0.01) 
    
    # make call to POETs -
    (EC,PC,RA) = estimate_ensemble(OF,NF,acceptance_probability_function,cooling_function,initial_parameter_array;rank_cutoff=rank_cutoff,maximum_number_of_iterations=maximum_number_of_iterations)    

    # return -
    return (EC,PC,RA)
end

# setup -
path_to_data_dir = "$(pwd())/data"
pV = readdlm("pvec_initial_v5.dat")
pV = neighbor_function(pV; sigma=0.25)
EC = 0
PC = 0
RA = 0

# execute -
number_of_trials = 20
for trial_index = 1:number_of_trials

    global pV
    global EC
    global PC
    global RA

   # do a local step -
    if (mod(trial_index,2) == 0)
        
        # find the lowest score pV -
        sum_error_array = sum(EC,dims=1)
        best_p_index = argmin(vec(sum_error_array))
        pV_best = PC[:,best_p_index]

        # local refine -
        pV = local_refienment_step(path_to_data_dir, pV_best; iteration_max=500) 
    end

    # main -
    (EC,PC,RA) = main(path_to_data_dir, vec(pV); rank_cutoff=3,maximum_number_of_iterations=100)

    # dump results to disk -
    fname = "./poets_ensemble_c1_restriction/RA_T$(trial_index).dat"
    writedlm(fname,RA)
    fname = "./poets_ensemble_c1_restriction/EC_T$(trial_index).dat"
    writedlm(fname,EC)
    fname = "./poets_ensemble_c1_restriction/PC_T$(trial_index).dat"
    writedlm(fname,PC)

    # find next (this will be overwritten by local search when index is even) -
    ECS = sum(EC,dims=1)
    idx_min = argmin(vec(ECS))
    pV = PC[:,idx_min]
end