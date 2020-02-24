# ----------------------------------------------------------------------------------- #
# Copyright (c) 2019 Varnerlab
# Robert Frederick Smith School of Chemical and Biomolecular Engineering
# Cornell University, Ithaca NY 14850
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ----------------------------------------------------------------------------------- #
#
# ----------------------------------------------------------------------------------- #
# Function: DataDictionary
# Description: Holds simulation and model parameters as key => value pairs in a Julia Dict()
# Generated on: 2019-10-02T16:15:21.453
#
# Input arguments:
# time_start::Float64 => Simulation start time value (scalar)
# time_stop::Float64 => Simulation stop time value (scalar)
# time_step::Float64 => Simulation time step (scalar)
#
# Output arguments:
# data_dictionary::Dict{String,Any} => Dictionary holding model and simulation parameters as key => value pairs
# ----------------------------------------------------------------------------------- #
function build_data_dictionary(time_span::Tuple{Float64,Float64,Float64}, path_to_biophysical_constants_file::String = "./Default.json", host_type::Symbol = :bacteria)::Dict{String,Any}

	# load the biophysical_constants dictionary
	biophysical_constants_dictionary = build_biophysical_dictionary(path_to_biophysical_constants_file, host_type)

	# stoichiometric_matrix and dilution_matrix -
	stoichiometric_matrix = readdlm("./Network.dat")

	# number of states, and rates -
	(number_of_states,number_of_rates) = size(stoichiometric_matrix)

	# array of species types -
	species_symbol_type_array = [
		:gene	;	# 1	clssrA
		:gene	;	# 2	deGFPssrA
		:gene	;	# 3	sigma_28
		:gene	;	# 4	sigma_70
		:mrna	;	# 5	mRNA_clssrA
		:mrna	;	# 6	mRNA_deGFPssrA
		:mrna	;	# 7	mRNA_sigma_28
		:mrna	;	# 8	mRNA_sigma_70
		:protein	;	# 9	protein_clssrA
		:protein	;	# 10	protein_deGFPssrA
		:protein	;	# 11	protein_sigma_28
		:protein	;	# 12	protein_sigma_70
	]

	# for some whacky reason, we need to add the species_symbol_type_array to the biophysical dictionary -
	biophysical_constants_dictionary["species_symbol_type_array"] = species_symbol_type_array

	# array of gene lengths -
	gene_coding_length_array = [
		1000.0	;	# 1	clssrA
		1000.0	;	# 2	deGFPssrA
		1000.0	;	# 3	sigma_28
		1000.0	;	# 4	sigma_70
	]

	# array of mRNA coding lengths -
	mRNA_coding_length_array = [
		gene_coding_length_array[1]	;	# 5	1	mRNA_clssrA
		gene_coding_length_array[2]	;	# 6	2	mRNA_deGFPssrA
		gene_coding_length_array[3]	;	# 7	3	mRNA_sigma_28
		gene_coding_length_array[4]	;	# 8	4	mRNA_sigma_70
	]

	# array of mRNA coding lengths -
	protein_coding_length_array = [
		round((0.33)*mRNA_coding_length_array[1])	;	# 9	1	protein_clssrA
		round((0.33)*mRNA_coding_length_array[2])	;	# 10	2	protein_deGFPssrA
		round((0.33)*mRNA_coding_length_array[3])	;	# 11	3	protein_sigma_28
		round((0.33)*mRNA_coding_length_array[4])	;	# 12	4	protein_sigma_70
	]

	# array of gene concentrations -
	gene_abundance_array = [
		5.0	;	# (nM) 1	clssrA
		5.0	;	# (nM) 2	deGFPssrA
		5.0	;	# (nM) 3	sigma_28
		0.0	;	# (nM) 4	sigma_70	# we have no sigma 70
	]

	# initial condition array -
	initial_condition_array = [
		gene_abundance_array[1]		;	# 1	clssrA
		gene_abundance_array[2]		;	# 2	deGFPssrA
		gene_abundance_array[3]		;	# 3	sigma_28
		gene_abundance_array[4]		;	# 4	sigma_70
		0.0							;	# 5	mRNA_clssrA
		0.0							;	# 6	mRNA_deGFPssrA
		0.0							;	# 7	mRNA_sigma_28
		0.0							;	# 8	mRNA_sigma_70
		0.0							;	# 9	protein_clssrA
		0.0							;	# 10	protein_deGFPssrA
		0.0							;	# 11	protein_sigma_28
		0.0							;	# 12	protein_sigma_70

		# extra state -
		100.0						;	# 13
	]

	binding_parameter_dictionary = Dict{String,Float64}()
	binding_parameter_dictionary["n_clssrA_sigma_28"] = 1.0
	binding_parameter_dictionary["K_clssrA_sigma_28"] = 0.05
	binding_parameter_dictionary["n_deGFPssrA_sigma_70"] = 1.0
	binding_parameter_dictionary["K_deGFPssrA_sigma_70"] = 0.05
	binding_parameter_dictionary["n_deGFPssrA_clssrA"] = 1.0
	binding_parameter_dictionary["K_deGFPssrA_clssrA"] = 0.05
	binding_parameter_dictionary["n_sigma_28_sigma_70"] = 1.0
	binding_parameter_dictionary["K_sigma_28_sigma_70"] = 0.05
	binding_parameter_dictionary["n_sigma_28_clssrA"] = 1.0
	binding_parameter_dictionary["K_sigma_28_clssrA"] = 0.05

	# Alias the control function parameters -
	control_parameter_dictionary = Dict{String,Float64}()
	control_parameter_dictionary["W_clssrA_RNAP"] = 0.001
	control_parameter_dictionary["W_clssrA_sigma_28"] = 1.0
	control_parameter_dictionary["W_deGFPssrA_RNAP"] = 0.001
	control_parameter_dictionary["W_deGFPssrA_sigma_70"] = 1.0
	control_parameter_dictionary["W_deGFPssrA_clssrA"] = 1.0
	control_parameter_dictionary["W_sigma_28_RNAP"] = 0.001
	control_parameter_dictionary["W_sigma_28_sigma_70"] = 1.0
	control_parameter_dictionary["W_sigma_28_clssrA"] = 1.0
	control_parameter_dictionary["W_sigma_70_RNAP"] = 0.0

	# degradation modifiers -
	degradation_modifier_array = [
		0.0	;	# 1	clssrA
		0.0	;	# 2	deGFPssrA
		0.0	;	# 3	sigma_28
		0.0	;	# 4	sigma_70
		1.0	;	# 5	mRNA_clssrA
		1.0	;	# 6	mRNA_deGFPssrA
		1.0	;	# 7	mRNA_sigma_28
		1.0	;	# 8	mRNA_sigma_70
		1.0	;	# 9	protein_clssrA
		1.0	;	# 10	protein_deGFPssrA
		1.0	;	# 11	protein_sigma_28
		1.0	;	# 12	protein_sigma_70
	]

	# time constant modifiers -
	time_constant_modifier_array = [
		0.0	;	# 1	clssrA
		0.0	;	# 2	deGFPssrA
		0.0	;	# 3	sigma_28
		0.0	;	# 4	sigma_70
		1.0	;	# 5	mRNA_clssrA
		1.0	;	# 6	mRNA_deGFPssrA
		1.0	;	# 7	mRNA_sigma_28
		1.0	;	# 8	mRNA_sigma_70
		1.0	;	# 9	protein_clssrA
		1.0	;	# 10	protein_deGFPssrA
		1.0	;	# 11	protein_sigma_28
		1.0	;	# 12	protein_sigma_70
	]

	# Dilution degrdation matrix -
	dilution_degradation_matrix = build_dilution_degradation_matrix(biophysical_constants_dictionary,species_symbol_type_array,degradation_modifier_array)

	# Precompute the translation parameters -
	translation_parameter_array = precompute_translation_parameter_array(biophysical_constants_dictionary, protein_coding_length_array, time_constant_modifier_array,host_type)

	# Precompute the kinetic limit of transcription -
	transcription_kinetic_limit_array = precompute_transcription_kinetic_limit_array(biophysical_constants_dictionary, gene_coding_length_array, gene_abundance_array, time_constant_modifier_array, host_type)

	# Parameter name index array -
	parameter_name_mapping_array = [
		"n_clssrA_sigma_28"	;	# 1
		"K_clssrA_sigma_28"	;	# 2
		"n_deGFPssrA_sigma_70"	;	# 3
		"K_deGFPssrA_sigma_70"	;	# 4
		"n_deGFPssrA_clssrA"	;	# 5
		"K_deGFPssrA_clssrA"	;	# 6
		"n_sigma_28_sigma_70"	;	# 7
		"K_sigma_28_sigma_70"	;	# 8
		"n_sigma_28_clssrA"	;	# 9
		"K_sigma_28_clssrA"	;	# 10
		"W_clssrA_RNAP"	;	# 11
		"W_clssrA_sigma_28"	;	# 12
		"W_deGFPssrA_RNAP"	;	# 13
		"W_deGFPssrA_sigma_70"	;	# 14
		"W_deGFPssrA_clssrA"	;	# 15
		"W_sigma_28_RNAP"	;	# 16
		"W_sigma_28_sigma_70"	;	# 17
		"W_sigma_28_clssrA"	;	# 18
		"W_sigma_70_RNAP"	;	# 19
		"rnapII_concentration"	;	# 20
		"ribosome_concentration"	;	# 21
		"degradation_constant_mRNA"	;	# 22
		"degradation_constant_protein"	;	# 23
		"kcat_transcription"	;	# 24
		"kcat_translation"	;	# 25
		"maximum_specific_growth_rate"	;	# 26
		"saturation_constant_transcription"	;	# 27
		"saturation_constant_translation"	;	# 28
	]

	# =============================== DO NOT EDIT BELOW THIS LINE ============================== #
	data_dictionary = Dict{String,Any}()
	data_dictionary["number_of_states"] = number_of_states
	data_dictionary["species_symbol_type_array"] = species_symbol_type_array
	data_dictionary["initial_condition_array"] = initial_condition_array
	data_dictionary["gene_coding_length_array"] = gene_coding_length_array
	data_dictionary["mRNA_coding_length_array"] = mRNA_coding_length_array
	data_dictionary["protein_coding_length_array"] = protein_coding_length_array
	data_dictionary["stoichiometric_matrix"] = stoichiometric_matrix
	data_dictionary["dilution_degradation_matrix"] = dilution_degradation_matrix
	data_dictionary["binding_parameter_dictionary"] = binding_parameter_dictionary
	data_dictionary["control_parameter_dictionary"] = control_parameter_dictionary
	data_dictionary["parameter_name_mapping_array"] = parameter_name_mapping_array
	data_dictionary["transcription_kinetic_limit_array"] = transcription_kinetic_limit_array
	data_dictionary["translation_parameter_array"] = translation_parameter_array
	data_dictionary["degradation_modifier_array"] = degradation_modifier_array
	data_dictionary["time_constant_modifier_array"] = time_constant_modifier_array
	data_dictionary["biophysical_constants_dictionary"] = biophysical_constants_dictionary
	data_dictionary["half_life_translation_capacity"] = 8.0

	# extra stuff -
	data_dictionary["R"] = 8.314 			# J mol^-1 K^-1
	data_dictionary["T_K"] = 273.15 + 29.0 	# K
	# =============================== DO NOT EDIT ABOVE THIS LINE ============================== #
	return data_dictionary
end
