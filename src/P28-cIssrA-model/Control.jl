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
# Function: calculate_transcription_control_array
# Description: Calculate the transcriptional control array at time t
# Generated on: 2019-10-02T16:15:21.68
#
# Input arguments:
# t::Float64 => Current time value (scalar)
# x::Array{Float64,1} => State array (number_of_species x 1)
# data_dictionary::Dict{String,Any} => Dictionary holding model parameters
#
# Output arguments:
# control_array::Array{Float64,1} => Transcriptional control array (number_of_genes x 1) at time t
# ----------------------------------------------------------------------------------- #
function calculate_transcription_control_array(t::Float64,x::Array{Float64,1},data_dictionary::Dict{String,Any})

	# initialize the control -
	control_array = zeros(4)

	# hack -
	x = abs.(x)

	# Alias the species -
	clssrA = x[1]
	deGFPssrA = x[2]
	sigma_28 = x[3]
	sigma_70 = x[4]
	mRNA_clssrA = x[5]
	mRNA_deGFPssrA = x[6]
	mRNA_sigma_28 = x[7]
	mRNA_sigma_70 = x[8]
	protein_clssrA = x[9]
	protein_deGFPssrA = x[10]
	protein_sigma_28 = x[11]
	protein_sigma_70 = x[12]

	# Alias the binding parameters -
	binding_parameter_dictionary = data_dictionary["binding_parameter_dictionary"]
	n_clssrA_sigma_28 = binding_parameter_dictionary["n_clssrA_sigma_28"]
	K_clssrA_sigma_28 = binding_parameter_dictionary["K_clssrA_sigma_28"]
	n_deGFPssrA_sigma_70 = binding_parameter_dictionary["n_deGFPssrA_sigma_70"]
	K_deGFPssrA_sigma_70 = binding_parameter_dictionary["K_deGFPssrA_sigma_70"]
	n_deGFPssrA_clssrA = binding_parameter_dictionary["n_deGFPssrA_clssrA"]
	K_deGFPssrA_clssrA = binding_parameter_dictionary["K_deGFPssrA_clssrA"]
	n_sigma_28_sigma_70 = binding_parameter_dictionary["n_sigma_28_sigma_70"]
	K_sigma_28_sigma_70 = binding_parameter_dictionary["K_sigma_28_sigma_70"]
	n_sigma_28_clssrA = binding_parameter_dictionary["n_sigma_28_clssrA"]
	K_sigma_28_clssrA = binding_parameter_dictionary["K_sigma_28_clssrA"]

	# Alias the control function parameters -
	control_parameter_dictionary = data_dictionary["control_parameter_dictionary"]
	W_clssrA_RNAP = control_parameter_dictionary["W_clssrA_RNAP"]
	W_clssrA_sigma_28 = control_parameter_dictionary["W_clssrA_sigma_28"]
	W_deGFPssrA_RNAP = control_parameter_dictionary["W_deGFPssrA_RNAP"]
	W_deGFPssrA_sigma_70 = control_parameter_dictionary["W_deGFPssrA_sigma_70"]
	W_deGFPssrA_clssrA = control_parameter_dictionary["W_deGFPssrA_clssrA"]
	W_sigma_28_RNAP = control_parameter_dictionary["W_sigma_28_RNAP"]
	W_sigma_28_sigma_70 = control_parameter_dictionary["W_sigma_28_sigma_70"]
	W_sigma_28_clssrA = control_parameter_dictionary["W_sigma_28_clssrA"]
	W_sigma_70_RNAP = control_parameter_dictionary["W_sigma_70_RNAP"]

	# Transfer function target:clssrA actor:sigma_28
	actor_set_clssrA_sigma_28 = [
		protein_sigma_28
	]
	actor = (prod(actor_set_clssrA_sigma_28))
	b_clssrA_sigma_28 = ((actor)^(n_clssrA_sigma_28))/((K_clssrA_sigma_28)^(n_clssrA_sigma_28)+(actor)^(n_clssrA_sigma_28))

	# Control function for clssrA -
	control_array[1] = (W_clssrA_RNAP+W_clssrA_sigma_28*b_clssrA_sigma_28)/(1+W_clssrA_RNAP+W_clssrA_sigma_28*b_clssrA_sigma_28)

	# Transfer function target:deGFPssrA actor:sigma_70
	actor_set_deGFPssrA_sigma_70 = [
		protein_sigma_70
	]
	actor = prod(actor_set_deGFPssrA_sigma_70)
	b_deGFPssrA_sigma_70 = ((actor)^(n_deGFPssrA_sigma_70))/((K_deGFPssrA_sigma_70)^(n_deGFPssrA_sigma_70)+(actor)^(n_deGFPssrA_sigma_70))

	# Transfer function target:deGFPssrA actor:clssrA
	actor_set_deGFPssrA_clssrA = [
		protein_clssrA
	]
	actor = prod(actor_set_deGFPssrA_clssrA)
	b_deGFPssrA_clssrA = (actor^(n_deGFPssrA_clssrA))/(K_deGFPssrA_clssrA^(n_deGFPssrA_clssrA)+actor^(n_deGFPssrA_clssrA))

	# Control function for deGFPssrA -
	control_array[2] = (W_deGFPssrA_RNAP+W_deGFPssrA_sigma_70*b_deGFPssrA_sigma_70)/(1+W_deGFPssrA_RNAP+W_deGFPssrA_sigma_70*b_deGFPssrA_sigma_70+W_deGFPssrA_clssrA*b_deGFPssrA_clssrA)

	# Transfer function target:sigma_28 actor:sigma_70
	actor_set_sigma_28_sigma_70 = [
		protein_sigma_70
	]
	actor = prod(actor_set_sigma_28_sigma_70)
	b_sigma_28_sigma_70 = (actor^(n_sigma_28_sigma_70))/(K_sigma_28_sigma_70^(n_sigma_28_sigma_70)+actor^(n_sigma_28_sigma_70))

	# Transfer function target:sigma_28 actor:clssrA
	actor_set_sigma_28_clssrA = [
		protein_clssrA
	]
	actor = prod(actor_set_sigma_28_clssrA)
	b_sigma_28_clssrA = (actor^(n_sigma_28_clssrA))/(K_sigma_28_clssrA^(n_sigma_28_clssrA)+actor^(n_sigma_28_clssrA))

	# Control function for sigma_28 -
	control_array[3] = (W_sigma_28_RNAP+W_sigma_28_sigma_70*b_sigma_28_sigma_70)/(1+W_sigma_28_RNAP+W_sigma_28_sigma_70*b_sigma_28_sigma_70+W_sigma_28_clssrA*b_sigma_28_clssrA)

	# Control function for sigma_70 -
	control_array[4] = (W_sigma_70_RNAP)/(1+W_sigma_70_RNAP)

	# return -
	return control_array
end

#
# ----------------------------------------------------------------------------------- #
# Function: calculate_translation_control_array
# Description: Calculate the translation control array at time t
# Generated on: 2019-10-02T16:15:21.684
#
# Input arguments:
# t::Float64 => Current time value (scalar)
# x::Array{Float64,1} => State array (number_of_species x 1)
# data_dictionary::Dict{String,Any} => Dictionary holding model parameters
#
# Output arguments:
# control_array::Array{Float64,1} => Translation control array (number_of_genes x 1) at time t
# ----------------------------------------------------------------------------------- #
function calculate_translation_control_array(t::Float64,x::Array{Float64,1},data_dictionary::Dict{String,Any})

	# initialize the control -
	control_array = ones(4)

	# for S70 - no translation
	control_array[4] = 0.0

	# correct for "resource"?
    correction_term = (x[13]/100.0)
    control_array = control_array*correction_term

	# return -
	return control_array
end
