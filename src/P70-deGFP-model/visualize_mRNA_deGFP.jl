include("Include.jl")

function plot_mRNA_deGFP_data(time_array, simulation_data_array, experimental_data_dictionary)

    (NR,NC) = size(time_array)    
    for index = 1:NC
        plot(time_array[:,index],1000.0*simulation_data_array[:,index],color="#7AD7F0",alpha=0.80,lw=0.5)
    end

    # Plot the region -
    # calculate the mean and std -
    μ = mean(1000.0*simulation_data_array,dims=2)
    σ = std(1000.0*simulation_data_array,dims=2)
    LB = μ .- (1.96/sqrt(1))*σ
    UB = μ .+ (1.96/sqrt(1))*σ
    fill_between(time_array[:,1], vec(UB), vec(LB), color="#d9d9d9")

    # show std error region -
    LB = μ .- (1.96/sqrt(NC))*σ
    UB = μ .+ (1.96/sqrt(NC))*σ
    fill_between(time_array[:,1], vec(UB), vec(LB), color="#fe9929")

    # Plot mean -
    plot(time_array[:,1],μ,"--",color="#fff7bc",lw=1.5)

    # plot the experimemtal data -
    TEXP = experimental_data_dictionary["mRNA_data_array"][:,1]
    DATA = experimental_data_dictionary["mRNA_data_array"][:,2]
    STD = (1.96/sqrt(3))*experimental_data_dictionary["mRNA_data_array"][:,3]
    yerr_array = transpose([STD STD])
    errorbar(TEXP, DATA,yerr=yerr_array,fmt="o",mfc="#252525",mec="#252525",color="#252525", lw=1)

    # labels -
    xlabel("Time (hr)",fontsize=16)
    ylabel("Concentration (nM)",fontsize=16)

     # set the face color -
     fig = gcf()
     ax = gca()
     ax.set_facecolor("#F5FCFF")
     fig.set_facecolor("#ffffff")
 
    # wwrite -
    savefig("./plots/mRNA-deGFP-Ensemble.pdf")
end


function main(path_to_simulation_dir::String, path_to_plot_file::String)

    # what index is mRNA deGFP?
    state_index = 3

    # load the experimemtal data -
    exp_dd = load_experimental_data_dictionary("./data")

    # how many files?
    file_name_array = searchdir(path_to_simulation_dir, ".dat")
    number_of_trials = length(file_name_array)

    # initialize -> hardcode the dimension for now 
    data_array = zeros(1601,number_of_trials)
    time_array = zeros(1601,number_of_trials)

    # read the simulation dir -
    for (file_index,file_name) in enumerate(file_name_array)
        
        # load -
        file_path = "$(path_to_simulation_dir)/$(file_name)"
        sim_data_array = readdlm(file_path)

        # what size?
        (NR,NC) = size(sim_data_array)
        for step_index = 1:NR
            data_array[step_index,file_index] = sim_data_array[step_index,(state_index+1)]
            time_array[step_index,file_index] = sim_data_array[step_index,1]
        end
    end

    # plot -
    plot_mRNA_deGFP_data(time_array, data_array, exp_dd)
end

path_to_simulation_dir = "./simulations"
path_to_plot_file = "./plots"
main(path_to_simulation_dir, path_to_plot_file)