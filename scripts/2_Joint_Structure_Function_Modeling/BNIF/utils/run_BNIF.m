function [f, P] = run_BNIF(subject,fmri_dir,tract_dir,out_dir)
    
    % Argument Parsing
    %p = inputParser;
    %addParameter(p, 'subject', @isstring);
    %addParameter(p, 'fmri_dir', @isstring);
    %addParameter(p, 'tract_dir', @isstring);
    %addParameter(p, 'out_dir', @isstring);
    %parse(p,subject,fmri_dir,tract_dir,out_dir);

    %subject = p.Results.subject;
    %fmri_dir = p.Results.fmri_dir;
    %tract_dir = p.Results.tract_dir;
    %out_dir = p.Results.out_dir;
    %contrast=int2str(contrast);

    addpath('/util/opt/IBM-ILOG-CPLEX/12.10/intel/19.0.1/cplex/matlab/x86-64_linux')

    % Define file paths
    connectome_fn = tract_dir + "/" + subject + "/" + "MMP1_structural_connectome_mu_scaled.csv";
    func_matrix_fn = fmri_dir + "/" + "R_matrix_sum.csv";

    % Read connectome data
    connectome = readmatrix(connectome_fn);
    connectome = rescale(connectome);

    % Read functional matrix data
    func_matrix = readmatrix(func_matrix_fn);
    func_matrix = rescale(func_matrix)
    
    c=parcluster('local');
    c.NumWorkers=16;
    %saveProfile(c);
    parpool(c,16);

    [f, P] = BNIF_optimized(connectome,func_matrix,'rho', 0.02);
    f_sum = sum(f, 3);
    f_max = max(f, [], 3);

    f_out = out_dir + "/" + "MMP1_structure_function_networks.mat";
    f_sum_out = out_dir + "/" + "MMP1_summed_structure_function_network.csv";
    f_max_out = out_dir + "/" + "MMP1_maxed_structure_function_network.csv";
    p_out = out_dir + "/" + "MMP1_underestimation_correction_matrix.csv";

    save(f_out,'f');
    writematrix(P,p_out);
    writematrix(f_sum,f_sum_out);
    writematrix(f_max,f_max_out);

end

