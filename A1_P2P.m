function [R_U2] = BridgeFunctionGenerate(R_V0, R_W0, copula_data, eta) 
%% 
% input: 
% 1. R_V0: double vector (n.*1): The correlation function of the first variables
% 2. R_W0: double vector (n.*1): The correlation function of the second variables
% 3. coupla_data: structure: copula_data_type = 'name of the copula', such
% as Gaussian, t, Gumbel;
% 4. copula.para: double vector or scalar;
% output: the bridge function corresponding to the two variables and the
% copula type
% 5. eta: double vector (n .*1): the sequence of the spatial interval 
% output:
% 1. R_U2: the bridge function;
copula_data.type = 'Gaussian';
copula_data.para = 0.5;
[RHO_Z2_map, rho_V2_map, rho_Z1_map] = AF1_r_mapping_step1(copula_data);
% RHO_Z2_map: the mapping function for interpalation
% rho_Z1_map: every coefficient value of the first variable for loop
% rho_V2_map: every coefficient value of the bridge function for loop
%% Generate the bridge function
R_V  = rho_Z1_map;
corr = RHO_Z2_map;
R_0  = rho_w0;
for ii = 1:length(eta) 
    for kk = 1:length( corr(:,1) )
        corr_W0(kk) = interp1(R_V, corr(kk,:), R_V0(ii), 'linear','extrap');
    end 
         W_max = max(corr_W0); W_min = min(corr_W0);

         if R_W0(ii) > W_max
             corr_U(ii) = interp1(corr_W0,  R_0 ,  R_W0(ii) , 'nearest' ,'extrap');
         else if R_W0(ii) < W_min
              corr_U(ii) = interp1(corr_W0,  R_0 ,  R_W0(ii) , 'nearest' ,'extrap');
             else
                  corr_U(ii) = interp1( corr_W0,  R_0 , R_W0(ii) , 'linear', 'extrap');
             end
          end
end
R_U2 = corr_U;
end