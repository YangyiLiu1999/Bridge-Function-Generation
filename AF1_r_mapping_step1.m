function [RHO_Z2_map, rho_V2_map, rho_Z1_map] = r_function(copula_data)
%========================================================================================
% INPUT:
%----------------------------------------------------------------------------------------
    % copula_para: a sturcture including distribution type distribution parameters 
    %              which is available to specific copula function types,
    %              including: 'Gaussian', 't', 'Frank', 'Clayton',
    %              'Gumbel';
    % rho_Z1_ana: the correlation coefficient vector of main variable
    % rho_Z2_ana: the correlation coefficient vector of subordinate variable
    % eta:        the distance interval vector corresponding to rho_Z1_ana and rho_Z2_ana
%-----------------------------------------------------------------------------------------
% OUTPUT:
%-----------------------------------------------------------------------------------------
    % rho_V2_ana: the correlation coefficient vector of processing variable
%-----------------------------------------------------------------------------------------
%% STEP 1: Construction of the inverse function of r
%% STEP 1-0: the scale of the coefficient of correlation
rho_V2_map = 0.999: -0.01: -0.999;     % variable in process
rho_Z1_map = 0.999: -0.01: -0.999;     % main variable

d_V2         = 0.1;
d_Z1         = 0.1;
v2_1         = -3: d_V2: 3;
v2_2         = -3: d_V2: 3;
z1_1         = -3: d_Z1: 3;
z1_2         = -3: d_Z1: 3;
[V2_1, V2_2] = meshgrid(v2_1, v2_2);
[Z1_1, Z1_2] = meshgrid(z1_1, z1_2);
%% STEP 1-1: evaluate the quadruple integral
for ii = 1: length(rho_V2_map)
    for jj = 1: length(rho_Z1_map)
        P_V2 = 1./(2*pi .* sqrt(1 - rho_V2_map(ii).^2)) .*...
                exp( -1./2 .* (V2_1.^2 - 2 .* rho_V2_map(ii) .* V2_1 .* V2_2 + V2_2.^2 )./ (1 - rho_V2_map(ii).^2));        % the pdf of variable in process 
        P_Z1 = 1./(2*pi .* sqrt(1 - rho_Z1_map(jj).^2)) .*...
                exp( -1./2 .* (Z1_1.^2 - 2 .* rho_Z1_map(jj) .* Z1_1 .* Z1_2 + Z1_2.^2 )./ (1 - rho_Z1_map(jj).^2));        % the pdf of subordibate variable
        for num1 =  1: length(v2_1)       
                part1(num1,:) = norminv( copula_transform(normcdf(z1_1),normcdf(v2_1(num1)), copula_data) );  
        end
        for num2 = 1: length(v2_2)
                part2(num2,:) = norminv( copula_transform(normcdf(z1_2),normcdf(v2_2(num2)), copula_data) );   
        end
        for num1 = 1: length(v2_1)
            for num2 = 1: length(v2_2)
                for num3 = 1: length(z1_1)
                    for num4 = 1: length(z1_2)
                        rho_Z2_ini(num1, num2, num3, num4) = part1(num1, num3) .* part2(num2, num4) .* P_V2(num2, num1) .* P_Z1(num4, num3); 
                    end
                end
            end
        end
%% STEP 1-2: the coefficient of correlation of Z2 (subordunate variable)
        RHO_Z2_map(ii,jj) = sum( sum( sum(sum(rho_Z2_ini)))) .* d_V2.^2 .* d_Z1.^2;
    end
end
end % function

function [U2] = copula_transform(U1, z, copula_data)
%===========================================================
% INPUT:
%-----------------------------------------------------------
    % U1:         the vector of mian variable cdf
    % z:          the value of process variable cdf
    % distr_para: a sturcture including distribution type
    %             distribution parameters
%-----------------------------------------------------------
% OUTPUT:
%-----------------------------------------------------------
    % U2:        the vector of subordinate variable cdf
%===========================================================
    % basic parameters
    copula_type = copula_data.type;
    copula_para = copula_data.para;   
    % 
    if sum(strcmpi(copula_type,{'Gaussian'})) >= 1
        rho = copula_para;
        U2 = normcdf( rho * norminv(U1) +  norminv(z) * sqrt(1 - rho^2) );
    elseif sum(strcmpi(copula_type,{'t'})) >= 1
        rho = copula_para(1); nu = copula_para(2);
        U2 = tcdf( tinv( z, nu+1) .* sqrt(( nu + (tinv(U1, nu)).^2 )...
            .* ( 1 - rho.^2) ./ (nu+1)) + rho .* tinv(U1,nu+1), nu);
    elseif sum(strcmpi(copula_type,{'Frank'})) >= 1
        theta = copula_para;
        U2 = -1 ./theta .* log((exp( -theta * U1) .*(1-z) +...
            z.* exp(-theta)) ./ (exp( -theta * U1) .*(1-z) + z));
    elseif sum(strcmpi(copula_type,{'Clayton'})) >= 1
        theta = copula_para;
        U2 = U1 .*( z.^(-theta./(1+theta)) + U1.^theta - 1).^(-1./theta);   
    elseif sum(strcmpi(copula_type,{'Gumbel'})) >= 1
        theta = copula_para;
         for ii  = 1: length(U1)
          u1 = U1(ii);
             myf = @(x) ( -log(u1)) .^(theta - 1) .*...
             exp( -( (-log(u1)).^theta + (-log(x)).^theta) .^ (1 ./theta)) ./...
            (u1 .*(  (-log(u1)).^theta + (-log(x)).^theta ).^ (1 - 1./theta))...
             - z;
             U2(ii) = fzero(myf, [0,1]);
         end % for
    end % if
end % function