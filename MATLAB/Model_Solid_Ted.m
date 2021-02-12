function s = Model_Solid_Ted(material, varargin)
% Manditory arguments:
% 'material': string of material
%
% Optional arguments
% 'mole_fraction_of_first_listed': number (default: 1)
% 'temperature': number (default: 300 K)
%
% Using parameters obtained from Material_Parameters.m, calculate absolute gamma point position,
% according to model-solid theoretic predictions of strain (on InP
% substrate). See Chuang p. 789 for calculation details.  
%
% Examples:
%	Model_Solid_Ted('AlInAs', 'mole_fraction_of_first_listed', 0.52, 'temperature', 250) would return the energy of Al_(0.52)In_(0.48)As at 250 K
%	Model_Solid_Ted('InAlAs', 'mole_fraction_of_first_listed', 0.53) would return the energy of Al_(0.47)In_(0.53)As at 300 K
%	Model_Solid_Ted('GaAs') would return the energy of GaAs at 300 K
%
% Notation:
%   x goes with FIRST atom 
%   'AlInAs' : Al(x)In(1 - x)As
%   'InAlAs' : In(x)Al(1 - x)As
%   'GaInAs' : Ga(x)In(1 - x)As
%   'InGaAs' : In(x)Ga(1 - x)As
%
% Material list is identical to that of Material_Parameters.m

%% Last updated: 6-8-20, Ted

%% Begin function

    % Defining default values
    mole_fraction_of_first_listed = 1;
    temperature = 300; % Kelvin

    % Checking for optional 
    for i1 = 1:2:length(varargin)
        switch varargin{i1}
            case 'mole_fraction_of_first_listed', mole_fraction_of_first_listed = varargin{i1 + 1};
            case 'temperature', temperature = varargin{i1 + 1};
        end
    end
    
    %% Calculating materials parameters
    
    matparams = Material_Parameters_Ted(material, 'mole_fraction_of_first_listed', mole_fraction_of_first_listed, ...
                                    'temperature', temperature); 
    InPparams = Material_Parameters_Ted('InP', 'mole_fraction_of_first_listed', 1, ...
                                    'temperature', temperature); 
    a0 = InPparams.lattice_constant;
    
    % Defining biaxial strains
    eps_par  = (a0 - matparams.lattice_constant) / a0;
    eps_perp = (-2 * matparams.C12 / matparams.C11) * eps_par;
    eps_tot  = 2 * eps_par + eps_perp;
    
    % Changes in band energies (Chuang p. 789)
    delta_Ec   = matparams.ac * eps_tot;
    delta_Evav = matparams.av * eps_tot;
    
    % Definition from Chuang
    Peps = -delta_Evav;
    Qeps = -matparams.b * (1 + 2 * (matparams.C12 / matparams.C11)) * eps_tot;
    Pc   = delta_Ec;
    
    % Changes in valence bands
    delta_EHH = matparams.delta_so / 3 - Qeps;
    delta_ELH = -matparams.delta_so / 6 + Qeps / 2  + (1 / 2) * sqrt(matparams.delta_so ^ 2 + ...
                                                                    2 * matparams.delta_so * Qeps + ...
                                                                    9 * Qeps ^ 2);
    delta_ESO = -matparams.delta_so / 6 + Qeps / 2  - (1 / 2) * sqrt(matparams.delta_so ^ 2 + ...
                                                                    2 * matparams.delta_so * Qeps + ...
                                                                    9 * Qeps ^ 2);
    
    % Defining final energy values
    final_Ec   = matparams.E_vav + matparams.delta_so / 3 + Peps + matparams.band_gap_gamma + Pc;
    final_Evav = matparams.E_vav + delta_Evav;
    final_EHH  = matparams.E_vav + delta_EHH;
    final_ELH  = matparams.E_vav + delta_ELH;
    final_ESO  = matparams.E_vav + delta_ESO;
    
    % Adding final energies to structure
    s      = matparams;
    s.Ec   = final_Ec;
    s.Evav = final_Evav;
    s.Elh  = final_ELH;
    s.Ehh  = final_EHH;
    s.Eso  = final_ESO;
    
end
    
     
    
    
    