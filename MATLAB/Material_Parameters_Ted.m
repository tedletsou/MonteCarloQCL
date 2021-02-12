function output = Material_Parameters_Ted(material, varargin)
% Manditory arguments:
% 'material': string of material
%
% Optional arguments
% 'mole_fraction_of_first_listed': number (default: 1)
% 'temperature': number (default: 300 K)
%
% Using linear interpolation and a quadratic bowing parameter, calculates band parameters.
%
% Unless otherwise specified, all parameters are from:
%   I. Vurgaftman, J. R. Meyer, and L.-R. Ram-Mohan, J. Appl. Phys. 89, 5815 (2001).
%
% Examples:
%	Material_Parameters('AlInAs', 'mole_fraction_of_first_listed', 0.52, 'temperature', 250) would return the parameters of Al_(0.52)In_(0.48)As at 250 K
%	Material_Parameters('InAlAs', 'mole_fraction_of_first_listed', 0.52) would return the parameters of Al_(0.48)In_(0.52)As at 300 K
%	Material_Parameters('GaAs', 'mole_fraction_of_first_listed', 1, 'temperature', 77) would return the parameters of GaAs at 77 K
%	Material_Parameters('InP') would return the parameters of InP at 300 K
%
% Material list:
%   GaAs
%   AlAs
%   InAs
%   InP
%   InGaAs
%   AlInAs
%
% Notation:
%   x goes with FIRST atom 
%   'AlInAs' : Al(x)In(1 - x)As
%   'InAlAs' : In(x)Al(1 - x)As
%   'GaInAs' : Ga(x)In(1 - x)As
%   'InGaAs' : In(x)Ga(1 - x)As

%% Last updated: 6-15-20, Ted

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

    %% BINARY COMPOUNDS
    
    %% Band parameters of GaAs, stored in struct

    GaAs.lattice_constant = 5.65325 + 3.88e-5 * (temperature - 300); % Angstroms
    GaAs.band_gap_gamma   = 1.519;                                   % eV
    GaAs.band_gap_X       = 1.981;                                   % eV
    GaAs.band_gap_L       = 1.815;                                   % eV
    GaAs.delta_so         = 0.341;                                   % eV
    GaAs.eff_mass_gamma   = 0.067;                                   % relative mass (unitless)
    GaAs.ac               = -7.17;                                   % eV
    GaAs.av               = -1.16;                                   % eV
    GaAs.b                = -2.0;                                    % eV
    GaAs.C11              = 1221;                                    % GPa
    GaAs.C12              = 566;                                     % GPa
    GaAs.E_vav            = -6.92;                                   % eV (Chuang)
    GaAs.E_p              = 28.8;                                    % eV
    GaAs.E_p              = 22.71;                                   % eV (Jerome override)
    GaAs.E_p              = 20.02;                                   % eV (An ab initio based approach to optical properties of semiconductor heterostructures)
    
    %% Band parameters of AlAs, stored in struct
    
    AlAs.lattice_constant = 5.6611 + 2.90e-5 * (temperature - 300); % Angstroms
    AlAs.band_gap_gamma   = 3.099;                                  % eV
    AlAs.band_gap_X       = 2.27;                                   % eV
    AlAs.band_gap_L       = 2.46;                                   % eV
    AlAs.delta_so         = 0.28;                                   % eV
    AlAs.eff_mass_gamma   = 0.15;                                   % relative mass (unitless)
    AlAs.ac               = -5.64;                                  % eV
    AlAs.av               = -2.47;                                  % eV
    AlAs.b                = -2.3;                                   % eV
    AlAs.C11              = 1250;                                   % GPa
    AlAs.C12              = 534;                                    % GPa
    AlAs.E_vav            = -7.49;                                  % eV (Chuang)
    AlAs.E_p              = 21.1;                                   % eV
    
    %% Band parameters of InAs, stored in struct
    
    InAs.lattice_constant = 6.0583 + 2.74e-5 * (temperature - 300);	% Angstroms
    InAs.band_gap_gamma   = 0.417;                                  % eV
    InAs.band_gap_X       = 1.433;                                  % eV
    InAs.band_gap_L       = 1.133;                                  % eV
    InAs.delta_so         = 0.39;                                   % eV
    InAs.eff_mass_gamma   = 0.026;                                  % relative mass (unitless)
    InAs.ac               = -5.08;                                  % eV
    InAs.av               = -1.0;                                   % eV
    InAs.b                = -1.8;                                   % eV
    InAs.C11              = 832.9;                                  % GPa
    InAs.C12              = 452.6;                                  % GPa
    InAs.E_vav            = -6.67;                                  % eV (Chuang)
    InAs.E_p              = 21.5;                                   % eV
    InAs.E_p              = 21.1;                                   % eV (Jerome override)

    %% Band parameters of InP, stored in struct  
    
    InP.lattice_constant = 5.8697 + 2.79e-5 * (temperature - 300); % Angstroms
    InP.band_gap_gamma   = 1.4236;                                 % eV
    InP.band_gap_X       = 2.384 - (3.7e-4 * temperature);         % eV, not a typo
    InP.band_gap_L       = 2.014;                                  % eV
    InP.delta_so         = 0.108;                                  % eV
    InP.eff_mass_gamma   = 0.0795;                                 % relative mass (unitless)
    InP.ac               = -6.0;                                   % eV
    InP.av               = -0.6;                                   % eV
    InP.b                = -2.0;                                   % eV
    InP.C11              = 1011;                                   % GPa
    InP.C12              = 561;                                    % GPa
    InP.E_vav            = -7.04;                                  % eV (Chuang)
    InP.E_p              = 20.7;                                   % eV
    
    %% TERNARY ALLOYS
    % Ternary bowing parameters (bowing term is defined as -C * x * (1 - x)
    
    %% Bowing band parameters of InGaAs, stored in struct      

    InGaAs.lattice_constant = 0;       % Angstroms
    InGaAs.band_gap_gamma   = 0.477;   % eV
    InGaAs.band_gap_X       = 1.4;     % eV
    InGaAs.band_gap_L       = 0.33;    % eV
    InGaAs.delta_so         = 0.15;    % eV
    InGaAs.eff_mass_gamma   = 0.0091;  % relative mass (unitless)
    InGaAs.ac               = 2.61;    % eV
    InGaAs.av               = 0;       % eV
    InGaAs.b                = 0;       % eV
    InGaAs.C11              = 0;       % GPa
    InGaAs.C12              = 0;       % GPa
    InGaAs.E_vav            = 0;       % eV
    InGaAs.E_p              = -1.48;   % eV

    %% Bowing band parameters of AlInAs, stored in struct 
    
    AlInAs.lattice_constant = 0;     % Angstroms
    AlInAs.band_gap_gamma   = 0.70;  % eV
    AlInAs.band_gap_X       = 0;     % eV
    AlInAs.band_gap_L       = 0;     % eV
    AlInAs.delta_so         = 0.15;  % eV
    AlInAs.eff_mass_gamma   = 0.049; % relative mass (unitless)
    AlInAs.ac               = -1.4;  % eV
    AlInAs.av               = 0;     % eV
    AlInAs.b                = 0;     % eV
    AlInAs.C11              = 0;     % GPa
    AlInAs.C12              = 0;     % GPa
    AlInAs.E_vav            = 0;     % eV
    AlInAs.E_p              = -4.81; % eV
    
    %% Bowing band parameters of AlGaAs, stored in struct 
    
    AlGaAs.lattice_constant = 0;     % Angstroms
    AlGaAs.band_gap_gamma   = -1.27 + 1.310 * mole_fraction_of_first_listed;  % eV
    AlGaAs.band_gap_X       = 0.055;     % eV
    AlGaAs.band_gap_L       = 0;     % eV
    AlGaAs.delta_so         = 0;  % eV
    AlGaAs.eff_mass_gamma   = 0;     % relative mass (unitless)
    AlGaAs.ac               = 0;     % eV
    AlGaAs.av               = 0;     % eV
    AlGaAs.b                = 0;     % eV
    AlGaAs.C11              = 0;     % GPa
    AlGaAs.C12              = 0;     % GPa
    AlGaAs.E_vav            = 0;     % eV
    AlGaAs.E_p              = 0;     % eV
    
    %% Create vectors containing all material parameters:
    
    GaAs_params   = [GaAs.lattice_constant, GaAs.band_gap_gamma, GaAs.band_gap_X,  ... 
                     GaAs.band_gap_L, GaAs.delta_so, GaAs.eff_mass_gamma, GaAs.ac, ...
                     GaAs.av, GaAs.b, GaAs.C11, GaAs.C12, GaAs.E_vav, GaAs.E_p];
    AlAs_params   = [AlAs.lattice_constant, AlAs.band_gap_gamma, AlAs.band_gap_X, ...
                     AlAs.band_gap_L, AlAs.delta_so, AlAs.eff_mass_gamma, AlAs.ac, ...
                     AlAs.av, AlAs.b, AlAs.C11, AlAs.C12, AlAs.E_vav, AlAs.E_p];
    InAs_params   = [InAs.lattice_constant, InAs.band_gap_gamma, InAs.band_gap_X, ...
                     InAs.band_gap_L, InAs.delta_so, InAs.eff_mass_gamma, InAs.ac, ...
                     InAs.av, InAs.b, InAs.C11, InAs.C12, InAs.E_vav, InAs.E_p];
    InP_params    = [InP.lattice_constant, InP.band_gap_gamma, InP.band_gap_X, ...
                     InP.band_gap_L, InP.delta_so, InP.eff_mass_gamma, InP.ac, ...
                     InP.av, InP.b, InP.C11, InP.C12, InP.E_vav, InP.E_p];
    InGaAs_params = [InGaAs.lattice_constant, InGaAs.band_gap_gamma, InGaAs.band_gap_X, ...
                     InGaAs.band_gap_L, InGaAs.delta_so, InGaAs.eff_mass_gamma, InGaAs.ac, ...
                     InGaAs.av, InGaAs.b, InGaAs.C11, InGaAs.C12, InGaAs.E_vav, InGaAs.E_p];
    AlInAs_params = [AlInAs.lattice_constant, AlInAs.band_gap_gamma, AlInAs.band_gap_X, ...
                     AlInAs.band_gap_L, AlInAs.delta_so, AlInAs.eff_mass_gamma, AlInAs.ac, ...
                     AlInAs.av, AlInAs.b, AlInAs.C11, AlInAs.C12, AlInAs.E_vav, AlInAs.E_p];  
    AlGaAs_params = [AlGaAs.lattice_constant, AlGaAs.band_gap_gamma, AlGaAs.band_gap_X, ...
                     AlGaAs.band_gap_L, AlGaAs.delta_so, AlGaAs.eff_mass_gamma, AlGaAs.ac, ...
                     AlGaAs.av, AlGaAs.b, AlGaAs.C11, AlGaAs.C12, AlGaAs.E_vav, AlGaAs.E_p];                

    %% Interpolation
        
    % Inline interpolater linearly interpolates (or extrapolates) between two points (x1, y1) and (x2, y2)
    function interp = InlineInterpolater(x1, y1, x2, y2, xval)
        interp = ((xval - x1)./(x2 - x1)) .* y2 + (1-((xval - x1) ./ (x2 - x1))) .* y1;
    end

    const = mole_fraction_of_first_listed * (1 - mole_fraction_of_first_listed);
    
    % Perform interpolation
    if (strcmp(material, 'InGaAs') == 1)
        output_params = InlineInterpolater(0, GaAs_params, 1, InAs_params, mole_fraction_of_first_listed) - const * InGaAs_params;
    elseif (strcmp(material, 'GaInAs') == 1)
        output_params = InlineInterpolater(0, InAs_params, 1, GaAs_params, mole_fraction_of_first_listed) - const * InGaAs_params;
    elseif (strcmp(material, 'InAlAs') == 1)
        output_params = InlineInterpolater(0, AlAs_params, 1, InAs_params, mole_fraction_of_first_listed) - const * AlInAs_params;
    elseif (strcmp(material, 'AlInAs') == 1)
        output_params = InlineInterpolater(0, InAs_params, 1, AlAs_params, mole_fraction_of_first_listed) - const * AlInAs_params;
    elseif (strcmp(material, 'AlGaAs') == 1)
        output_params = InlineInterpolater(0, GaAs_params, 1, AlAs_params, mole_fraction_of_first_listed) - const * AlGaAs_params;
    elseif (strcmp(material, 'GaAlAs') == 1)
        output_params = InlineInterpolater(0, AlAs_params, 1, GaAs_params, mole_fraction_of_first_listed) - const * AlGaAs_params;
    elseif (strcmp(material, 'GaAs') == 1)
        output_params = GaAs_params;
    elseif (strcmp(material, 'AlAs') == 1)
        output_params = AlAs_params;
    elseif (strcmp(material, 'InAs') == 1)
        output_params = InAs_params;
    elseif (strcmp(material, 'InP') == 1)
        output_params = InP_params;
    else
        fprintf('Material selection invalid!\n')
        return
    end

    %% Output
    
    % Create output structure from vector
    clear output

    output.lattice_constant = output_params(1);  % Angstroms
    output.band_gap_gamma   = output_params(2);  % eV
    output.band_gap_X       = output_params(3);  % eV
    output.band_gap_L       = output_params(4);  % eV
    output.delta_so         = output_params(5);  % eV
    output.eff_mass_gamma   = output_params(6);  % relative mass (unitless)
    output.ac               = output_params(7);  % eV
    output.av               = output_params(8);  % eV
    output.b                = output_params(9);  % eV
    output.C11              = output_params(10); % GPa
    output.C12              = output_params(11); % GPa
    output.E_vav            = output_params(12); % eV
    output.E_p              = output_params(13); % eV

end