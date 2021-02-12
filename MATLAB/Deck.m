function s = Deck(field_in)
% Manditory arguments:
%   'field in': Applied field in V/m
%
% Input deck for QCL solver

%% Last updated: 1-30-21, Ted

%% Begin function

    % Layer type definitions. Odd (1) indicates barrier.  Even (2) indicates
    % well.
    b = 1;
    w = 2;

    % Lattice constant for InP
    InPparams = Material_Parameters_Ted('InP'); 
    a0 = InPparams.lattice_constant;

    % Lattice constant for active region
    sb = Model_Solid_Ted('InAlAs', 'mole_fraction_of_first_listed', 0.52, 'temperature', 300);
    b0 = sb.lattice_constant;    
        
    %% Structure definitions 
    
    % Band structure inputs and parameters in terms of (1 / 2) * lattice constant (monolayers)
    s.num_modules          = 2;
    s.layer_thicknesses    = [38, 21, 5, 60, 6, 55, 7, 45, 10, 39, 11, 38, 17, 37, 24, 35] / (b0 / 2); %% 11.5 um
    s.doping_profile       = 1 * [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.21e17, 2.21e17, 0];
    s.doping_type = 'block';
%     s.layer_thicknesses    = [200, 100, 200] / (b0 / 2); 
%     s.doping_profile       = [0, 0, 0];
    %s.layer_thicknesses    = [37, 31, 27, 75, 7, 58, 15, 52, 18, 41, 15, 38, 16, 35, 17, 34, 20, 34, 23, 34, 28, 33] / (b0 / 2);
    %s.layer_thicknesses    = [40, 19, 8, 56, 10, 51, 11, 42, 13, 32, 15, 32, 20, 31, 29, 30] / (b0 / 2); %% 9.4 um
    %s.layer_thicknesses    = [39, 22, 8, 60, 9, 59, 10, 52, 13, 43, 14, 38, 15, 36, 16, 34, 19, 33, 23, 32, 25, 32, 29, 31] / (b0 / 2);
    %s.layer_thickness      = [40, 15, 10, 48, 12, 47, 13, 42, 15, 32, 17, 30, 18, 28, 23, 26, 34, 24] / (b0 / 2); %% 7.3 um
    %s.doping_profile       = 1 * [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.53e17, 2.53e17, 0];
    
    % Defining layers
    s.layer_types          = s.layer_thicknesses;
    s.layer_types(1:2:end) = b;
    s.layer_types(2:2:end) = w;
    
    % Extra layers left (monolayers)
    s.extra_layers_left      = [];
    s.extra_layers_left_type = [];
    
    % Extra layers right (monolayers)
    s.extra_layers_right      = [] / (b0 / 2);
    s.extra_layers_right_type = []; 
    
    s.extra_energy = 0;
    
    s.module_thickness = sum(s.layer_thicknesses);

    % Simulation parametes
    s.applied_field       = field_in; % V/m
    s.nodes_per_monolayer = 10;      % Unitless
   
    %% Material definitions 
    
    % Material parameters
    sw = Model_Solid_Ted('InGaAs', 'mole_fraction_of_first_listed', 0.53, 'temperature', 300);
    sb = Model_Solid_Ted('InAlAs', 'mole_fraction_of_first_listed', 0.52, 'temperature', 300);
    sw_text = "('InGaAs', 'mole_fraction_of_first_listed', 0.53, 'temperature', 300)";
    sb_text = "('InAlAs', 'mole_fraction_of_first_listed', 0.52, 'temperature', 300)";
    
    % Material constants
    s.lattice_constant = sb.lattice_constant;                    % Angstrom
    s.effective_masses = [sb.eff_mass_gamma, sw.eff_mass_gamma]; % unitless
    
    % Material energies (relative to well conduction band)
    s.conduction_band_energies = [sb.Ec, sw.Ec] - sw.Ec; 
    % Valence band from (3.2.53) Jerome's book 
    s.valence_band_energies    = (2 / 3) * [sb.Elh, sw.Elh] + (1 / 3) * [sb.Eso, sw.Eso] - sw.Ec; 
    
    % Individual valence bands for non-parabolicity
    s.light_hole               = [sb.Elh, sw.Elh] - sw.Ec;
    s.split_off                = [sb.Eso, sw.Eso] - sw.Ec;
    s.delta_so                 = [sb.delta_so, sw.delta_so];

    % Kane energies
    s.Ep                       = [sb.E_p, sw.E_p];

    % CBO override : to match literature values
    % Jerome's group, Hugi et al (2009): estimate 519 meV for In_0.52Al_0.48As / In_0.53Ga_0.47As
    % Capasso's group, Pflugl et al (2010): states 520 meV
    CBO = 0.5195; % eV
    MS_CBO = abs(diff(s.conduction_band_energies));
    %CBO = MS_CBO * 0.70;
    
    s.conduction_band_energies(1) = s.conduction_band_energies(1) - MS_CBO + CBO;
    s.valence_band_energies(1)    = s.valence_band_energies(1) - MS_CBO;  
    s.energy_gap = [sb.band_gap_gamma, sw.band_gap_gamma];
    
    % Kane energies calculated using Dave's masters thesis and page 125 of
    % Chuang.  With the approximation made in Chuang (me << m0), this
    % reduces to Dave's theis.  However, I am not making that
    % approximation.
    
    % Dave's equation
    % s.Ep = (3 ./ s.effective_masses) .* (2 ./ s.energy_gap + (1 ./ (s.energy_gap + s.delta_so))) .^ (-1);
    
    % Chuang reduced
    % s.Ep = ((1 ./ s.effective_masses)  .* (s.energy_gap .* (s.energy_gap + s.delta_so) ./ (s.energy_gap + (2 / 3) * s.delta_so)));
    
    % Chuang full
    s.Ep = ((1 ./ s.effective_masses)  .* (1 - s.effective_masses) .* (s.energy_gap .* (s.energy_gap + s.delta_so) ./ (s.energy_gap + (2 / 3) * s.delta_so)));
    
    s.barrier = b;
    s.well    = w;
    
    s.barrier_mat = sb_text;
    s.well_mat    = sw_text;
end
    
    
    

    
    