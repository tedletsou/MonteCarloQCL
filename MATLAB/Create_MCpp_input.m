% Takes an input deck from Deck.m and creates input (.dat) file for C++ Monte Carlo code

%% Last updated: 2-10-21, Ted

%% Begin program 

% Call Deck([input fields])
s = Deck(0);

% Create blank input file to write to
fid = fopen('mcpp_input.dat', 'wt');

% Header for .dat file
fprintf(fid, '// Monte-Carlo code (C++) input deck \n');
fprintf(fid, '\n');

% Total number of modules
fprintf(fid, '// Total number of modules\n');
fprintf(fid, 'nummod %d\n', s.num_modules);
fprintf(fid, '\n');

% Barrier and well materials
fprintf(fid, '// Barrier and well materials\n');
fprintf(fid, 'barrier %s\n', s.barrier_mat);
fprintf(fid, 'well %s\n', s.well_mat);
fprintf(fid, '\n');

% Layer thicknesses in one module (in monolayers)
fprintf(fid, '// Layer thicknesses in one module (in monolayers)\n');
for i1 = 1:length(s.layer_thicknesses)
    fprintf(fid, 'laythick %d\t', s.layer_thicknesses(i1));
    fprintf(fid, '\n');
end
fprintf(fid, '\n');

% Doping values in one module (in cm^-3)
fprintf(fid, '// Doping values in one module (in cm^-3)\n');
for i1 = 1:length(s.doping_profile)
    fprintf(fid, 'laydop %d\n', s.doping_profile(i1));
end
fprintf(fid, '\n');

% Type of each layer (1: barrier, 2: well)
fprintf(fid, '// Type of each layer (1: barrier, 2: well)\n');
for i1 = 1:length(s.layer_types)
    fprintf(fid, 'laytype %d\n', s.layer_types(i1));
end
fprintf(fid, '\n');

% Total thickness of one module (in monolayers)
fprintf(fid, '// Total thickness of one module (in monolayers)\n');
fprintf(fid, 'modthick %d\n', s.module_thickness);
fprintf(fid, '\n');

% Applied fields (in V/m)
fprintf(fid, '// Applied fields (in V/m)\n');
for i1 = 1:length(s.applied_field)
    fprintf(fid, 'field_vals %d\n', s.applied_field(i1));
end
fprintf(fid, '\n');

% Grid meshing density
fprintf(fid, '// Grid meshing density\n');
fprintf(fid, 'meshden %d\n', s.nodes_per_monolayer);
fprintf(fid, '\n');

% Lattive constant (in angstroms)
fprintf(fid, '// Lattive constant (in angstroms)\n');
fprintf(fid, 'a_lat %d\n', s.lattice_constant);
fprintf(fid, '\n');

% Effective masses (relative to m0)
fprintf(fid, '// Effective masses (relative to m0)\n');
fprintf(fid, 'mstar_bar %d\n', s.effective_masses(1));
fprintf(fid, 'mstar_well %d\n', s.effective_masses(2));
fprintf(fid, '\n');

% Conduction band energies (relative to conduction band of well, eV)
fprintf(fid, '// Conduction band energies (relative to conduction band of well, eV)\n');
fprintf(fid, 'cband_bar %d\n', s.conduction_band_energies(1));
fprintf(fid, 'cband_well %d\n', s.conduction_band_energies(2));
fprintf(fid, '\n');

% Valence band energies (relative to conduction band well, eV)
fprintf(fid, '// Valence band energies (relative to conduction band well, eV)\n');
fprintf(fid, 'vband_bar %d\n', s.valence_band_energies(1));
fprintf(fid, 'vband_well %d\n', s.valence_band_energies(2));
fprintf(fid, '\n');

% Light hole band energies (relative to conduction band well, eV)
fprintf(fid, '// Light hole band energies (relative to conduction band well, eV)\n');
fprintf(fid, 'lhole_bar %d\n', s.light_hole(1));
fprintf(fid, 'lhole_well %d\n', s.light_hole(2));
fprintf(fid, '\n');

% Split-off band energies (relative to conduction band well, eV)
fprintf(fid, '// Split-off band energies (relative to conduction band well, eV)\n');
fprintf(fid, 'sploff_bar %d\n', s.split_off(1));
fprintf(fid, 'sploff_well %d\n', s.split_off(2));
fprintf(fid, '\n');

% Delta split-off energies, eV
fprintf(fid, '// Delta split-off energies, eV\n');
fprintf(fid, 'delso_bar %d\n', s.delta_so(1));
fprintf(fid, 'delso_well %d\n', s.delta_so(2));
fprintf(fid, '\n');

% Kane energies (derived from band edge, eV)
fprintf(fid, '// Kane energies (derived from band edge, eV)\n');
fprintf(fid, 'Ep_bar %d\n', s.Ep(1));
fprintf(fid, 'Ep_well %d\n', s.Ep(2));
fprintf(fid, '\n');

% Band gaps, eV
fprintf(fid, '// Band gaps, eV\n');
fprintf(fid, 'Eg_bar %d\n', s.energy_gap(1));
fprintf(fid, 'Eg_well %d\n', s.energy_gap(2));
fprintf(fid, '\n');

fclose(fid);
