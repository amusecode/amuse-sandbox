
int get_mass(int * index_of_the_particle, double * mass, int length);

int commit_particles();

int set_unit_length(double code_length_unit);

int get_tree_domain_update_frequency(double * tree_domain_update_frequency);

int get_time(double * time);

int get_epsilon_dm_part(int * index_of_the_particle, double * radius, int length);

int set_nsmooth(int nsmooth);

int set_mass(int * index_of_the_particle, double * mass, int length);

int get_info_file(char * * info_file);

int get_time_limit_cpu(double * time_limit_cpu);

int get_index_of_first_particle(int * index_of_the_particle);

int get_total_radius(double * radius);

int get_potential_at_point(double eps, double x, double y, double z, double * phi);

int get_periodic_boundaries_flag(int * value);

int get_softening_gas_max_phys(double * softening_gas_max_phys);

int set_min_gas_hsmooth_fractional(double min_gas_hsmooth_fractional);

int set_periodic_boundaries_flag(int value);

int get_gdgtol(double * gdgtol);

int get_gadget_output_directory(char * * gadget_output_directory);

int get_n_neighbours(int * index_of_the_particle, double * num_neighbours, int length);

int get_interpret_heat_as_feedback_flag(int * value);

int get_total_mass(double * mass);

int evolve_model(double time);

int set_state_sph(int * index_of_the_particle, double * mass, double * x, double * y, double * z, double * vx, double * vy, double * vz, double * u, int length);

int get_nsmooth(int * nsmooth);

int get_epsgas(double * gas_epsilon);

int set_type_of_timestep_criterion(int type_of_timestep_criterion);

int get_redshift_begin(double * redshift_begin);

int get_omega_lambda(double * omega_lambda);

int get_unit_velocity(double * code_velocity_unit);

int set_eps2(double epsilon_squared);

int set_time_limit_cpu(double time_limit_cpu);

int get_begin_time(double * time);

int get_eps2(double * epsilon_squared);

int get_epsilon(double * epsilon);

int get_index_of_next_particle(int index_of_the_particle, int * index_of_the_next_particle);

int get_alpha_visc(int * index_of_the_particle, double * alpha, int length);

int new_sph_particle(int * index_of_the_particle, double mass, double x, double y, double z, double vx, double vy, double vz, double u);

int delete_particle(int index_of_the_particle);

int set_omega_lambda(double omega_lambda);

int set_gadget_output_directory(char * gadget_output_directory);

int get_isotherm(int * isothermal_flag);

int set_interpret_heat_as_feedback_flag(int value);

int set_max_size_timestep(double max_size_timestep);

int get_dtalpha_visc(int * index_of_the_particle, double * dtalpha, int length);

int get_unit_time(double * code_time_unit);

int get_box_size(double * value);

int get_unit_length(double * code_length_unit);

int set_omega_zero(double omega_zero);

int get_potential(int * index_of_the_particle, double * Potential, int length);

int set_gdgop(int gadget_cell_opening_flag);

int get_unit_mass(double * code_mass_unit);

int synchronize_model();

int set_epsgas(double gas_epsilon);

int set_state(int * index_of_the_particle, double * mass, double * x, double * y, double * z, double * vx, double * vy, double * vz, int length);

int set_epsilon(double epsilon);

int get_state(int * index_of_the_particle, double * mass, double * x, double * y, double * z, double * vx, double * vy, double * vz, int length);

int set_min_size_timestep(double min_size_timestep);

int set_redshift_begin(double redshift_begin);

int get_gdgop(int * gadget_cell_opening_flag);

int get_time_step(double * time_step);

int recommit_particles();

int set_box_size(double value);

int set_gdgtol(double gdgtol);

int get_kinetic_energy(double * kinetic_energy);

int get_number_of_particles(int * number_of_particles);

int get_redshift_max(double * redshift_max);

int set_internal_energy(int * index_of_the_particle, double * u, int length);

int get_hubble_param(double * hubble_param);

int set_alpha(double alpha);

int set_unit_time(double code_time_unit);

int get_time_between_statistics(double * time_between_statistics);

int get_internal_energy(int * index_of_the_particle, double * u, int length);

int get_interpret_kicks_as_feedback_flag(int * value);

int set_unit_mass(double code_mass_unit);

int set_acceleration(int index_of_the_particle, double ax, double ay, double az);

int get_center_of_mass_position(double * x, double * y, double * z);

int set_time_step(double time_step);

int set_err_tol_int_accuracy(double err_tol_int_accuracy);

int get_alpha(double * alpha);

int get_comoving_integration_flag(int * comoving_integration_flag);

int get_hydro_state_at_point(double x, double y, double z, double vx, double vy, double vz, double * rho, double * rhovx, double * rhovy, double * rhovz, double * rhoe);

int set_courant(double courant);

int get_timings_file(char * * timings_file);

int get_center_of_mass_velocity(double * vx, double * vy, double * vz);

int get_radius(int index_of_the_particle, double * radius);

int get_bh_tol(double * bh_tol);

int get_omega_zero(double * omega_zero);

int get_omega_baryon(double * omega_baryon);

int get_energy_file(char * * energy_file);

int get_courant(double * courant);

int get_cpu_file(char * * cpu_file);

int get_smoothing_length(int * index_of_the_particle, double * h_smooth, int length);

int get_d_internal_energy_dt(int * index_of_the_particle, double * du_dt, int length);

int get_density(int * index_of_the_particle, double * rho, int length);

int set_begin_time(double time);

int set_nsmtol(double n_neighbour_tol);

int set_cpu_file(char * cpu_file);

int get_epsilon_gas_part(int * index_of_the_particle, double * radius, int length);

int evolve_to_redshift(double redshift);

int new_dm_particle(int * index_of_the_particle, double mass, double x, double y, double z, double vx, double vy, double vz);

int set_min_gas_temp(double min_gas_temp);

int get_err_tol_int_accuracy(double * err_tol_int_accuracy);

int get_softening_halo_max_phys(double * softening_halo_max_phys);

int get_gamma(double * gamma);

int set_energy_file(char * energy_file);

int get_type_of_timestep_criterion(int * type_of_timestep_criterion);

int set_radius(int index_of_the_particle, double radius);

int set_bh_tol(double bh_tol);

int set_softening_halo_max_phys(double softening_halo_max_phys);

int set_interpret_kicks_as_feedback_flag(int value);

int cleanup_code();

int get_min_size_timestep(double * min_size_timestep);

int get_thermal_energy(double * thermal_energy);

int set_comoving_integration_flag(int comoving_integration_flag);

int recommit_parameters();

int initialize_code();

int get_redshift(double * redshift);

int set_tree_domain_update_frequency(double tree_domain_update_frequency);

int get_pressure(int * index_of_the_particle, double * pressure, int length);

int get_potential_energy(double * potential_energy);

int get_min_gas_temp(double * min_gas_temp);

int get_state_sph(int * index_of_the_particle, double * mass, double * x, double * y, double * z, double * vx, double * vy, double * vz, double * u, int length);

int get_nogravity(int * no_gravity_flag);

int get_gravity_at_point(double eps, double x, double y, double z, double * ax, double * ay, double * az);

int set_redshift_max(double redshift_max);

int get_velocity(int * index_of_the_particle, double * vx, double * vy, double * vz, int length);

int get_min_gas_hsmooth_fractional(double * min_gas_hsmooth_fractional);

int set_hubble_param(double hubble_param);

int set_info_file(char * info_file);

int set_softening_gas_max_phys(double softening_gas_max_phys);

int get_nsmtol(double * n_neighbour_tol);

int get_time_max(double * time_max);

int get_max_size_timestep(double * max_size_timestep);

int get_position(int * index_of_the_particle, double * x, double * y, double * z, int length);

int set_position(int * index_of_the_particle, double * x, double * y, double * z, int length);

int set_time_max(double time_max);

int set_time_between_statistics(double time_between_statistics);

int get_eps_is_h(int * eps_is_h_flag);

int get_acceleration(int * index_of_the_particle, double * ax, double * ay, double * az, int length);

int set_unit_velocity(double code_velocity_unit);

int commit_parameters();

int set_omega_baryon(double omega_baryon);

int set_timings_file(char * timings_file);

int set_velocity(int * index_of_the_particle, double * vx, double * vy, double * vz, int length);

