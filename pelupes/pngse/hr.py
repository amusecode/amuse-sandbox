from amuse.units import units

from pngse import PNGSE, dimm_star

from amuse.datamodel import Particle

def simulate_evolution_tracks(
    stellar_evolution,
    masses=[0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0] | units.MSun,
    name_of_the_figure="HR.png"
):

    number_of_stars = len(masses)
    all_tracks_luminosity = []
    all_tracks_temperature = []
    all_tracks_stellar_type = []

    print(
        "The evolution across the Hertzsprung-Russell diagram of ",
        str(number_of_stars),
        " stars with\nvarying masses will be simulated..."
    )

    for j in range(number_of_stars):
        star = Particle()
        star.mass = masses[j]
        print("Created new star with mass: ", star.mass)

        star = stellar_evolution.particles.add_particle(star)
        stellar_evolution.commit_particles()

        luminosity_at_time = [] | units.LSun
        temperature_at_time = [] | units.K
        #~ stellar_type_at_time = [] | units.stellar_type

        stopped_evolving = False
        # Evolve this star until it changes into a compact stellar remnant
        # (white dwarf, neutron star, or black hole)
        while not dimm_star(star) and not stopped_evolving and stellar_evolution.model_time<20.|units.Gyr:
            luminosity_at_time.append(star.luminosity)
            temperature_at_time.append(star.temperature)
            #~ stellar_type_at_time.append(star.stellar_type)
            previous_age = star.age
            try:
                stellar_evolution.evolve_model(stellar_evolution.model_time+star.timestep)
                # Check whether the age has stopped increasing
                stopped_evolving = (star.age == previous_age)
            except Exception as ex:
                print(str(ex))
                stopped_evolving = True
            print(stellar_evolution.model_time)
        if stopped_evolving:
            print("Age did not increase during timestep. Aborted evolving...")
        else:
            #~ stellar_type_at_time.append(star.stellar_type)
            # Fudged: final stellar type annotation at previous (Teff, L);
            # BHs and neutron stars would otherwise fall off the chart.
            luminosity_at_time.append(luminosity_at_time[-1])
            temperature_at_time.append(temperature_at_time[-1])
        print(" ... evolved model to t = " + \
            str(star.age.as_quantity_in(units.Myr)))
        #~ print(
            #~ "Star has now become a: ",
            #~ star.stellar_type,
            #~ "(stellar_type: "
            #~ + str(
                #~ star.stellar_type.value_in(units.stellar_type)
            #~ )
            #~ + ")"
        #~ )
        print()
        all_tracks_luminosity.append(luminosity_at_time)
        all_tracks_temperature.append(temperature_at_time)
        #~ all_tracks_stellar_type.append(stellar_type_at_time)

#        Remove the star before creating the next one. See comments at the top.
        stellar_evolution.particles.remove_particle(star)
        stellar_evolution.reset()
        
    stellar_evolution.stop()

    plot_HR_diagram(
        masses,
        all_tracks_luminosity,
        all_tracks_temperature,
        all_tracks_stellar_type,
        name_of_the_figure
    )

    print("All done!")         


def plot_HR_diagram(
        masses,
        luminosity_tracks,
        temperature_tracks,
        stellar_type_tracks,
        plotfile):
        
        # This removes the need for ssh -X to be able to do plotting
        #~ import matplotlib
        #~ matplotlib.use("Agg")

        from matplotlib import pyplot
        print("Plotting the data...")

        pyplot.figure(figsize=(7, 8))
        pyplot.title('Hertzsprung-Russell diagram', fontsize=12)
        pyplot.xlabel('Effective Temperature (K)')
        pyplot.ylabel('Luminosity (solar luminosity)')

        # Define some strings for line formatting (colors, symbols, etc.), used
        # recurrently when many stars are simulated
        plot_format_strings_lines = ["r-", "y-", "c-", "b-", "m-"]
        len_fmt_str_lin = len(plot_format_strings_lines)
        plot_format_strings_symbols = [
            "r^", "y^", "c^", "b^", "m^", "rs", "ys", "cs", "bs", "ms"]
        len_fmt_str_sym = len(plot_format_strings_symbols)

        number_of_stars = len(masses)
        for j in range(number_of_stars):
            # Plot track of the current star j
            x_values = temperature_tracks[j].value_in(units.K)
            y_values = luminosity_tracks[j].value_in(units.LSun)
            pyplot.loglog(x_values, y_values,
                          plot_format_strings_lines[j % len_fmt_str_lin])
  
        pyplot.axis([300000., 2500., 1.e-2, 1.e6])
        # Or use these axes to also view neutron stars and black holes:
        # pyplot.axis([1.e7, 2500., 1.e-11, 1.e6])
        pyplot.savefig(plotfile)
        pyplot.show()

if __name__=="__main__":
    simulate_evolution_tracks(PNGSE(32),masses=[1,2,3,5,8] | units.MSun)
