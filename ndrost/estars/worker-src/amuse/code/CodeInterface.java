package amuse.code;

public interface CodeInterface {

    int set_color(int[] index_of_the_particle, double[] red, double[] green, double[] blue, double[] alpha, int npoints);
   
    //first parameter is output parameter!
    int new_particle(int[] index_of_the_particle, int[] type, double[] x, double[] y, double[] z, double[] radius, double[] red, double[] green, double[] blue, double[] alpha, int npoints);

    int set_type(int[] index_of_the_particle, int[] type, int npoints);

    int set_radius(int[] index_of_the_particle, double[] radius, int npoints);

    int store_view(double time);

    int cleanup_code();

    int recommit_parameters();

    int initialize_code();

    int set_position(int[] index_of_the_particle, double[] x, double[] y, double[] z, int npoints);

    int commit_parameters();

    int commit_particles();
}
