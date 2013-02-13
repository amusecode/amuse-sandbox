package amuse.code;

import java.nio.*;

public interface CodeInterface {

    int set_color(IntBuffer index_of_the_particle, DoubleBuffer red, DoubleBuffer green, DoubleBuffer blue, DoubleBuffer alpha, int npoints);
   
    //first parameter is output parameter!
    int new_particle(IntBuffer index_of_the_particle, IntBuffer type, DoubleBuffer x, DoubleBuffer y, DoubleBuffer z, DoubleBuffer radius, DoubleBuffer red, DoubleBuffer green, DoubleBuffer blue, DoubleBuffer alpha, int npoints);

    int set_type(IntBuffer index_of_the_particle, IntBuffer type, int npoints);

    int set_radius(IntBuffer index_of_the_particle, DoubleBuffer radius, int npoints);

    int store_view(double time);

    int cleanup_code();

    int recommit_parameters();

    int initialize_code();

    int set_position(IntBuffer index_of_the_particle, DoubleBuffer x, DoubleBuffer y, DoubleBuffer z, int npoints);

    int commit_parameters();

    int commit_particles();
}