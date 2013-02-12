package nl.esciencecenter.estars;

import java.nio.DoubleBuffer;
import java.nio.IntBuffer;

import amuse.code.CodeInterface;

public class Code implements CodeInterface {

    @Override
    public int set_color(IntBuffer index_of_the_particle, DoubleBuffer red, DoubleBuffer green, DoubleBuffer blue,
            int npoints) {
        System.err.println("set_color");
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public int new_particle(IntBuffer index_of_the_particle, IntBuffer type, DoubleBuffer x, DoubleBuffer y,
            DoubleBuffer z, DoubleBuffer red, DoubleBuffer green, DoubleBuffer blue, int npoints) {
        System.err.println("new_particle");
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public int set_type(IntBuffer index_of_the_particle, IntBuffer type, int npoints) {
        System.err.println("set_type");
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public int store_view(double time) {
        System.err.println("store view @ " + time);
        return 0;
    }

    @Override
    public int cleanup_code() {
        System.err.println("cleanup_code!");
        return 0;
    }

    @Override
    public int recommit_parameters() {
        System.err.println("recommit_particles!");
        return 0;
    }

    @Override
    public int initialize_code() {
        System.err.println("initialize code!");
        return 0;
    }

    @Override
    public int set_position(IntBuffer index_of_the_particle, DoubleBuffer x, DoubleBuffer y, DoubleBuffer z, int npoints) {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public int commit_parameters() {
        System.err.println("commit particles!");
        return 0;
    }

}
