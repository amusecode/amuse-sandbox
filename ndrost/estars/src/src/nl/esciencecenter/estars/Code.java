package nl.esciencecenter.estars;

import java.util.ArrayList;

import nl.esciencecenter.visualization.amuse.planetformation.AmuseLib;
import nl.esciencecenter.visualization.amuse.planetformation.glue.*;

import amuse.code.CodeInterface;

public class Code implements CodeInterface {

    private final AmuseLib amuseLib;

    private final ArrayList<CodeParticle> particles;

    public Code() {
        amuseLib = new AmuseLib();

        particles = new ArrayList<CodeParticle>();

    }

    public int internal__get_message_polling_interval(int[] interval) {
        System.err.println("get_message_polling_interval! (does nothing)");
        return 0;
    }

    public int internal__set_message_polling_interval(int interval) {
        System.err.println("get_message_polling_interval! (does nothing)");
        return 0;
    }

    private CodeParticle findParticle(int index) {
        if (index >= particles.size()) {
            System.err.println("No particle with index " + index + " exists");
            return null;
        }
        CodeParticle particle = particles.get(index);

        if (particle == null) {
            System.err.println("No particle with index " + index + " exists");
        }

        // particle could still be null here
        return particle;

    }

    @Override
    public int initialize_code() {
        System.err.println("initialize code! (does nothing)");
        return 0;
    }

    @Override
    public int commit_parameters() {
        System.err.println("commit parameters! (does nothing)");
        return 0;
    }

    @Override
    // first parameter is output parameter!
    public int new_particle(int[] index_of_the_particle, int[] type, double[] x, double[] y, double[] z,
            double[] radius, double[] red, double[] green, double[] blue, double[] alpha, int npoints) {
        System.err.println("new_particle");
        for (int i = 0; i < npoints; i++) {
            System.err
                    .printf("new particle of type %d at position %.4f,%.4f,%.4f and radius %.4f with rgba %.4f,%.4f,%.4f,%.4f\n",
                            type[i], x[i], y[i], z[i], radius[i], red[i], green[i], blue[i], alpha[i]);

            index_of_the_particle[i] = particles.size();
            CodeParticle particle = new CodeParticle(index_of_the_particle[i], type[i], x[i], y[i], z[i], radius[i],
                    red[i], green[i], blue[i], alpha[i]);

            particles.add(particle);
        }

        return 0;
    }

    @Override
    public int commit_particles() {
        System.err.println("commit particles! (does nothing)");
        return 0;
    }

    @Override
    public int set_type(int[] index_of_the_particle, int[] type, int npoints) {
        for (int i = 0; i < npoints; i++) {
            CodeParticle particle = findParticle(index_of_the_particle[i]);

            if (particle == null) {
                return 1;
            }

            particle.setType(type[i]);
        }
        return 0;
    }

    @Override
    public int set_position(int[] index_of_the_particle, double[] x, double[] y, double[] z, int npoints) {
        for (int i = 0; i < npoints; i++) {
            CodeParticle particle = findParticle(index_of_the_particle[i]);

            if (particle == null) {
                return 1;
            }

            particle.setPosition(x[i], y[i], z[i]);
        }
        return 0;
    }

    @Override
    public int set_radius(int[] index_of_the_particle, double[] radius, int npoints) {
        for (int i = 0; i < npoints; i++) {
            CodeParticle particle = findParticle(index_of_the_particle[i]);

            if (particle == null) {
                return 1;
            }

            particle.setRadius(radius[i]);
        }
        return 0;
    }

    @Override
    // radius is output parameter!
    public int get_radius(int[] index_of_the_particle, double[] radius, int npoints) {

        for (int i = 0; i < npoints; i++) {
            CodeParticle particle = findParticle(index_of_the_particle[i]);

            if (particle == null) {
                return 1;
            }

            radius[i] = particle.getRadius();
        }
        return 0;
    }

    @Override
    public int set_color(int[] index_of_the_particle, double[] red, double[] green, double[] blue, int npoints) {
        for (int i = 0; i < npoints; i++) {
            CodeParticle particle = findParticle(index_of_the_particle[i]);

            if (particle == null) {
                return 1;
            }

            particle.setColor(red[i], green[i], blue[i]);
        }
        return 0;
    }

    @Override
    public int store_view(String description) {
        System.err.println("store view with description " + description);

        // Sphere[] spheres, Star[] stars, Planet[] planets, SPHGas[] sphGas,
        // PointGas[] pointGas

        ArrayList<Sphere> spheres = new ArrayList<Sphere>();
        ArrayList<Star> stars = new ArrayList<Star>();
        ArrayList<Planet> planets = new ArrayList<Planet>();
        ArrayList<SPHGas> sphGas = new ArrayList<SPHGas>();
        ArrayList<PointGas> pointGas = new ArrayList<PointGas>();

        for (CodeParticle particle : particles) {
            if (particle != null) {
                switch (particle.getType()) {
                case CodeParticle.TYPE_SPHERE:
                    System.err.println("sphere added: " + particle);
                    spheres.add(particle.asSphere());
                    break;
                case CodeParticle.TYPE_STAR:
                    System.err.println("star added: " + particle);
                    stars.add(particle.asStar());
                    break;
                case CodeParticle.TYPE_PLANET:
                    System.err.println("planet added: " + particle);
                    planets.add(particle.asPlanet());
                    break;
                case CodeParticle.TYPE_SPH_GAS:
                    sphGas.add(particle.asSphGas());
                    System.err.println("sph gass added: " + particle);
                    break;
                case CodeParticle.TYPE_POINT_GAS:
                    pointGas.add(particle.asPointGas());
                    System.err.println("point gass added: " + particle);
                    break;
                default:
                    // FIXME:we actually want to throw an error
                    System.err.println("Unknown particle type " + particle.getType());
                    return 1;
                }
            }

        }
        amuseLib.addScene(new Scene(description, spheres.toArray(new Sphere[spheres.size()]),
                stars.toArray(new Star[stars.size()]), planets.toArray(new Planet[planets.size()]), sphGas
                        .toArray(new SPHGas[sphGas.size()]), pointGas.toArray(new PointGas[pointGas.size()])));

        // succes!
        return 0;
    }

    @Override
    public int cleanup_code() {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public int delete_particle(int[] index_of_the_particle, int npoints) {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    // parameter "red" is an output parameter!
    // parameter "green" is an output parameter!
    // parameter "blue" is an output parameter!
    public int get_color(int[] index_of_the_particle, double[] red, double[] green, double[] blue, int npoints) {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    // parameter "alpha" is an output parameter!
    public int get_opacity(int[] index_of_the_particle, double[] alpha, int npoints) {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    // parameter "x" is an output parameter!
    // parameter "y" is an output parameter!
    // parameter "z" is an output parameter!
    public int get_position(int[] index_of_the_particle, double[] x, double[] y, double[] z, int npoints) {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public int get_type(int[] index_of_the_particle, double[] type, int npoints) {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public int recommit_parameters() {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public int set_opacity(int[] index_of_the_particle, double[] alpha, int npoints) {
        // TODO Auto-generated method stub
        return 0;
    }
}
