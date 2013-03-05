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
    
    private CodeParticle findParticle(int index) {
        if (index >= particles.size()) {
            System.err.println("No particle with index " + index + " exists");
            return null;
        }
        CodeParticle particle = particles.get(index);

        if (particle == null) {
            System.err.println("No particle with index " + index + " exists");
        }
        
        //particle could still be null here
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
            System.err.printf(
                    "new particle of type %d at position %.4f,%.4f,%.4f and radius %.4f with rgba %.4f,%.4f,%.4f,%.4f\n",
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
    public int set_color(int[] index_of_the_particle, double[] red, double[] green, double[] blue, double[] alpha,
            int npoints) {
        for (int i = 0; i < npoints; i++) {
            CodeParticle particle = findParticle(index_of_the_particle[i]);
            
            if (particle == null) {
                return 1;
            }
            
            particle.setColor(red[i], green[i], blue[i], alpha[i]);
        }
        return 0;
    }

    @Override
    public int store_view(double time) {
        System.err.println("store view @ " + time);

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
                    spheres.add(particle.asSphere());
                    break;
                case CodeParticle.TYPE_STAR:
                    stars.add(particle.asStar());
                    break;
                case CodeParticle.TYPE_PLANET:
                    planets.add(particle.asPlanet());
                    break;
                case CodeParticle.TYPE_SPH_GAS:
                    sphGas.add(particle.asSphGas());
                    break;
                case CodeParticle.TYPE_POINT_GAS:
                    pointGas.add(particle.asPointGas());
                    break;
                default:
                    // FIXME:we actually want to throw an error
                    System.err.println("Unknown particle type " + particle.getType());
                    return 1;
                }
            }

        }
        amuseLib.addScene(new Scene(time, spheres.toArray(new Sphere[spheres.size()]), stars.toArray(new Star[stars
                .size()]), planets.toArray(new Planet[planets.size()]), sphGas.toArray(new SPHGas[sphGas.size()]),
                pointGas.toArray(new PointGas[pointGas.size()])));

        //succes!
        return 0;
    }

    @Override
    public int cleanup_code() {
        System.err.println("cleanup_code! (does nothing)");
        return 0;
    }

    @Override
    public int recommit_parameters() {
        System.err.println("recommit_particles! (does nothing)");
        return 0;
    }

    
   

    

    

}
