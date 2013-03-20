package nl.esciencecenter.estars;

import java.util.ArrayList;

import nl.esciencecenter.visualization.amuse.planetformation.AmuseLib;
import nl.esciencecenter.visualization.amuse.planetformation.glue.*;

import amuse.code.CodeInterface;

public class Code implements CodeInterface {

    private final AmuseLib amuseLib;

    private final ArrayList<CodeParticle> particles;

    private boolean useOctreeForGas = false;
    private boolean useStarShader = true;

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

    private int findEmptyIndex(int start) {
        for (int i = start; i < particles.size(); i++) {
            if (particles.get(i) == null) {
                return i;
            }
        }
        // append the list with a new (empty) particle
        int result = particles.size();
        particles.add(null);
        return result;
    }

    // first parameter is output parameter!
    private int newParticle(int[] index_of_the_particle, double[] x, double[] y, double[] z, double[] radius,
            double[] red, double[] green, double[] blue, double[] opacity, int npoints, int type) {
        System.err.println("new_particle: " + npoints + " new particles being added");
        int newIndex = 0;
        for (int i = 0; i < npoints; i++) {
            newIndex = findEmptyIndex(newIndex);

            // System.err
            // .printf("new particle of type %d at position %.4f,%.4f,%.4f and radius %.4f with rgba %.4f,%.4f,%.4f,%.4f\n",
            // type[i], x[i], y[i], z[i], radius[i], red[i], green[i], blue[i],
            // alpha[i]);

            index_of_the_particle[i] = newIndex;
            CodeParticle particle = new CodeParticle(index_of_the_particle[i], x[i], y[i], z[i], radius[i], red[i],
                    green[i], blue[i], opacity[i], type);

            particles.set(newIndex, particle);
        }

        // System.err.println("new particle result = " +
        // Arrays.toString(index_of_the_particle));

        return 0;
    }

    @Override
    public int new_star_particle(int[] index_of_the_particle, double[] x, double[] y, double[] z, double[] radius,
            double[] red, double[] green, double[] blue, double[] opacity, int npoints) {
        if (useStarShader) {
            return newParticle(index_of_the_particle, x, y, z, radius, red, green, blue, opacity, npoints,
                    CodeParticle.TYPE_STAR);
        } else {
            return newParticle(index_of_the_particle, x, y, z, radius, red, green, blue, opacity, npoints,
                    CodeParticle.TYPE_SPHERE);
        }
    }

    @Override
    public int new_gas_particle(int[] index_of_the_particle, double[] x, double[] y, double[] z, double[] radius,
            double[] red, double[] green, double[] blue, double[] opacity, int npoints) {
        if (useOctreeForGas) {
            return newParticle(index_of_the_particle, x, y, z, radius, red, green, blue, opacity, npoints,
                    CodeParticle.TYPE_STAR);
        } else {
            return newParticle(index_of_the_particle, x, y, z, radius, red, green, blue, opacity, npoints,
                    CodeParticle.TYPE_SPHERE);
        }
    }

    @Override
    public int new_sphere_particle(int[] index_of_the_particle, double[] x, double[] y, double[] z, double[] radius,
            double[] red, double[] green, double[] blue, double[] opacity, int npoints) {
        return newParticle(index_of_the_particle, x, y, z, radius, red, green, blue, opacity, npoints,
                CodeParticle.TYPE_SPHERE);
    }

    @Override
    public int new_marker_particle(int[] index_of_the_particle, double[] x, double[] y, double[] z, double[] radius,
            double[] red, double[] green, double[] blue, double[] opacity, int npoints) {
        return newParticle(index_of_the_particle, x, y, z, radius, red, green, blue, opacity, npoints,
                CodeParticle.TYPE_MARKER);
    }

    @Override
    public int commit_particles() {
        System.err.println("commit particles! (does nothing)");
        return 0;
    }

    @Override
    public int recommit_particles() {
        System.err.println("recommit particles! (does nothing)");
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
            // System.err.printf("particle %d position updated to %.4f,%.4f,%.4f\n",
            // index_of_the_particle[i], x[i], y[i], z[i]);
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
        ArrayList<Sphere> spheres = new ArrayList<Sphere>();
        ArrayList<Star> stars = new ArrayList<Star>();
        ArrayList<SPHGas> sphGas = new ArrayList<SPHGas>();
        ArrayList<PointGas> pointGas = new ArrayList<PointGas>();

        for (CodeParticle particle : particles) {
            if (particle != null) {
                switch (particle.getType()) {
                case CodeParticle.TYPE_STAR:
                    // System.err.println("star added: " + particle);
                    stars.add(particle.asStar());
                    break;
                case CodeParticle.TYPE_POINT_GAS:
                    pointGas.add(particle.asPointGas());
                    // System.err.println("point gass added: " + particle);
                    break;
                case CodeParticle.TYPE_OCTREE_GAS:
                    sphGas.add(particle.asSphGas());
                    // System.err.println("sph gass added: " + particle);
                    break;
                case CodeParticle.TYPE_SPHERE:
                    // System.err.println("sphere added: " + particle);
                    spheres.add(particle.asSphere());
                    break;
                case CodeParticle.TYPE_MARKER:
                    // System.err.println("sphere added: " + particle);
                    spheres.add(particle.asSphere());
                    break;

                default:
                    // FIXME:we actually want to throw an error
                    System.err.println("Unknown particle type " + particle.getType());
                    return 1;
                }
            }

        }

        System.err.printf(
                "storing view with %d spheres, %d stars, %d sph_gas and %d pointgas particles. Description = \"%s\"\n",
                spheres.size(), stars.size(), sphGas.size(), pointGas.size(), description);

        amuseLib.addScene(new Scene(description, spheres.toArray(new Sphere[spheres.size()]), stars
                .toArray(new Star[stars.size()]), sphGas.toArray(new SPHGas[sphGas.size()]), pointGas
                .toArray(new PointGas[pointGas.size()])));

        // succes!
        return 0;
    }

    @Override
    public int cleanup_code() {
        System.err.println("cleanup code! (does nothing)");
        return 0;
    }

    @Override
    public int delete_particle(int[] index_of_the_particle, int npoints) {
        for (int i = 0; i < npoints; i++) {
            int index = index_of_the_particle[i];
            if (index >= particles.size() || particles.get(index) == null) {
                System.err.println("No particle with index " + index + " exists");
                return 1;
            }

            particles.set(index, null);
        }
        return 0;
    }

    @Override
    // parameter "red" is an output parameter!
    // parameter "green" is an output parameter!
    // parameter "blue" is an output parameter!
    public int get_color(int[] index_of_the_particle, double[] red, double[] green, double[] blue, int npoints) {
        for (int i = 0; i < npoints; i++) {
            CodeParticle particle = findParticle(index_of_the_particle[i]);

            if (particle == null) {
                return 1;
            }

            red[index_of_the_particle[i]] = particle.getRed();
            green[index_of_the_particle[i]] = particle.getGreen();
            green[index_of_the_particle[i]] = particle.getBlue();

        }
        return 0;
    }

    @Override
    // parameter "alpha" is an output parameter!
    public int get_opacity(int[] index_of_the_particle, double[] opacity, int npoints) {
        for (int i = 0; i < npoints; i++) {
            CodeParticle particle = findParticle(index_of_the_particle[i]);

            if (particle == null) {
                return 1;
            }

            opacity[index_of_the_particle[i]] = particle.getOpacity();
        }
        return 0;
    }

    @Override
    // parameter "x" is an output parameter!
    // parameter "y" is an output parameter!
    // parameter "z" is an output parameter!
    public int get_position(int[] index_of_the_particle, double[] x, double[] y, double[] z, int npoints) {
        for (int i = 0; i < npoints; i++) {
            CodeParticle particle = findParticle(index_of_the_particle[i]);

            if (particle == null) {
                return 1;
            }

            x[index_of_the_particle[i]] = particle.getX();
            y[index_of_the_particle[i]] = particle.getY();
            z[index_of_the_particle[i]] = particle.getZ();

        }
        return 0;
    }

    @Override
    public int recommit_parameters() {
        System.err.println("recommit_parameters! (does nothing)");
        return 0;
    }

    @Override
    public int set_opacity(int[] index_of_the_particle, double[] opacity, int npoints) {
        for (int i = 0; i < npoints; i++) {
            CodeParticle particle = findParticle(index_of_the_particle[i]);

            if (particle == null) {
                return 1;
            }

            particle.setOpacity(opacity[i]);
        }
        return 0;
    }

    @Override
    public int set_use_star_shader_flag(int use_star_shader_flag) {
        if (use_star_shader_flag == 1) {
            useStarShader = true;
        } else {
            useStarShader = false;
        }
        return 0;
    }

    @Override
    public int get_use_star_shader_flag(int[] use_star_shader_flag) {
        if (useStarShader) {
            use_star_shader_flag[0] = 1;
        } else {
            use_star_shader_flag[0] = 0;
        }
        return 0;
    }

    @Override
    public int set_use_octree_for_gas_flag(int use_octree_for_gas_flag) {
        if (use_octree_for_gas_flag == 1) {
            useOctreeForGas = true;
        } else {
            useOctreeForGas = false;
        }
        return 0;
    }

    @Override
    public int get_use_octree_for_gas_flag(int[] use_octree_for_gas_flag) {
        if (useOctreeForGas) {
            use_octree_for_gas_flag[0] = 1;
        } else {
            use_octree_for_gas_flag[0] = 0;
        }
        return 0;
    }

}
