package nl.esciencecenter.estars;

import nl.esciencecenter.visualization.amuse.planetformation.glue.*;

public class CodeParticle {

    public static final int TYPE_SPHERE = 0;
    public static final int TYPE_STAR = 1;
    public static final int TYPE_SPH_GAS = 2;
    public static final int TYPE_POINT_GAS = 3;
    public static final int TYPE_PLANET = 4;

    private int index;
    private int type;
    private double x;
    private double y;
    private double z;
    private double radius;
    private double red;
    private double green;
    private double blue;
    private double alpha;

    CodeParticle(int index, int type, double x, double y, double z, double radius, double red, double green,
            double blue, double alpha) {
        this.index = index;
        this.type = type;
        this.x = x;
        this.y = y;
        this.z = z;
        this.radius = radius;
        this.red = red;
        this.green = green;
        this.blue = blue;
        this.alpha = alpha;

    }

    public int getIndex() {
        return index;
    }

    public void setIndex(int index) {
        this.index = index;
    }

    public int getType() {
        return type;
    }

    Sphere asSphere() {
        float[] coordinates = new float[] { (float) x, (float) y, (float) z };

        float[] color = new float[] { (float) red, (float) green, (float) blue, (float) alpha };

        return new Sphere(index, coordinates, (float) radius, color);
    }

    Star asStar() {
        float[] coordinates = new float[] { (float) x, (float) y, (float) z };

        float[] color = new float[] { (float) red, (float) green, (float) blue, (float) alpha };
        
        return new Star(index, coordinates, (float) radius, color);
    }

    SPHGas asSphGas() {
        float[] coordinates = new float[] { (float) x, (float) y, (float) z };

        float[] color = new float[] { (float) red, (float) green, (float) blue, (float) alpha };

        return new SPHGas(index, coordinates, color);
    }

    PointGas asPointGas() {
        float[] coordinates = new float[] { (float) x, (float) y, (float) z };

        float[] color = new float[] { (float) red, (float) green, (float) blue, (float) alpha };

        return new PointGas(index, coordinates, color);
    }

    Planet asPlanet() {
        float[] coordinates = new float[] { (float) x, (float) y, (float) z };

        float[] color = new float[] { (float) red, (float) green, (float) blue, (float) alpha };

        return new Planet(index, coordinates, (float) radius, color);
    }

    public void setColor(double red, double green, double blue, double alpha) {
        this.red = red;
        this.green = green;
        this.blue = blue;
        this.alpha = alpha;
    }

    public void setType(int type) {
        this.type = type;
        
    }

    public void setPosition(double x, double y, double z) {
        this.x = x;
        this.y = y;
        this.z = z;
        
    }

    public void setRadius(double radius) {
        this.radius = radius;
    }
}
