def binary_with_neighbours_initial_conditions():
    """
    ad-hoc binary at origin and two stars far away on nearly circular orbit
    orbital period: roughly 200
    """
    from amuse.support.data import core
    from amuse.support.units import nbody_system
    a = 10
    e = 0.9
    G = 1
    m1 = 0.5
    m2 = 0.5
    mu = G * (m1 + m2)
    h = -0.5 * mu / a
    X = 0.5*a*(1-e)
    V = 0.5*math.sqrt(-h)
    #V = 0.5 * math.sqrt( 2 * mu * ( 1 / (a * (1 - e))  - 1 / (2 * a) ) )
    x = []
    x.append([+X,  0,  0])
    x.append([-X,  0,  0])
    v = []
    v.append([ 0, -V,  0])
    v.append([ 0, +V,  0])
    x.append([+20*X,  0,  0])
    x.append([-20*X,  0,  0])
    v.append([ 0, -3*V,  0])
    v.append([ 0, +3*V,  0])
    stars = core.Stars(len(x))
    for i in range(len(stars)):
        star = stars[i]
        star.mass = 0.5 | nbody_system.mass
        star.position = x[i] | nbody_system.length
        star.velocity = v[i] | ( nbody_system.length / nbody_system.time )
        star.radius = 0. | nbody_system.length # obligatory for Hermite() integrator
    return stars
