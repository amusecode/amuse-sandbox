from amuse.community.aarsethzare.interface import AarsethZare 


if __name__ == '__main__':
    instance = AarsethZare()

    t = [0.,0.,0.]
    m = [1.1, 1.3, 1.0]
    x = [-12.578787679, -58.4511049545, 89.8231028877]
    y = [-198.157193278, -189.734002337, 464.627115644]
    z = [-13.4221707428, 11.3572213977, 0.0]
    vx = [0.0660907760075, -0.0201875848322, -0.0464559933263]
    vy = [0.0723351192845, -0.0927973717159, 0.0410679520178]
    vz = [0.0561101265997, -0.0474777994305, 0.0]
    t, x, y, z, vx, vy, vz, error = instance.call_aarseth_zare(t, m, x, y, z, vx, vy, vz) 
    print 't=', t[0], 'x=', x, y, z, vx, vy, vz, error
    instance.stop()

