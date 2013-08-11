from amuse.lab import *
from amuse.support.literature import TrackLiteratureReferences

GRAVITY_CODES = (
    BHTree,
    Hermite,
    PhiGRAPE,
    Octgrav,
    TwoBody,
    Huayno,
    ph4,
    Bonsai,
    Pikachu,
    AarsethZare,
    Adaptb,
    Hacs64,
    HiGPUs,
    Kepler,
    Mercury,
    MI6,
    Mikkola,
    SmallN
)

HYDRO_CODES=(
    Fi,
    Gadget2,
    Athena,
    Capreole,
    MpiAmrVac,
)

RADIATIVE_CODES=(
    SimpleX,
    Mocassin,
    SPHRay
)

STAR_CODES=(
    SSE,
    BSE,
    SeBa,
    EVtwin,
    MESA,
    MakeMeAMassiveStar
)

MISC_CODES=(
    Hop, 
)

stopping_conditions = (
    'collision_detection',
    'pair_detection',
    'escaper_detection',
    'timeout_detection',
    'number_of_steps_detection',
    'out_of_box_detection',
    'density_limit_detection', 
    'internal_energy_limit_detection', 
    'interaction_over_detection'
)

def report(codes, title):
    print title
    print '=' * 80
    print 'number of codes:', len(codes)
    print '=' * 80
    for code in codes:
        try:
            instance = code(redirection="null")
        except Exception as ex:
            print code.__name__, '***' , 'NOT AVAILABLE', '***'
            continue
        if not hasattr(instance, 'stopping_conditions'):
            print code.__name__, '***' , 'NO SUPPORT', '***'
            continue
        print code.__name__, ':',
        instance.initialize_code()
        for name in stopping_conditions:
            sc = getattr(instance.stopping_conditions, name)
            if sc.is_supported():
                print name,
        print 
    print '=' * 80
    
if __name__ == '__main__':
    TrackLiteratureReferences.suppress_output()
    report(GRAVITY_CODES, 'gravity codes')
    report(HYDRO_CODES, 'hydrodynamics codes')
    report(RADIATIVE_CODES, 'radiative transfer codes')
    report(STAR_CODES, 'stellar evolution codes')