
#ifndef NOMPI
    #include <mpi.h>
#endif
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#ifdef WIN32
	#include <winsock2.h>
#else
	#include <sys/socket.h>
	#include <netinet/in.h>
	#include <netdb.h>
	#include <unistd.h>
	#include <netinet/tcp.h>
#endif

#include "interface.h"
#include "worker_code.h"
#include "stopcond.h"

static bool NEEDS_MPI = true;

static int MAX_INTS_IN = 2;
static int MAX_INTS_OUT = 3;
static int MAX_STRINGS_IN = 1;
static int MAX_STRINGS_OUT = 1;
static int MAX_DOUBLES_IN = 8;
static int MAX_DOUBLES_OUT = 8;
static int MAX_LONGS_IN = 0;
static int MAX_LONGS_OUT = 0;
static int MAX_BOOLEANS_IN = 1;
static int MAX_BOOLEANS_OUT = 1;
static int MAX_FLOATS_IN = 0;
static int MAX_FLOATS_OUT = 0;


static int ERROR_FLAG = 256;
static int HEADER_SIZE = 10; //integers

static int HEADER_FLAGS = 0;
static int HEADER_CALL_ID = 1;
static int HEADER_FUNCTION_ID = 2;
static int HEADER_CALL_COUNT = 3;
static int HEADER_INTEGER_COUNT = 4;
static int HEADER_LONG_COUNT = 5;
static int HEADER_FLOAT_COUNT = 6;
static int HEADER_DOUBLE_COUNT = 7;
static int HEADER_BOOLEAN_COUNT = 8;
static int HEADER_STRING_COUNT = 9;

static bool TRUE_BYTE = 1;
static bool FALSE_BYTE = 0;

static bool mpiIntercom = false;

static int socketfd = 0;

static int * header_in;
static int * header_out;

static int * ints_in;
static int * ints_out;

static long long int * longs_in;
static long long int * longs_out;

static float * floats_in;
static float * floats_out;

static double * doubles_in;
static double * doubles_out;

static int * booleans_in;
static int * booleans_out;

/* sizes of strings */
static int * string_sizes_in;
static int * string_sizes_out;

/* pointers to input and output strings (contents not stored here) */
static char * * strings_in;
static char * * strings_out;

/* actual string data */
static char * characters_in = 0;
static char * characters_out = 0;


static int polling_interval = 0;
#ifndef NOMPI
#define MAX_COMMUNICATORS 2048
static char portname_buffer[MPI_MAX_PORT_NAME+1];
static MPI_Comm communicators[MAX_COMMUNICATORS];
static int lastid = -1;
static int activeid = -1;
static int id_to_activate = -1;
#else
static const char * empty_string = "";
#endif

int internal__get_message_polling_interval(int * outval)
{
    *outval = polling_interval;
    
    return 0;
}

int internal__set_message_polling_interval(int inval)
{
    polling_interval = inval;
    
    return 0;
}
int internal__open_port(char ** output)
{
#ifndef NOMPI
    MPI_Open_port(MPI_INFO_NULL, portname_buffer);
    *output = portname_buffer;
#else
    *output = empty_string;
#endif
    return 0;
}
int internal__accept_on_port(char * port_identifier, int * comm_identifier)
{
#ifndef NOMPI
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    lastid++;
    if(lastid >= MAX_COMMUNICATORS) {
        lastid--;
        return -1;
    }
    if(rank == 0){
        MPI_Comm merged;
        MPI_Comm communicator;
        MPI_Comm_accept(port_identifier, MPI_INFO_NULL, 0,  MPI_COMM_SELF, &communicator);
        MPI_Intercomm_merge(communicator, 0, &merged);
        MPI_Intercomm_create(MPI_COMM_WORLD,0,merged, 1, 65, &communicators[lastid]);
        MPI_Comm_disconnect(&merged);
        MPI_Comm_disconnect(&communicator);
    } else {
        MPI_Intercomm_create(MPI_COMM_WORLD,0, MPI_COMM_NULL, 1, 65, &communicators[lastid]);
    }
    *comm_identifier = lastid;
#else
    *comm_identifier = -1;
#endif
    return 0;
}


int internal__connect_to_port(char * port_identifier, int * comm_identifier)
{
#ifndef NOMPI
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    lastid++;
    if(lastid >= MAX_COMMUNICATORS) {
        lastid--;
        return -1;
    }
    if(rank == 0){
        MPI_Comm merged;
        MPI_Comm communicator;
        MPI_Comm_connect(port_identifier, MPI_INFO_NULL, 0,  MPI_COMM_SELF, &communicator);
        MPI_Intercomm_merge(communicator, 1, &merged);
        MPI_Intercomm_create(MPI_COMM_WORLD, 0, merged, 0, 65, &communicators[lastid]);
        MPI_Comm_disconnect(&merged);
        MPI_Comm_disconnect(&communicator);
    } else {
        MPI_Intercomm_create(MPI_COMM_WORLD, 0, MPI_COMM_NULL, 1, 65, &communicators[lastid]);
    }
    *comm_identifier = lastid;
#else
    *comm_identifier = -1;
#endif
    return 0;
}

int internal__activate_communicator(int comm_identifier){
#ifndef NOMPI
    if(comm_identifier < 0 || comm_identifier > lastid) {
        return -1;
    }
    id_to_activate = comm_identifier;
#endif
    return 0;
}



#ifndef NOMPI
#include <unistd.h>

int mpi_recv_header(MPI_Comm & parent)
{
    MPI_Request header_request;
    MPI_Status request_status;
   
    MPI_Irecv(header_in, HEADER_SIZE, MPI_INT, 0, 989, parent, &header_request);
        
    if(polling_interval > 0)
    {
        int is_finished = 0;
        MPI_Test(&header_request, &is_finished, &request_status);
        while(!is_finished) {
            usleep(polling_interval);
            MPI_Test(&header_request, &is_finished, &request_status);
        }
        MPI_Wait(&header_request, &request_status);
    } else {
        MPI_Wait(&header_request, &request_status);
    }
    return 0;
}
#endif


bool handle_call() {
  int call_count = header_in[HEADER_CALL_COUNT];
  
  switch(header_in[HEADER_FUNCTION_ID]) {
    case 0:
      return false;
      break;
    case 20284990:
      ints_out[0] = get_mass(
        &ints_in[0] ,
        &doubles_out[0] ,
        call_count
      );
      for (int i = 1 ; i < call_count; i++){
        ints_out[i] = ints_out[0];
      }
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 20920053:
      ints_out[0] = commit_particles();
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 26035633:
      ints_out[0] = set_unit_length(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 27307077:
      ints_out[0] = get_tree_domain_update_frequency(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 37921492:
      ints_out[0] = set_stopping_condition_timeout_parameter(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 44188957:
      ints_out[0] = get_time(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 74884513:
      ints_out[0] = get_epsilon_dm_part(
        &ints_in[0] ,
        &doubles_out[0] ,
        call_count
      );
      for (int i = 1 ; i < call_count; i++){
        ints_out[i] = ints_out[0];
      }
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 95026274:
      ints_out[0] = set_nsmooth(
        ints_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 104547857:
      ints_out[0] = set_mass(
        &ints_in[0] ,
        &doubles_in[0] ,
        call_count
      );
      for (int i = 1 ; i < call_count; i++){
        ints_out[i] = ints_out[0];
      }
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 107212404:
      ints_out[0] = get_stopping_condition_maximum_density_parameter(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 115569384:
      ints_out[0] = get_info_file(
        &strings_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_STRING_COUNT] = 1 * call_count;
      break;
    
    case 123510229:
      ints_out[0] = get_time_limit_cpu(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 128926247:
      ints_out[0] = get_index_of_first_particle(
        &ints_out[1]
      );
      header_out[HEADER_INTEGER_COUNT] = 2 * call_count;
      break;
    
    case 159095171:
      for (int i = 0 ; i < call_count; i++){
        ints_out[i] = enable_stopping_condition(
          ints_in[i]
        );
      }
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 205426934:
      ints_out[0] = get_total_radius(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 205548127:
      ints_out[0] = get_stopping_condition_maximum_internal_energy_parameter(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 210995141:
      for (int i = 0 ; i < call_count; i++){
        ints_out[i] = get_potential_at_point(
          doubles_in[i] ,
          doubles_in[( 1 * call_count) + i] ,
          doubles_in[( 2 * call_count) + i] ,
          doubles_in[( 3 * call_count) + i] ,
          &doubles_out[i]
        );
      }
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 211853413:
      ints_out[0] = get_periodic_boundaries_flag(
        &booleans_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_BOOLEAN_COUNT] = 1 * call_count;
      break;
    
    case 212446545:
      ints_out[0] = get_softening_gas_max_phys(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 219539042:
      for (int i = 0 ; i < call_count; i++){
        ints_out[i] = get_number_of_stopping_conditions_set(
          &ints_out[( 1 * call_count) + i]
        );
      }
      header_out[HEADER_INTEGER_COUNT] = 2 * call_count;
      break;
    
    case 242056790:
      ints_out[0] = set_min_gas_hsmooth_fractional(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 271440754:
      for (int i = 0 ; i < call_count; i++){
        ints_out[i] = is_stopping_condition_set(
          ints_in[i] ,
          &ints_out[( 1 * call_count) + i]
        );
      }
      header_out[HEADER_INTEGER_COUNT] = 2 * call_count;
      break;
    
    case 287868353:
      ints_out[0] = set_periodic_boundaries_flag(
        booleans_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 300486370:
      ints_out[0] = get_gdgtol(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 335500686:
      ints_out[0] = get_gadget_output_directory(
        &strings_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_STRING_COUNT] = 1 * call_count;
      break;
    
    case 350645236:
      ints_out[0] = get_n_neighbours(
        &ints_in[0] ,
        &doubles_out[0] ,
        call_count
      );
      for (int i = 1 ; i < call_count; i++){
        ints_out[i] = ints_out[0];
      }
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 358946551:
      ints_out[0] = get_interpret_heat_as_feedback_flag(
        &booleans_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_BOOLEAN_COUNT] = 1 * call_count;
      break;
    
    case 383112453:
      ints_out[0] = internal__connect_to_port(
        strings_in[0] ,
        &ints_out[1]
      );
      header_out[HEADER_INTEGER_COUNT] = 2 * call_count;
      break;
    
    case 384567015:
      ints_out[0] = get_total_mass(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 416721977:
      ints_out[0] = evolve_model(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 440694980:
      ints_out[0] = set_state_sph(
        &ints_in[0] ,
        &doubles_in[0] ,
        &doubles_in[( 1 * call_count)] ,
        &doubles_in[( 2 * call_count)] ,
        &doubles_in[( 3 * call_count)] ,
        &doubles_in[( 4 * call_count)] ,
        &doubles_in[( 5 * call_count)] ,
        &doubles_in[( 6 * call_count)] ,
        &doubles_in[( 7 * call_count)] ,
        call_count
      );
      for (int i = 1 ; i < call_count; i++){
        ints_out[i] = ints_out[0];
      }
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 461785899:
      ints_out[0] = get_nsmooth(
        &ints_out[1]
      );
      header_out[HEADER_INTEGER_COUNT] = 2 * call_count;
      break;
    
    case 466707249:
      ints_out[0] = get_epsgas(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 508605261:
      ints_out[0] = set_stopping_condition_out_of_box_parameter(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 508859371:
      ints_out[0] = set_type_of_timestep_criterion(
        ints_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 513205695:
      ints_out[0] = get_redshift_begin(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 515213424:
      ints_out[0] = get_omega_lambda(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 525520780:
      ints_out[0] = get_unit_velocity(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 542058817:
      ints_out[0] = set_eps2(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 576132811:
      ints_out[0] = set_time_limit_cpu(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 632979349:
      ints_out[0] = set_stopping_condition_number_of_steps_parameter(
        ints_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 639951605:
      ints_out[0] = get_stopping_condition_timeout_parameter(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 653682513:
      ints_out[0] = get_begin_time(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 658631024:
      ints_out[0] = get_eps2(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 662348285:
      ints_out[0] = get_stopping_condition_minimum_internal_energy_parameter(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 666749878:
      ints_out[0] = get_epsilon(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 678380482:
      ints_out[0] = get_index_of_next_particle(
        ints_in[0] ,
        &ints_out[1]
      );
      header_out[HEADER_INTEGER_COUNT] = 2 * call_count;
      break;
    
    case 686376645:
      ints_out[0] = get_alpha_visc(
        &ints_in[0] ,
        &doubles_out[0] ,
        call_count
      );
      for (int i = 1 ; i < call_count; i++){
        ints_out[i] = ints_out[0];
      }
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 726484084:
      for (int i = 0 ; i < call_count; i++){
        ints_out[i] = new_sph_particle(
          &ints_out[( 1 * call_count) + i] ,
          doubles_in[i] ,
          doubles_in[( 1 * call_count) + i] ,
          doubles_in[( 2 * call_count) + i] ,
          doubles_in[( 3 * call_count) + i] ,
          doubles_in[( 4 * call_count) + i] ,
          doubles_in[( 5 * call_count) + i] ,
          doubles_in[( 6 * call_count) + i] ,
          doubles_in[( 7 * call_count) + i]
        );
      }
      header_out[HEADER_INTEGER_COUNT] = 2 * call_count;
      break;
    
    case 727361823:
      ints_out[0] = internal__set_message_polling_interval(
        ints_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 728786188:
      for (int i = 0 ; i < call_count; i++){
        ints_out[i] = delete_particle(
          ints_in[i]
        );
      }
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 733749514:
      for (int i = 0 ; i < call_count; i++){
        ints_out[i] = is_stopping_condition_enabled(
          ints_in[i] ,
          &ints_out[( 1 * call_count) + i]
        );
      }
      header_out[HEADER_INTEGER_COUNT] = 2 * call_count;
      break;
    
    case 748943991:
      ints_out[0] = set_omega_lambda(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 751474078:
      ints_out[0] = set_gadget_output_directory(
        strings_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 763726023:
      ints_out[0] = get_isotherm(
        &ints_out[1]
      );
      header_out[HEADER_INTEGER_COUNT] = 2 * call_count;
      break;
    
    case 771423968:
      ints_out[0] = set_interpret_heat_as_feedback_flag(
        booleans_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 780854995:
      ints_out[0] = set_max_size_timestep(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 783846686:
      ints_out[0] = get_dtalpha_visc(
        &ints_in[0] ,
        &doubles_out[0] ,
        call_count
      );
      for (int i = 1 ; i < call_count; i++){
        ints_out[i] = ints_out[0];
      }
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 809889639:
      ints_out[0] = get_unit_time(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 810771933:
      ints_out[0] = get_box_size(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 812400739:
      ints_out[0] = get_unit_length(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 818436780:
      ints_out[0] = set_omega_zero(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 835969050:
      ints_out[0] = get_potential(
        &ints_in[0] ,
        &doubles_out[0] ,
        call_count
      );
      for (int i = 1 ; i < call_count; i++){
        ints_out[i] = ints_out[0];
      }
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 847059843:
      ints_out[0] = set_gdgop(
        ints_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 869411910:
      ints_out[0] = get_unit_mass(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 887125873:
      ints_out[0] = synchronize_model();
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 898908835:
      ints_out[0] = set_epsgas(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 919768251:
      ints_out[0] = internal__get_message_polling_interval(
        &ints_out[1]
      );
      header_out[HEADER_INTEGER_COUNT] = 2 * call_count;
      break;
    
    case 929181341:
      ints_out[0] = set_state(
        &ints_in[0] ,
        &doubles_in[0] ,
        &doubles_in[( 1 * call_count)] ,
        &doubles_in[( 2 * call_count)] ,
        &doubles_in[( 3 * call_count)] ,
        &doubles_in[( 4 * call_count)] ,
        &doubles_in[( 5 * call_count)] ,
        &doubles_in[( 6 * call_count)] ,
        call_count
      );
      for (int i = 1 ; i < call_count; i++){
        ints_out[i] = ints_out[0];
      }
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 953558780:
      ints_out[0] = get_stopping_condition_minimum_density_parameter(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 965901567:
      ints_out[0] = set_epsilon(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 967950880:
      ints_out[0] = get_state(
        &ints_in[0] ,
        &doubles_out[0] ,
        &doubles_out[( 1 * call_count)] ,
        &doubles_out[( 2 * call_count)] ,
        &doubles_out[( 3 * call_count)] ,
        &doubles_out[( 4 * call_count)] ,
        &doubles_out[( 5 * call_count)] ,
        &doubles_out[( 6 * call_count)] ,
        call_count
      );
      for (int i = 1 ; i < call_count; i++){
        ints_out[i] = ints_out[0];
      }
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 7 * call_count;
      break;
    
    case 989934113:
      ints_out[0] = set_min_size_timestep(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 1000173217:
      ints_out[0] = set_redshift_begin(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 1018099458:
      ints_out[0] = get_gdgop(
        &ints_out[1]
      );
      header_out[HEADER_INTEGER_COUNT] = 2 * call_count;
      break;
    
    case 1024680297:
      ints_out[0] = get_time_step(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 1050085724:
      ints_out[0] = recommit_particles();
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 1062056798:
      ints_out[0] = set_box_size(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 1068248944:
      ints_out[0] = set_gdgtol(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 1071152125:
      ints_out[0] = get_kinetic_energy(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 1082457792:
      ints_out[0] = get_number_of_particles(
        &ints_out[1]
      );
      header_out[HEADER_INTEGER_COUNT] = 2 * call_count;
      break;
    
    case 1098838617:
      ints_out[0] = get_stopping_condition_number_of_steps_parameter(
        &ints_out[1]
      );
      header_out[HEADER_INTEGER_COUNT] = 2 * call_count;
      break;
    
    case 1121294785:
      ints_out[0] = get_redshift_max(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 1133889020:
      ints_out[0] = set_internal_energy(
        &ints_in[0] ,
        &doubles_in[0] ,
        call_count
      );
      for (int i = 1 ; i < call_count; i++){
        ints_out[i] = ints_out[0];
      }
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 1137702459:
      for (int i = 0 ; i < call_count; i++){
        ints_out[i] = disable_stopping_condition(
          ints_in[i]
        );
      }
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 1139107325:
      ints_out[0] = get_hubble_param(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 1142015244:
      ints_out[0] = set_alpha(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 1157309227:
      ints_out[0] = set_unit_time(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 1175408297:
      ints_out[0] = get_time_between_statistics(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 1178602732:
      ints_out[0] = get_internal_energy(
        &ints_in[0] ,
        &doubles_out[0] ,
        call_count
      );
      for (int i = 1 ; i < call_count; i++){
        ints_out[i] = ints_out[0];
      }
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 1193548493:
      ints_out[0] = get_interpret_kicks_as_feedback_flag(
        &booleans_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_BOOLEAN_COUNT] = 1 * call_count;
      break;
    
    case 1198317578:
      ints_out[0] = set_unit_mass(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 1221190952:
      for (int i = 0 ; i < call_count; i++){
        ints_out[i] = set_acceleration(
          ints_in[i] ,
          doubles_in[i] ,
          doubles_in[( 1 * call_count) + i] ,
          doubles_in[( 2 * call_count) + i]
        );
      }
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 1231060790:
      for (int i = 0 ; i < call_count; i++){
        ints_out[i] = get_center_of_mass_position(
          &doubles_out[i] ,
          &doubles_out[( 1 * call_count) + i] ,
          &doubles_out[( 2 * call_count) + i]
        );
      }
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 3 * call_count;
      break;
    
    case 1236118821:
      ints_out[0] = set_time_step(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 1253068953:
      ints_out[0] = set_err_tol_int_accuracy(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 1254238607:
      ints_out[0] = get_alpha(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 1261964575:
      ints_out[0] = get_comoving_integration_flag(
        &booleans_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_BOOLEAN_COUNT] = 1 * call_count;
      break;
    
    case 1273715782:
      for (int i = 0 ; i < call_count; i++){
        ints_out[i] = get_hydro_state_at_point(
          doubles_in[i] ,
          doubles_in[( 1 * call_count) + i] ,
          doubles_in[( 2 * call_count) + i] ,
          doubles_in[( 3 * call_count) + i] ,
          doubles_in[( 4 * call_count) + i] ,
          doubles_in[( 5 * call_count) + i] ,
          &doubles_out[i] ,
          &doubles_out[( 1 * call_count) + i] ,
          &doubles_out[( 2 * call_count) + i] ,
          &doubles_out[( 3 * call_count) + i] ,
          &doubles_out[( 4 * call_count) + i]
        );
      }
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 5 * call_count;
      break;
    
    case 1278969821:
      ints_out[0] = set_courant(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 1305591623:
      ints_out[0] = get_timings_file(
        &strings_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_STRING_COUNT] = 1 * call_count;
      break;
    
    case 1315680918:
      for (int i = 0 ; i < call_count; i++){
        ints_out[i] = get_center_of_mass_velocity(
          &doubles_out[i] ,
          &doubles_out[( 1 * call_count) + i] ,
          &doubles_out[( 2 * call_count) + i]
        );
      }
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 3 * call_count;
      break;
    
    case 1317242279:
      for (int i = 0 ; i < call_count; i++){
        ints_out[i] = get_radius(
          ints_in[i] ,
          &doubles_out[i]
        );
      }
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 1320443362:
      ints_out[0] = get_bh_tol(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 1339484635:
      ints_out[0] = get_omega_zero(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 1341990430:
      ints_out[0] = get_omega_baryon(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 1371060478:
      ints_out[0] = get_energy_file(
        &strings_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_STRING_COUNT] = 1 * call_count;
      break;
    
    case 1377054868:
      ints_out[0] = get_courant(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 1393481104:
      ints_out[0] = set_stopping_condition_minimum_internal_energy_parameter(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 1431402304:
      ints_out[0] = get_cpu_file(
        &strings_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_STRING_COUNT] = 1 * call_count;
      break;
    
    case 1432951376:
      ints_out[0] = get_smoothing_length(
        &ints_in[0] ,
        &doubles_out[0] ,
        call_count
      );
      for (int i = 1 ; i < call_count; i++){
        ints_out[i] = ints_out[0];
      }
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 1486820981:
      ints_out[0] = get_d_internal_energy_dt(
        &ints_in[0] ,
        &doubles_out[0] ,
        call_count
      );
      for (int i = 1 ; i < call_count; i++){
        ints_out[i] = ints_out[0];
      }
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 1496005460:
      ints_out[0] = get_density(
        &ints_in[0] ,
        &doubles_out[0] ,
        call_count
      );
      for (int i = 1 ; i < call_count; i++){
        ints_out[i] = ints_out[0];
      }
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 1508430886:
      ints_out[0] = set_begin_time(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 1512123058:
      ints_out[0] = set_nsmtol(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 1515168189:
      ints_out[0] = set_cpu_file(
        strings_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 1542847027:
      ints_out[0] = get_epsilon_gas_part(
        &ints_in[0] ,
        &doubles_out[0] ,
        call_count
      );
      for (int i = 1 ; i < call_count; i++){
        ints_out[i] = ints_out[0];
      }
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 1544727346:
      ints_out[0] = set_stopping_condition_minimum_density_parameter(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 1545587048:
      ints_out[0] = evolve_to_redshift(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 1552833343:
      for (int i = 0 ; i < call_count; i++){
        ints_out[i] = new_dm_particle(
          &ints_out[( 1 * call_count) + i] ,
          doubles_in[i] ,
          doubles_in[( 1 * call_count) + i] ,
          doubles_in[( 2 * call_count) + i] ,
          doubles_in[( 3 * call_count) + i] ,
          doubles_in[( 4 * call_count) + i] ,
          doubles_in[( 5 * call_count) + i] ,
          doubles_in[( 6 * call_count) + i]
        );
      }
      header_out[HEADER_INTEGER_COUNT] = 2 * call_count;
      break;
    
    case 1556511325:
      ints_out[0] = set_min_gas_temp(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 1558032891:
      ints_out[0] = get_err_tol_int_accuracy(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 1582712896:
      ints_out[0] = get_softening_halo_max_phys(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 1583464598:
      ints_out[0] = get_gamma(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 1616493870:
      ints_out[0] = set_energy_file(
        strings_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 1623514167:
      ints_out[0] = get_type_of_timestep_criterion(
        &ints_out[1]
      );
      header_out[HEADER_INTEGER_COUNT] = 2 * call_count;
      break;
    
    case 1623630901:
      for (int i = 0 ; i < call_count; i++){
        ints_out[i] = set_radius(
          ints_in[i] ,
          doubles_in[i]
        );
      }
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 1625942852:
      for (int i = 0 ; i < call_count; i++){
        ints_out[i] = has_stopping_condition(
          ints_in[i] ,
          &ints_out[( 1 * call_count) + i]
        );
      }
      header_out[HEADER_INTEGER_COUNT] = 2 * call_count;
      break;
    
    case 1626438772:
      ints_out[0] = set_bh_tol(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 1633912880:
      ints_out[0] = set_softening_halo_max_phys(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 1640685562:
      ints_out[0] = set_interpret_kicks_as_feedback_flag(
        booleans_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 1644113439:
      ints_out[0] = cleanup_code();
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 1646011510:
      ints_out[0] = get_min_size_timestep(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 1653324494:
      ints_out[0] = get_thermal_energy(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 1655137210:
      ints_out[0] = set_stopping_condition_maximum_density_parameter(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 1700863285:
      ints_out[0] = set_comoving_integration_flag(
        booleans_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 1706735752:
      ints_out[0] = internal__activate_communicator(
        ints_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 1744145122:
      ints_out[0] = recommit_parameters();
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 1768994498:
      ints_out[0] = initialize_code();
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 1822251412:
      ints_out[0] = get_redshift(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 1829880540:
      ints_out[0] = internal__open_port(
        &strings_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_STRING_COUNT] = 1 * call_count;
      break;
    
    case 1836308614:
      ints_out[0] = set_tree_domain_update_frequency(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 1851305433:
      ints_out[0] = get_pressure(
        &ints_in[0] ,
        &doubles_out[0] ,
        call_count
      );
      for (int i = 1 ; i < call_count; i++){
        ints_out[i] = ints_out[0];
      }
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 1852958273:
      ints_out[0] = get_potential_energy(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 1859127366:
      ints_out[0] = get_min_gas_temp(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 1861923976:
      ints_out[0] = get_state_sph(
        &ints_in[0] ,
        &doubles_out[0] ,
        &doubles_out[( 1 * call_count)] ,
        &doubles_out[( 2 * call_count)] ,
        &doubles_out[( 3 * call_count)] ,
        &doubles_out[( 4 * call_count)] ,
        &doubles_out[( 5 * call_count)] ,
        &doubles_out[( 6 * call_count)] ,
        &doubles_out[( 7 * call_count)] ,
        call_count
      );
      for (int i = 1 ; i < call_count; i++){
        ints_out[i] = ints_out[0];
      }
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 8 * call_count;
      break;
    
    case 1862635568:
      ints_out[0] = get_nogravity(
        &ints_out[1]
      );
      header_out[HEADER_INTEGER_COUNT] = 2 * call_count;
      break;
    
    case 1877555269:
      for (int i = 0 ; i < call_count; i++){
        ints_out[i] = get_gravity_at_point(
          doubles_in[i] ,
          doubles_in[( 1 * call_count) + i] ,
          doubles_in[( 2 * call_count) + i] ,
          doubles_in[( 3 * call_count) + i] ,
          &doubles_out[i] ,
          &doubles_out[( 1 * call_count) + i] ,
          &doubles_out[( 2 * call_count) + i]
        );
      }
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 3 * call_count;
      break;
    
    case 1891888088:
      ints_out[0] = set_redshift_max(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 1891926464:
      ints_out[0] = internal__accept_on_port(
        strings_in[0] ,
        &ints_out[1]
      );
      header_out[HEADER_INTEGER_COUNT] = 2 * call_count;
      break;
    
    case 1892689129:
      ints_out[0] = get_velocity(
        &ints_in[0] ,
        &doubles_out[0] ,
        &doubles_out[( 1 * call_count)] ,
        &doubles_out[( 2 * call_count)] ,
        call_count
      );
      for (int i = 1 ; i < call_count; i++){
        ints_out[i] = ints_out[0];
      }
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 3 * call_count;
      break;
    
    case 1895641480:
      ints_out[0] = get_min_gas_hsmooth_fractional(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 1911756772:
      ints_out[0] = set_hubble_param(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 1918746284:
      ints_out[0] = set_info_file(
        strings_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 1937183958:
      ints_out[0] = get_stopping_condition_out_of_box_parameter(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 1938857737:
      ints_out[0] = set_softening_gas_max_phys(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 1952713508:
      ints_out[0] = get_nsmtol(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 1979728989:
      ints_out[0] = get_time_max(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 2006446212:
      ints_out[0] = get_max_size_timestep(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 2010900811:
      ints_out[0] = get_position(
        &ints_in[0] ,
        &doubles_out[0] ,
        &doubles_out[( 1 * call_count)] ,
        &doubles_out[( 2 * call_count)] ,
        call_count
      );
      for (int i = 1 ; i < call_count; i++){
        ints_out[i] = ints_out[0];
      }
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 3 * call_count;
      break;
    
    case 2016681522:
      ints_out[0] = set_stopping_condition_maximum_internal_energy_parameter(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 2026192840:
      ints_out[0] = set_position(
        &ints_in[0] ,
        &doubles_in[0] ,
        &doubles_in[( 1 * call_count)] ,
        &doubles_in[( 2 * call_count)] ,
        call_count
      );
      for (int i = 1 ; i < call_count; i++){
        ints_out[i] = ints_out[0];
      }
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 2032046814:
      ints_out[0] = set_time_max(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 2033880761:
      ints_out[0] = set_time_between_statistics(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 2039500071:
      ints_out[0] = get_eps_is_h(
        &ints_out[1]
      );
      header_out[HEADER_INTEGER_COUNT] = 2 * call_count;
      break;
    
    case 2046699524:
      for (int i = 0 ; i < call_count; i++){
        ints_out[i] = get_stopping_condition_info(
          ints_in[i] ,
          &ints_out[( 1 * call_count) + i] ,
          &ints_out[( 2 * call_count) + i]
        );
      }
      header_out[HEADER_INTEGER_COUNT] = 3 * call_count;
      break;
    
    case 2061473599:
      ints_out[0] = get_acceleration(
        &ints_in[0] ,
        &doubles_out[0] ,
        &doubles_out[( 1 * call_count)] ,
        &doubles_out[( 2 * call_count)] ,
        call_count
      );
      for (int i = 1 ; i < call_count; i++){
        ints_out[i] = ints_out[0];
      }
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 3 * call_count;
      break;
    
    case 2064347450:
      ints_out[0] = set_unit_velocity(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 2069478464:
      ints_out[0] = commit_parameters();
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 2112574469:
      ints_out[0] = set_omega_baryon(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 2129795713:
      for (int i = 0 ; i < call_count; i++){
        ints_out[i] = get_stopping_condition_particle_index(
          ints_in[i] ,
          ints_in[( 1 * call_count) + i] ,
          &ints_out[( 1 * call_count) + i]
        );
      }
      header_out[HEADER_INTEGER_COUNT] = 2 * call_count;
      break;
    
    case 2143800670:
      ints_out[0] = set_timings_file(
        strings_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 2144268908:
      ints_out[0] = set_velocity(
        &ints_in[0] ,
        &doubles_in[0] ,
        &doubles_in[( 1 * call_count)] ,
        &doubles_in[( 2 * call_count)] ,
        call_count
      );
      for (int i = 1 ; i < call_count; i++){
        ints_out[i] = ints_out[0];
      }
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    default:
      header_out[HEADER_FLAGS] = header_out[HEADER_FLAGS] | ERROR_FLAG;
      strings_out[0] = new char[100];
      sprintf(strings_out[0], "unknown function id: %d\n", header_in[HEADER_FUNCTION_ID]);
      fprintf(stderr, "unknown function id: %d\n", header_in[HEADER_FUNCTION_ID]);
      header_out[HEADER_STRING_COUNT] = 1;
  }
  return true;
}

void onexit_mpi(void) {
#ifndef NOMPI
    int flag = 0;
    MPI_Finalized(&flag);
    
    if(!flag) {
        MPI_Comm parent;
        MPI_Comm_get_parent(&parent);
        
        int rank = 0;
        
        MPI_Comm_rank(parent, &rank);
        
        header_out[HEADER_FLAGS] = ERROR_FLAG;

        header_out[HEADER_CALL_ID] = 0;
        header_out[HEADER_FUNCTION_ID] = 0;
        header_out[HEADER_CALL_COUNT] = 0;
        header_out[HEADER_INTEGER_COUNT] = 0;
        header_out[HEADER_LONG_COUNT] = 0;
        header_out[HEADER_FLOAT_COUNT] = 0;
        header_out[HEADER_DOUBLE_COUNT] = 0;
        header_out[HEADER_BOOLEAN_COUNT] = 0;
        header_out[HEADER_STRING_COUNT] = 0;

        MPI_Send(header_out, HEADER_SIZE, MPI_INT, 0, 999, parent);
        
        for(int i = 0; i < lastid + 1; i++) {
            MPI_Comm_disconnect(&communicators[i]);
        }
        
        MPI_Finalize();
    }
#endif
}

void onexit_sockets(void) {
#ifdef WIN32
	closesocket(socketfd);
#else
	close(socketfd);
#endif
}

void send_array_sockets(void *buffer, int length, int file_descriptor, int rank) {
    int total_written = 0;
    int bytes_written;

    if (rank != 0) {
        return;
    }

    while (total_written < length) {
    
#ifdef WIN32
        bytes_written = send(file_descriptor, ((char *) buffer) + total_written,
                        length - total_written, 0);
#else
        bytes_written = write(file_descriptor, ((char *) buffer) + total_written,
                        length - total_written);
#endif


        if (bytes_written == -1) {
            perror("could not write data");
            exit(1);
        }

        total_written = total_written + bytes_written;
    }
}

void receive_array_sockets(void *buffer, int length, int file_descriptor, int rank) {
    int total_read = 0;
    int bytes_read;

    if (rank != 0) {
        return;
    }

    while (total_read < length) {
    
#ifdef WIN32
        bytes_read = recv(file_descriptor, ((char *) buffer) + total_read,
                        length - total_read, 0);
#else
        bytes_read = read(file_descriptor, ((char *) buffer) + total_read,
                        length - total_read);
#endif

        if (bytes_read == -1) {
            perror("could not read data");
            exit(1);
        }

        total_read = total_read + bytes_read;
    }
}

void new_arrays(int max_call_count) {
  ints_in = new int[ max_call_count * MAX_INTS_IN];
  ints_out = new int[ max_call_count * MAX_INTS_OUT];

  longs_in = new long long int[ max_call_count * MAX_LONGS_IN];
  longs_out = new long long int[ max_call_count * MAX_LONGS_OUT];

  floats_in = new float[ max_call_count * MAX_FLOATS_IN];
  floats_out = new float[ max_call_count * MAX_FLOATS_OUT];

  doubles_in = new double[ max_call_count * MAX_DOUBLES_IN];
  doubles_out = new double[ max_call_count * MAX_DOUBLES_OUT];
  
  booleans_in = new int[ max_call_count * MAX_BOOLEANS_IN];
  booleans_out = new int[ max_call_count * MAX_BOOLEANS_OUT];
  
  string_sizes_in = new int[ max_call_count * MAX_STRINGS_IN];
  string_sizes_out = new int[ max_call_count * MAX_STRINGS_OUT];

  strings_in = new char *[ max_call_count * MAX_STRINGS_IN];
  strings_out = new char *[ max_call_count * MAX_STRINGS_OUT];
}

void delete_arrays() {
  delete[] ints_in;
  delete[] ints_out;
  delete[] longs_in;
  delete[] longs_out;
  delete[] floats_in;
  delete[] floats_out;
  delete[] doubles_in;
  delete[] doubles_out;
  delete[] booleans_in;
  delete[] booleans_out;
  delete[] string_sizes_in;
  delete[] string_sizes_out;
  delete[] strings_in;
  delete[] strings_out;
}

void run_mpi(int argc, char *argv[]) {
#ifndef NOMPI
  int provided;
  int rank = 0;
  
  mpiIntercom = true;

  //fprintf(stderr, "C worker: running in mpi mode\n");
  
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  MPI_Comm parent;
  MPI_Comm_get_parent(&communicators[0]);
  lastid += 1;
  activeid = 0;
  parent = communicators[activeid];
  MPI_Comm_rank(parent, &rank);
  atexit(onexit_mpi);
  
  bool must_run_loop = true;
  
  int max_call_count = 10;
  
  header_in = new int[HEADER_SIZE];
  header_out = new int[HEADER_SIZE];

  new_arrays(max_call_count);  

  while(must_run_loop) {
    //fprintf(stderr, "receiving header\n");
    if(id_to_activate >= 0 && id_to_activate != activeid){
        activeid = id_to_activate;
        id_to_activate = -1;
        parent = communicators[activeid];
        MPI_Comm_rank(parent, &rank);
    }
    
    mpi_recv_header(parent);
    
    //fprintf(stderr, "C worker code: got header %d %d %d %d %d %d %d %d %d %d\n", header_in[0], header_in[1], header_in[2], header_in[3], header_in[4], header_in[5], header_in[6], header_in[7], header_in[8], header_in[9]);
    
    int call_count = header_in[HEADER_CALL_COUNT];

    if (call_count > max_call_count) {
      delete_arrays();
      max_call_count = call_count + 255;
      new_arrays(max_call_count);
    }
    
    if(header_in[HEADER_INTEGER_COUNT] > 0) {
      MPI_Bcast(ints_in, header_in[HEADER_INTEGER_COUNT] , MPI_INT, 0, parent);
    }
    
    if(header_in[HEADER_LONG_COUNT] > 0) {
      MPI_Bcast(longs_in, header_in[HEADER_LONG_COUNT], MPI_LONG_LONG_INT, 0, parent);
    }
    
    if(header_in[HEADER_FLOAT_COUNT] > 0) {
      MPI_Bcast(floats_in, header_in[HEADER_FLOAT_COUNT], MPI_FLOAT, 0, parent);
    }
    
    if(header_in[HEADER_DOUBLE_COUNT] > 0) {
      MPI_Bcast(doubles_in, header_in[HEADER_DOUBLE_COUNT], MPI_DOUBLE, 0, parent);
    }
    
    if(header_in[HEADER_BOOLEAN_COUNT] > 0) {
      MPI_Bcast(booleans_in, header_in[HEADER_BOOLEAN_COUNT], MPI_INTEGER, 0, parent);
    }
    
    if(header_in[HEADER_STRING_COUNT] > 0) {
      MPI_Bcast(string_sizes_in, header_in[HEADER_STRING_COUNT], MPI_INTEGER, 0, parent);
      
      int total_string_size = 0;
      for (int i = 0; i < header_in[HEADER_STRING_COUNT];i++) {
        total_string_size += string_sizes_in[i] + 1;
      }
      
      characters_in = new char[total_string_size];
      MPI_Bcast(characters_in, total_string_size, MPI_CHARACTER, 0, parent);

      int offset = 0;
      for (int i = 0 ; i <  header_in[HEADER_STRING_COUNT];i++) {
          strings_in[i] = characters_in + offset;
          offset += string_sizes_in[i] + 1;
      } 
    }

    header_out[HEADER_FLAGS] = 0;
    header_out[HEADER_CALL_ID] = header_in[HEADER_CALL_ID];
    header_out[HEADER_FUNCTION_ID] = header_in[HEADER_FUNCTION_ID];
    header_out[HEADER_CALL_COUNT] = call_count;
    header_out[HEADER_INTEGER_COUNT] = 0;
    header_out[HEADER_LONG_COUNT] = 0;
    header_out[HEADER_FLOAT_COUNT] = 0;
    header_out[HEADER_DOUBLE_COUNT] = 0;
    header_out[HEADER_BOOLEAN_COUNT] = 0;
    header_out[HEADER_STRING_COUNT] = 0;

    //fprintf(stderr, "c worker mpi: handling call\n");
    
    must_run_loop = handle_call();
    
    //fprintf(stderr, "c worker mpi: call handled\n");
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if(rank == 0) {
      MPI_Send(header_out, HEADER_SIZE, MPI_INT, 0, 999, parent);
      
      if(header_out[HEADER_INTEGER_COUNT] > 0) {
        MPI_Send(ints_out, header_out[HEADER_INTEGER_COUNT], MPI_INT, 0, 999, parent);
      }
      if(header_out[HEADER_LONG_COUNT] > 0) {
        MPI_Send(longs_out, header_out[HEADER_LONG_COUNT], MPI_LONG_LONG_INT, 0, 999, parent);
      }
      if(header_out[HEADER_FLOAT_COUNT] > 0) {
        MPI_Send(floats_out, header_out[HEADER_FLOAT_COUNT], MPI_FLOAT, 0, 999, parent);
      }
      if(header_out[HEADER_DOUBLE_COUNT] > 0) {
        MPI_Send(doubles_out, header_out[HEADER_DOUBLE_COUNT], MPI_DOUBLE, 0, 999, parent);
      }
      if(header_out[HEADER_BOOLEAN_COUNT] > 0) {
        MPI_Send(booleans_out, header_out[HEADER_BOOLEAN_COUNT], MPI_INTEGER, 0, 999, parent);
      }
      if(header_out[HEADER_STRING_COUNT] > 0) {
        int offset = 0;
        for( int i = 0; i < header_out[HEADER_STRING_COUNT] ; i++) {
          
          int length = strlen(strings_out[i]);
          string_sizes_out[i] = length;
          offset += length + 1;
        }
        
        characters_out = new char[offset + 1];
        offset = 0;
        
        for( int i = 0; i < header_out[HEADER_STRING_COUNT]  ; i++) {
          strcpy(characters_out+offset, strings_out[i]);
          offset += string_sizes_out[i] + 1;
        }
        MPI_Send(string_sizes_out, header_out[HEADER_STRING_COUNT], MPI_INTEGER, 0, 999, parent);
        MPI_Send(characters_out, offset, MPI_BYTE, 0, 999, parent);
      }
    
    }
    
    if (characters_in) { 
        delete[] characters_in;
        characters_in = 0;
    }
    
    if (characters_out) {
        delete[] characters_out;
        characters_out = 0;
    }
    //fprintf(stderr, "call done\n");
  }
  delete_arrays();
  
    for(int i = 0; i < lastid + 1; i++) {
        MPI_Comm_disconnect(&communicators[i]);
    }
    
    MPI_Finalize();
  //fprintf(stderr, "mpi finalized\n");
#else
  fprintf(stderr, "mpi support not compiled into worker\n");
  exit(1);
#endif
}

void run_sockets_mpi(int argc, char *argv[], int port, char *host) {
#ifndef NOMPI
 bool must_run_loop = true;
  int max_call_count = 10;
  struct sockaddr_in serv_addr;
  struct hostent *server;
  int on = 1;
  
  mpiIntercom = false;

  MPI::Init_thread(argc, argv, MPI_THREAD_MULTIPLE);
  int rank = MPI::COMM_WORLD.Get_rank();
  
  if (rank == 0) {  
    //fprintf(stderr, "C worker: running in sockets+mpi mode\n");
  
   
    socketfd = socket(AF_INET, SOCK_STREAM, 0);
    
    if (socketfd < 0) {
      perror("ERROR opening socket");
      //fprintf(stderr, "cannot open socket\n");
      exit(1);
    }

    //turn on no-delay option in tcp for huge speed improvement
    setsockopt (socketfd, IPPROTO_TCP, TCP_NODELAY, &on, sizeof (on));
    
    server = gethostbyname(host);
    
    memset((char *) &serv_addr, '\0', sizeof(serv_addr));
    serv_addr.sin_family = AF_INET;
    memcpy((char *) &serv_addr.sin_addr.s_addr, (char *) server->h_addr, server->h_length);
    serv_addr.sin_port = htons(port);
  
    if (connect(socketfd, (struct sockaddr *) &serv_addr, sizeof(serv_addr)) < 0) {
      perror("ERROR connecting socket");
      //fprintf(stderr, "cannot connect socket\n");
      exit(1);
    }
    
    //fprintf(stderr, "sockets_mpi: finished initializing code\n");
  
    atexit(onexit_sockets);
  
  }
  
  header_in = new int[HEADER_SIZE];
  header_out = new int[HEADER_SIZE];

  new_arrays(max_call_count);  
  
  while(must_run_loop) {
    //fprintf(stderr, "sockets_mpi: receiving header\n");
    receive_array_sockets(header_in, HEADER_SIZE * sizeof(int), socketfd, rank);
    MPI::COMM_WORLD.Bcast(header_in, HEADER_SIZE, MPI_INT, 0);
    
    //fprintf(stderr, "C sockets_mpi worker code: got header %d %d %d %d %d %d %d %d %d %d\n", header_in[0], header_in[1], header_in[2], header_in[3], header_in[4], header_in[5], header_in[6], header_in[7], header_in[8], header_in[9]);
    
    int call_count = header_in[HEADER_CALL_COUNT];

    if (call_count > max_call_count) {
      delete_arrays();
      max_call_count = call_count + 255;
      new_arrays(max_call_count);
    }
    
    if (header_in[HEADER_INTEGER_COUNT] > 0) {
      receive_array_sockets(ints_in, header_in[HEADER_INTEGER_COUNT] * sizeof(int), socketfd, rank);
      MPI::COMM_WORLD.Bcast(ints_in, header_in[HEADER_INTEGER_COUNT], MPI_INTEGER, 0);
    }
     
    if (header_in[HEADER_LONG_COUNT] > 0) {
      receive_array_sockets(longs_in, header_in[HEADER_LONG_COUNT] * sizeof(long long int), socketfd, rank);
      MPI::COMM_WORLD.Bcast(longs_in, header_in[HEADER_LONG_COUNT], MPI_LONG_LONG_INT, 0);
    }
    
    if(header_in[HEADER_FLOAT_COUNT] > 0) {
      receive_array_sockets(floats_in, header_in[HEADER_FLOAT_COUNT] * sizeof(float), socketfd, rank);
      MPI::COMM_WORLD.Bcast(floats_in, header_in[HEADER_FLOAT_COUNT], MPI_FLOAT, 0);
    }
    
    if(header_in[HEADER_DOUBLE_COUNT] > 0) {
      receive_array_sockets(doubles_in, header_in[HEADER_DOUBLE_COUNT] * sizeof(double), socketfd, rank);
      MPI::COMM_WORLD.Bcast(doubles_in, header_in[HEADER_DOUBLE_COUNT], MPI_DOUBLE, 0);
    }
    
    if(header_in[HEADER_BOOLEAN_COUNT] > 0) {
      //receive_array_sockets(booleans_in, header_in[HEADER_BOOLEAN_COUNT], socketfd , rank);
      for (int i = 0; i < header_in[HEADER_BOOLEAN_COUNT]; i++) {
        booleans_in[i] = 0;
        receive_array_sockets(&booleans_in[i], 1, socketfd , rank);
      }
      MPI::COMM_WORLD.Bcast(booleans_in, header_in[HEADER_BOOLEAN_COUNT], MPI_INTEGER, 0);
    }
    
    if(header_in[HEADER_STRING_COUNT] > 0) {
      receive_array_sockets(string_sizes_in, header_in[HEADER_STRING_COUNT] * sizeof(int), socketfd, rank);
      MPI::COMM_WORLD.Bcast(string_sizes_in, header_in[HEADER_STRING_COUNT], MPI_INT, 0);
      for (int i = 0; i < header_in[HEADER_STRING_COUNT]; i++) {
        strings_in[i] = new char[string_sizes_in[i] + 1];
        receive_array_sockets(strings_in[i], string_sizes_in[i], socketfd, rank);
        MPI::COMM_WORLD.Bcast(strings_in[i], string_sizes_in[i], MPI_CHARACTER, 0);
        strings_in[i][string_sizes_in[i]] = '\0';
      }
    }
    
    header_out[HEADER_FLAGS] = 0;
    header_out[HEADER_CALL_ID] = header_in[HEADER_CALL_ID];
    header_out[HEADER_FUNCTION_ID] = header_in[HEADER_FUNCTION_ID];
    header_out[HEADER_CALL_COUNT] = call_count;
    header_out[HEADER_INTEGER_COUNT] = 0;
    header_out[HEADER_LONG_COUNT] = 0;
    header_out[HEADER_FLOAT_COUNT] = 0;
    header_out[HEADER_DOUBLE_COUNT] = 0;
    header_out[HEADER_BOOLEAN_COUNT] = 0;
    header_out[HEADER_STRING_COUNT] = 0;

    //fprintf(stderr, "c worker sockets_mpi: handling call\n");
    
    must_run_loop = handle_call();
    
    //fprintf(stderr, "c worker sockets_mpi: call handled\n");
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (rank == 0) {

      send_array_sockets(header_out, HEADER_SIZE * sizeof(int), socketfd, 0);
          
      if(header_out[HEADER_INTEGER_COUNT] > 0) {
        send_array_sockets(ints_out, header_out[HEADER_INTEGER_COUNT] * sizeof(int), socketfd, 0);
      }
          
      if(header_out[HEADER_LONG_COUNT] > 0) {
        send_array_sockets(longs_out, header_out[HEADER_LONG_COUNT] * sizeof(long long int), socketfd, 0);
      }
          
      if(header_out[HEADER_FLOAT_COUNT] > 0) {
        send_array_sockets(floats_out, header_out[HEADER_FLOAT_COUNT] * sizeof(float), socketfd, 0);
      }
          
      if(header_out[HEADER_DOUBLE_COUNT] > 0) {
        send_array_sockets(doubles_out, header_out[HEADER_DOUBLE_COUNT] * sizeof(double), socketfd, 0);
      }
          
      if(header_out[HEADER_BOOLEAN_COUNT] > 0) {
          for (int i = 0; i < header_out[HEADER_BOOLEAN_COUNT]; i++) {
            if (booleans_out[i]) {
              send_array_sockets(&TRUE_BYTE, 1, socketfd, 0);
            } else {
              send_array_sockets(&FALSE_BYTE, 1, socketfd, 0);
            }
         }
      }
          
      if(header_out[HEADER_STRING_COUNT] > 0) {
        for (int i = 0; i < header_out[HEADER_STRING_COUNT]; i++) {
          string_sizes_out[i] = strlen(strings_out[i]);
        }
        send_array_sockets(string_sizes_out, header_out[HEADER_STRING_COUNT] * sizeof(int), socketfd, 0);
          
        for (int i = 0; i < header_out[HEADER_STRING_COUNT]; i++) {
          send_array_sockets(strings_out[i], string_sizes_out[i] * sizeof(char), socketfd, 0);
        }
      }
        
      //fprintf(stderr, "sockets_mpicall done\n");
    }

  }
  delete_arrays();
  
  if (rank == 0) {
  
#ifdef WIN32
	closesocket(socketfd);
#else
	close(socketfd);
#endif
  }
  
  MPI_Finalize();
  
  //fprintf(stderr, "sockets_mpi done\n");
#else
  fprintf(stderr, "mpi support not compiled into worker\n");
  exit(1);
#endif
}

void run_sockets(int port, char *host) {
  bool must_run_loop = true;
  int max_call_count = 10;
  struct sockaddr_in serv_addr;
  struct hostent *server;
  int on = 1;

#ifdef WIN32
	WSADATA wsaData;
	int iResult;

	// Initialize Winsock
	iResult = WSAStartup(MAKEWORD(2,2), &wsaData);
	if (iResult != 0) {
	printf("WSAStartup failed: %d\n", iResult);
	exit(1);
	}
#endif
  mpiIntercom = false;

  //fprintf(stderr, "C worker: running in sockets mode\n");
   
  socketfd = socket(AF_INET, SOCK_STREAM, 0);
    
  if (socketfd < 0) {
    fprintf(stderr, "cannot open socket\n");
    exit(1);
  }
  
  //turn on no-delay option in tcp for huge speed improvement
  setsockopt (socketfd, IPPROTO_TCP, TCP_NODELAY, (const char *)&on, sizeof (on));
    
  server = gethostbyname(host);
    
  memset((char *) &serv_addr, '\0', sizeof(serv_addr));
  serv_addr.sin_family = AF_INET;
  memcpy((char *) &serv_addr.sin_addr.s_addr, (char *) server->h_addr, server->h_length);
  serv_addr.sin_port = htons(port);
  
  if (connect(socketfd, (struct sockaddr *) &serv_addr, sizeof(serv_addr)) < 0) {
    fprintf(stderr, "cannot connect socket to host %s, port %d\n", host, port);

    perror("ERROR connecting socket");
    //fprintf(stderr, "cannot connect socket\n");
    exit(1);
  }
    
  //fprintf(stderr, "sockets: finished initializing code\n");
  
  atexit(onexit_sockets);
  
  header_in = new int[HEADER_SIZE];
  header_out = new int[HEADER_SIZE];

  new_arrays(max_call_count);  
  
  while(must_run_loop) {
    //fprintf(stderr, "sockets: receiving header\n");
    receive_array_sockets(header_in, HEADER_SIZE * sizeof(int), socketfd, 0);
    //fprintf(stderr, "C sockets worker code: got header %d %d %d %d %d %d %d %d %d %d\n", header_in[0], header_in[1], header_in[2], header_in[3], header_in[4], header_in[5], header_in[6], header_in[7], header_in[8], header_in[9]);
    
    int call_count = header_in[HEADER_CALL_COUNT];

    if (call_count > max_call_count) {
      delete_arrays();
      max_call_count = call_count + 255;
      new_arrays(max_call_count);
    }
    
    if (header_in[HEADER_INTEGER_COUNT] > 0) {
      receive_array_sockets(ints_in, header_in[HEADER_INTEGER_COUNT] * sizeof(int), socketfd, 0);
    }
     
    if (header_in[HEADER_LONG_COUNT] > 0) {
      receive_array_sockets(longs_in, header_in[HEADER_LONG_COUNT] * sizeof(long long int), socketfd, 0);
    }
    
    if(header_in[HEADER_FLOAT_COUNT] > 0) {
      receive_array_sockets(floats_in, header_in[HEADER_FLOAT_COUNT] * sizeof(float), socketfd, 0);
    }
    
    if(header_in[HEADER_DOUBLE_COUNT] > 0) {
      receive_array_sockets(doubles_in, header_in[HEADER_DOUBLE_COUNT] * sizeof(double), socketfd, 0);
    }
    
    if(header_in[HEADER_BOOLEAN_COUNT] > 0) {
      //receive_array_sockets(booleans_in, header_in[HEADER_BOOLEAN_COUNT], socketfd , 0);
      for (int i = 0; i < header_in[HEADER_BOOLEAN_COUNT]; i++) {
        booleans_in[i] = 0;
        receive_array_sockets(&booleans_in[i], 1, socketfd , 0);
      }
    }
    
    if(header_in[HEADER_STRING_COUNT] > 0) {
      receive_array_sockets(string_sizes_in, header_in[HEADER_STRING_COUNT] * sizeof(int), socketfd, 0);
      for (int i = 0; i < header_in[HEADER_STRING_COUNT]; i++) {
        strings_in[i] = new char[string_sizes_in[i] + 1];
        receive_array_sockets(strings_in[i], string_sizes_in[i], socketfd, 0);
        strings_in[i][string_sizes_in[i]] = '\0';
      }
    }
    
    header_out[HEADER_FLAGS] = 0;
    header_out[HEADER_CALL_ID] = header_in[HEADER_CALL_ID];
    header_out[HEADER_FUNCTION_ID] = header_in[HEADER_FUNCTION_ID];
    header_out[HEADER_CALL_COUNT] = call_count;
    header_out[HEADER_INTEGER_COUNT] = 0;
    header_out[HEADER_LONG_COUNT] = 0;
    header_out[HEADER_FLOAT_COUNT] = 0;
    header_out[HEADER_DOUBLE_COUNT] = 0;
    header_out[HEADER_BOOLEAN_COUNT] = 0;
    header_out[HEADER_STRING_COUNT] = 0;

    //fprintf(stderr, "c worker sockets: handling call\n");
    
    must_run_loop = handle_call();
    
    //fprintf(stderr, "c worker sockets: call handled\n");

    send_array_sockets(header_out, HEADER_SIZE * sizeof(int), socketfd, 0);
      
    if(header_out[HEADER_INTEGER_COUNT] > 0) {
      send_array_sockets(ints_out, header_out[HEADER_INTEGER_COUNT] * sizeof(int), socketfd, 0);
    }
      
    if(header_out[HEADER_LONG_COUNT] > 0) {
      send_array_sockets(longs_out, header_out[HEADER_LONG_COUNT] * sizeof(long long int), socketfd, 0);
    }
      
    if(header_out[HEADER_FLOAT_COUNT] > 0) {
      send_array_sockets(floats_out, header_out[HEADER_FLOAT_COUNT] * sizeof(float), socketfd, 0);
    }
      
    if(header_out[HEADER_DOUBLE_COUNT] > 0) {
      send_array_sockets(doubles_out, header_out[HEADER_DOUBLE_COUNT] * sizeof(double), socketfd, 0);
    }
      
    if(header_out[HEADER_BOOLEAN_COUNT] > 0) {
        for (int i = 0; i < header_out[HEADER_BOOLEAN_COUNT]; i++) {
          if (booleans_out[i]) {
            send_array_sockets(&TRUE_BYTE, 1, socketfd, 0);
          } else {
            send_array_sockets(&FALSE_BYTE, 1, socketfd, 0);
          }
       }
    }
      
    if(header_out[HEADER_STRING_COUNT] > 0) {
      for (int i = 0; i < header_out[HEADER_STRING_COUNT]; i++) {
        string_sizes_out[i] = strlen(strings_out[i]);
      }
      send_array_sockets(string_sizes_out, header_out[HEADER_STRING_COUNT] * sizeof(int), socketfd, 0);
      
      for (int i = 0; i < header_out[HEADER_STRING_COUNT]; i++) {
        send_array_sockets(strings_out[i], string_sizes_out[i] * sizeof(char), socketfd, 0);
      }
    }
    
    //fprintf(stderr, "call done\n");
  }
  delete_arrays();
  
#ifdef WIN32
	closesocket(socketfd);
#else
	close(socketfd);
#endif
  //fprintf(stderr, "sockets done\n");
}
 
int main(int argc, char *argv[]) {
  int port;
  bool use_mpi;
  char *host;
  
  //for(int i = 0 ; i < argc; i++) {
  //  fprintf(stderr, "argument %d is %s\n", i, argv[i]);
  //}

  if (argc == 1) {
    run_mpi(argc, argv);
  } else if (argc == 4) {
    port = atoi(argv[1]);
    host = argv[2];
    
    if (strcmp(argv[3], "true") == 0) {
      use_mpi = true;
    } else if (strcmp(argv[3], "false") == 0) {
      use_mpi = false;
    } else {
      fprintf(stderr, "mpi enabled setting must be either 'true' or 'false', not %s\n", argv[2]);
      fprintf(stderr, "usage: %s [PORT HOST MPI_ENABLED]\n", argv[0]);
      exit(1);
    }    
    
    if (use_mpi) {
      run_sockets_mpi(argc, argv, port, host);
    } else {
      run_sockets(port, host);
    }
  } else {
    fprintf(stderr, "%s need either 0 or 4 arguments, not %d\n", argv[0], argc);
    fprintf(stderr, "usage: %s [PORT HOST MPI_ENABLED]\n", argv[0]);
    exit(1);
  }

  return 0;
}   


