
#ifndef NOMPI
    #include <mpi.h>
#endif
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <unistd.h>
#include <netinet/tcp.h>

#include "src/code.h"

static bool NEEDS_MPI = true;

static int MAX_INTS_IN = 1;
static int MAX_INTS_OUT = 3;
static int MAX_STRINGS_IN = 0;
static int MAX_STRINGS_OUT = 0;
static int MAX_DOUBLES_IN = 14;
static int MAX_DOUBLES_OUT = 13;
static int MAX_LONGS_IN = 0;
static int MAX_LONGS_OUT = 0;
static int MAX_BOOLEANS_IN = 0;
static int MAX_BOOLEANS_OUT = 0;
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


bool handle_call() {
  int call_count = header_in[HEADER_CALL_COUNT];
  
  switch(header_in[HEADER_FUNCTION_ID]) {
    case 0:
      return false;
      break;
    case 61072875:
      ints_out[0] = set_f_25PN_2(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 111364973:
      ints_out[0] = evolve(
        doubles_in[0] ,
        doubles_in[1] ,
        doubles_in[2] ,
        doubles_in[3] ,
        doubles_in[4] ,
        doubles_in[5] ,
        doubles_in[6] ,
        doubles_in[7] ,
        doubles_in[8] ,
        doubles_in[9] ,
        doubles_in[10] ,
        doubles_in[11] ,
        doubles_in[12] ,
        doubles_in[13] ,
        &doubles_out[0] ,
        &doubles_out[1] ,
        &doubles_out[2] ,
        &doubles_out[3] ,
        &doubles_out[4] ,
        &doubles_out[5] ,
        &doubles_out[6] ,
        &doubles_out[7] ,
        &doubles_out[8] ,
        &doubles_out[9] ,
        &doubles_out[10] ,
        &doubles_out[11] ,
        &doubles_out[12] ,
        &ints_out[1] ,
        &ints_out[2]
      );
      header_out[HEADER_INTEGER_COUNT] = 3 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 13 * call_count;
      break;
    
    case 213726058:
      ints_out[0] = get_f_25PN_2(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 675097567:
      ints_out[0] = set_f_1PN_2(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 727361823:
      ints_out[0] = internal__set_message_polling_interval(
        ints_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 907189400:
      ints_out[0] = get_f_1PN_2(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 919768251:
      ints_out[0] = internal__get_message_polling_interval(
        &ints_out[1]
      );
      header_out[HEADER_INTEGER_COUNT] = 2 * call_count;
      break;
    
    case 1321960859:
      ints_out[0] = set_f_1PN_1(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 1357132498:
      ints_out[0] = get_f_1PN_1(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 1700088751:
      ints_out[0] = set_f_25PN_1(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 1783340332:
      ints_out[0] = get_f_25PN_1(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 1912826050:
      ints_out[0] = set_f_oct(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 1918429639:
      ints_out[0] = get_f_1PN_12(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 1968407667:
      ints_out[0] = get_relative_tolerance(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 2094040131:
      ints_out[0] = get_f_oct(
        &doubles_out[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      header_out[HEADER_DOUBLE_COUNT] = 1 * call_count;
      break;
    
    case 2101755716:
      ints_out[0] = set_f_1PN_12(
        doubles_in[0]
      );
      header_out[HEADER_INTEGER_COUNT] = 1 * call_count;
      break;
    
    case 2147132108:
      ints_out[0] = set_relative_tolerance(
        doubles_in[0]
      );
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
        
        MPI_Comm_disconnect(&parent);
        
        MPI_Finalize();
    }
#endif
}

void onexit_sockets(void) {
    close(socketfd);
}

void send_array_sockets(void *buffer, int length, int file_descriptor, int rank) {
    int total_written = 0;
    int bytes_written;

    if (rank != 0) {
        return;
    }

    while (total_written < length) {
        bytes_written = write(file_descriptor, ((char *) buffer) + total_written,
                        length - total_written);

        if (bytes_written == -1) {
            fprintf(stderr, "could not write data\n");
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
        bytes_read = read(file_descriptor, ((char *) buffer) + total_read,
                        length - total_read);

        if (bytes_read == -1) {
            fprintf(stderr, "could not read data\n");
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
  MPI_Comm parent;
  int rank = 0;
  
  mpiIntercom = true;

  //fprintf(stderr, "C worker: running in mpi mode\n");
  
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  MPI_Comm_get_parent(&parent);
  MPI_Comm_rank(parent, &rank);
  atexit(onexit_mpi);
  
  bool must_run_loop = true;
  
  int max_call_count = 10;
  
  header_in = new int[HEADER_SIZE];
  header_out = new int[HEADER_SIZE];

  new_arrays(max_call_count);  

  while(must_run_loop) {
    //fprintf(stderr, "receiving header\n");
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
  
  MPI_Comm_disconnect(&parent);
  MPI_Finalize();
  //fprintf(stderr, "mpi finalized\n");
#else
  fprintf(stderr, "mpi support not compiled into worker\n");
  exit(1);
#endif
}

void run_sockets_mpi(int argc, char *argv[], int port) {
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
      fprintf(stderr, "cannot open socket\n");
      exit(1);
    }

    //turn on no-delay option in tcp for huge speed improvement
    setsockopt (socketfd, IPPROTO_TCP, TCP_NODELAY, &on, sizeof (on));
    
    server = gethostbyname("localhost");
    
    bzero((char *) &serv_addr, sizeof(serv_addr));
    serv_addr.sin_family = AF_INET;
    bcopy((char *) server->h_addr, (char *) &serv_addr.sin_addr.s_addr, server->h_length);
    serv_addr.sin_port = htons(port);
  
    if (connect(socketfd, (struct sockaddr *) &serv_addr, sizeof(serv_addr)) < 0) {
      fprintf(stderr, "cannot connect socket\n");
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
    close(socketfd);
  }
  
  MPI_Finalize();
  
  //fprintf(stderr, "sockets_mpi done\n");
#else
  fprintf(stderr, "mpi support not compiled into worker\n");
  exit(1);
#endif
}

void run_sockets(int port) {
  bool must_run_loop = true;
  int max_call_count = 10;
  struct sockaddr_in serv_addr;
  struct hostent *server;
  int on = 1;
  
  mpiIntercom = false;

  //fprintf(stderr, "C worker: running in sockets mode\n");
   
  socketfd = socket(AF_INET, SOCK_STREAM, 0);
    
  if (socketfd < 0) {
    fprintf(stderr, "cannot open socket\n");
    exit(1);
  }
  
  //turn on no-delay option in tcp for huge speed improvement
  setsockopt (socketfd, IPPROTO_TCP, TCP_NODELAY, &on, sizeof (on));
    
  server = gethostbyname("localhost");
    
  bzero((char *) &serv_addr, sizeof(serv_addr));
  serv_addr.sin_family = AF_INET;
  bcopy((char *) server->h_addr, (char *) &serv_addr.sin_addr.s_addr, server->h_length);
  serv_addr.sin_port = htons(port);
  
  if (connect(socketfd, (struct sockaddr *) &serv_addr, sizeof(serv_addr)) < 0) {
    fprintf(stderr, "cannot connect socket\n");
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
  
  close(socketfd);
  //fprintf(stderr, "sockets done\n");
}
 
int main(int argc, char *argv[]) {
  int port;
  bool use_mpi = NEEDS_MPI;
  
  for(int i = 0 ; i < argc; i++) {
    //fprintf(stderr, "argument %d is %s\n", i, argv[i]);
  }

  if (argc == 1) {
    run_mpi(argc, argv);
  } else {
    port = atoi(argv[1]);
    
    if (argc >= 3) {
       if (strcmp(argv[2], "no_mpi") == 0) {
         use_mpi = false;
       } else if (strcmp(argv[2], "mpi") == 0) {
         use_mpi = true;
       } else if (strcmp(argv[2], "auto") == 0) {
         use_mpi = NEEDS_MPI;
       }

    }
    
    if (use_mpi) {
      run_sockets_mpi(argc, argv, port);
    } else {
      run_sockets(port);
    }
  }  

  return 0;
}   


