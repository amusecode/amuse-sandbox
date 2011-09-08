#include <stdio.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <mpi.h>
#include <unistd.h>
#include <stdlib.h>
#include <netdb.h>
#include <string.h>

#include "amuse-rpc.h"

int socketfd = 0;

void log(char *message) {
	fprintf(stderr, "socketrpc: %s\n", message);
}

void error(char *message) {
	fprintf(stderr, "socketrpc: ERROR: %s\n", message);
}

void socketrpc_send(void *buffer, int length, int file_descriptor) {
    int total_written = 0;
    int bytes_written;

    while (total_written < length) {
        bytes_written = write(file_descriptor, ((char *) buffer) + total_written,
                        length - total_written);

        if (bytes_written == -1) {
            error("could not write data");
            exit(1);
        }

        total_written = total_written + bytes_written;
    }
}

void socketrpc_receive(void *buffer, int length, int file_descriptor) {
    int total_read = 0;
    int bytes_read;

    while (total_read < length) {
        bytes_read = read(file_descriptor, ((char *) buffer) + total_read,
                        length - total_read);

        if (bytes_read == -1) {
            error("could not read data");
            exit(1);
        }

        total_read = total_read + bytes_read;
    }
}

//handle calls one by one until a "termination" call is received
void handle_calls() {
//while (not finished) {
//receive message
//broadcast mpi
//call code
//send result
//end while
}

int main(int argc, char *argv[]) {
	int portno;
	struct sockaddr_in serv_addr;
	struct hostent *server;
	int size;
	int rank;

	log("initializing MPI");

	//init MPI (if applicable)
	//MPI::Init_thread(argc, argv, MPI_THREAD_MULTIPLE);
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	log("connecting to socket");

	//connect socket to port
	portno = atoi(argv[1]);
	socketfd = socket(AF_INET, SOCK_STREAM, 0);

	if (socketfd < 0) {
		error("cannot open socket");
		exit(1);
	}

	server = gethostbyname("localhost");

	bzero((char *) &serv_addr, sizeof(serv_addr));
	serv_addr.sin_family = AF_INET;
	bcopy((char *) server->h_addr, (char *) &serv_addr.sin_addr.s_addr,
					server->h_length);
	serv_addr.sin_port = htons(portno);

	if (connect(socketfd, (struct sockaddr *) &serv_addr, sizeof(serv_addr)) < 0) {
		error("cannot connect socket");
		exit(1);
	}

	log("initialization complete, now handling calls");

	handle_calls();

	log("closing socket");

	close(socketfd);

	log("finalizing MPI");

	MPI_Finalize();

	return 0;
}
