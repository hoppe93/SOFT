/* Various MPI helper functions */
#include <mpi.h>
#include <stdio.h>
#include <limits.h>

#include "smpi.h"

/**
 * Sends an 'output ready' signal to the next
 * MPI process.
 */
void smpi_sor(int tag) {
	int rank, nprocesses;
	MPI_Comm_size(MPI_COMM_WORLD, &nprocesses);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	/* Is this the last process in line? */
	if (rank == nprocesses-1) return;	/* Then we're done! */

	int msg = SMPI_OUTPUT_READY, err;
	if ((err=MPI_Send(&msg, 1, MPI_INT, rank+1, tag, MPI_COMM_WORLD)) != MPI_SUCCESS) {
		fprintf(stderr, "MPI ERROR %d: Failed to send MPI message 'output ready'. Aborting...\n", err);
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
}
/**
 * Waits for an 'output ready' signal from
 * process below (i.e. with rank = myRank-1).
 * The 'output ready' signal is sent by the
 * lower process when it has finished writing
 * to the output file, indicating that the file
 * is ready to receive data from the next MPI process.
 */
void smpi_wor(int tag) {
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	/* If this is the first process there's no use in waiting! */
	if (rank == 0) return;

	int msg, err;
	if ((err=MPI_Recv(&msg, 1, MPI_INT, rank-1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE)) != MPI_SUCCESS) {
		fprintf(stderr, "MPI ERROR %d: Failed to receive MPI message 'output ready'. Aborting...\n", err);
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
}

long long int smpi_receive_integer(int sender, int tag) {
	int err;
	long long int val;

	if ((err=MPI_Recv(&val, 1, MPI_LONG_LONG_INT, sender, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE))!=MPI_SUCCESS) {
		fprintf(stderr, "MPI ERROR %d: Failed to receive integer from process %d over MPI. Aborting...\n", err, sender);
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	return val;
}

/**
 * Receive a matrix (double array) from
 * the process with index 'destination' over MPI.
 *
 * matrix: Target to store matrix in.
 * elements: Number of elements to read.
 * sender: MPI ID of process sending the data.
 */
void smpi_receive_matrix(double *matrix, size_t elements, int sender, int tag) {
	int err;
	int chunksize = INT_MAX, toread;
	size_t read = 0;

	while (read < elements) {
		if (elements-read > (size_t)chunksize) toread = chunksize;
		else toread = elements-read;

		if ((err=MPI_Recv(matrix+read, toread, MPI_DOUBLE, sender, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE))!=MPI_SUCCESS) {
			fprintf(stderr, "MPI ERROR %d: Failed to receive matrix from process %d over MPI. Aborting...\n", err, sender);
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		read += toread;
	}
}

/**
 * Send integer value to process
 * with index 'destination' over MPI.
 *
 * value: Integer value to send
 * destination: Destination process ID
 * tag: MPI tag to use
 */
void smpi_send_integer(long long int value, int destination, int tag) {
	int err;

	if ((err=MPI_Send(&value, 1, MPI_LONG_LONG_INT, destination, tag, MPI_COMM_WORLD))!=MPI_SUCCESS) {
		fprintf(stderr, "MPI ERROR %d: Faield to send integer over MPI to process %d. Aborting...\n", err, destination);
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
}

/**
 * Send a matrix (double array) to process
 * with index 'destination' over MPI.
 *
 * matrix: Pointer to beginning of matrix. The matrix
 *   is assumed to be continuous in memory.
 * elements: Total number of elements in 'matrix'. First
 *   element is assumed to be matrix[0], last element
 *   matrix[elements-1], and all other elements to lie
 *   between these two.
 * destination: MPI rank of receing process.
 */
void smpi_send_matrix(double *matrix, size_t elements, int destination, int tag) {
	int err, chunksize = INT_MAX, tosend;
	size_t sent = 0;

	while (sent < elements) {
		if (elements-sent > (size_t)chunksize) tosend = chunksize;
		else tosend = elements-sent;

		if ((err=MPI_Send(matrix+sent, tosend, MPI_DOUBLE, destination, tag, MPI_COMM_WORLD))!=MPI_SUCCESS) {
			fprintf(stderr, "MPI ERROR %d: Failed to send matrix over MPI to process %d. Aborting...\n", err, destination);
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		sent += tosend;
	}
}

