#ifndef _SMPI_H
#define _SMPI_H

#define SMPI_OUTPUT_READY 255

void smpi_sor(int);
void smpi_wor(int);
void smpi_receive_matrix(double*, size_t, int, int);
void smpi_send_matrix(double*, size_t, int, int);

#endif/*_SMPI_H*/
