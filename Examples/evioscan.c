/*
 * evioscan - stream through an evio data file (vintage 2000-2001)
 *            and look for unusual or unexpected banks.
 *
 * author: richard.t.jones at uconn.edu
 * version: april 4, 2013
 */

#include <stdio.h>
#include <stdlib.h>

void swap_endian(int *buffer, int length);

int main(int argc, char *argv[]) {
   FILE *evtfile;
   int block_size;
   int length;
   int max_length = 9999;
   int event[max_length];
   int total_words_read = 0;

   // open the input file listed on the command line
   if (argc != 2) {
      printf("Usage: evioscan <datafile.evt>\n");
      exit(1);
   }
   evtfile = fopen(argv[1],"r");
   if (evtfile == 0) {
      printf("Error: unable to open input file %s\n",argv[1]);
      exit(2);
   }

   fread(event,sizeof(int),8,evtfile);
   total_words_read += 8;
   swap_endian(event,1);
   block_size = event[0];

   while (fread(&length,sizeof(int),1,evtfile) > 0) {
      total_words_read += 1;

      swap_endian(&length,1);

      if (total_words_read > block_size) {
         fread(event,sizeof(int),7,evtfile);
         fread(&length,sizeof(int),1,evtfile);
         swap_endian(&length,1);
         total_words_read = 9;
      }
      if (length > max_length) {
         printf("Error: input event overflows buffer: %d > %d\n",
                length,max_length);
         exit(3);
      }
      else if (fread(event,sizeof(int),length,evtfile) != length) {
         printf("Error: read error on input file: requested %d\n",length);
            {
               fread(event,4,20,evtfile);
               printf("0x%x 0x%x 0x%x 0x%x 0x%x 0x%x 0x%x 0x%x 0x%x 0x%x\n",
                      event[0],event[1],event[2],
                      event[2],event[3],event[4],
                      event[5],event[6],event[7],event[8]);
               printf("0x%x 0x%x 0x%x 0x%x 0x%x 0x%x 0x%x 0x%x 0x%x 0x%x\n",
                      event[9],event[10],event[11],
                      event[12],event[13],event[14],
                      event[15],event[16],event[17],event[18]);
            }
         printf("total words read was %d\n",total_words_read);
         exit(4);
      }
      else if (total_words_read + length > block_size) {
         int i;
         fread(&event[length],sizeof(int),8,evtfile);
         for (i = block_size - total_words_read; i < length; ++i) {
            event[i] = event[i+8];
         }
         total_words_read += 8 - block_size;
      }
      total_words_read += length;

      swap_endian(event,length);

      if (length > 0) {
         int event_tag = (*event >> 16) & 0xffff;
         int event_dtype = (*event >> 8) & 0xff;
         int event_num = *event & 0xff;
         printf("event size %d tag 0x%x dtype %d num 0x%x\n",
                length,event_tag,event_dtype,event_num);
         if (event_tag < 0x11) {
            int *bankptr;
            int *bankptr_begin;
            int *bankptr_end;
            bankptr_begin = event + 1;
            bankptr_end = event + length;
            for (bankptr = bankptr_begin; bankptr < bankptr_end;) {
               int banksize = *(bankptr++);
               int bank_tag = (*bankptr >> 16) & 0xffff;
               int bank_dtype = (*bankptr >> 8) & 0xff;
               int bank_num = *bankptr & 0xff;
               if (bank_tag == 0xc000) {
                  printf("   * event number %d\n",bankptr[1]);
               }
               printf("   - tag 0x%x dtype %d num 0x%x size %d\n",
                      bank_tag,bank_dtype,bank_num,banksize);
               bankptr += banksize;
            }
            if (bankptr == bankptr_end) {
               continue;
            }
            else {
               printf("Error: event does not terminate at end of bank!\n");
               exit(5);
            }
         }
      }
   }
   fclose(evtfile);
   exit(0);
}
      
void swap_endian(int *buffer, int length) {
   int n;
   int oneword;
   char *dptr = (char*)&oneword;
   for (n = 0; n < length; ++n) {
      char *sptr = (char*)&buffer[n];
      dptr[3] = sptr[0];
      dptr[2] = sptr[1];
      dptr[1] = sptr[2];
      dptr[0] = sptr[3];
      buffer[n] = oneword;
   }
}
