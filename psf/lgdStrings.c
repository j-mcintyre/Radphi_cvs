/* * * * * * * * * * * * * * * *
* lgdStrings - September 2, 2002
*
* Searches for the longest string of events in each run
* that contains a given channel hit or not hit.  The results
* are printed at the end of each run in 3-column ascii format
* into a file named according to the format string fmt[].
*/

#define HIT_THRESHOLD 0.010

void lgdStrings(char *fmt) {
  
   int maxStringOn[LGD_MATRIX_ROWS*LGD_MATRIX_COLS];
   int maxStringOff[LGD_MATRIX_ROWS*LGD_MATRIX_COLS];
   int stringOn[LGD_MATRIX_ROWS*LGD_MATRIX_COLS];
   int stringOff[LGD_MATRIX_ROWS*LGD_MATRIX_COLS];
   int totalOn[LGD_MATRIX_ROWS*LGD_MATRIX_COLS];

   int run=0;
   int chan;
 
   while (1 > 0) {
      int ev = eventGetNext();
      if ((run != runNumber) || (ev < 0)) {
         if (run) {
            FILE *runFile;
            char fileName[300];
            if ((sprintf(fileName,fmt,run)) &&
                (runFile=fopen(fileName,"w+"))) {
               fprintf(runFile,"EventNo %d\n",event->eventNo);
               for (chan=0;chan<ndim;chan++) {
                  if (stringOff[chan] > maxStringOff[chan]) {
	             maxStringOff[chan] = stringOff[chan];
                  }
                  if (stringOn[chan] > maxStringOn[chan]) {
	             maxStringOn[chan] = stringOn[chan];
                  }
                  fprintf(runFile,"%d\t%d\t%d\t%d\n",chan,maxStringOn[chan],
                                                          maxStringOff[chan],
                                                          totalOn[chan]);
               }
               fclose(runFile);
            }
            else {
               fprintf(stderr,"ERROR - cannot open output file %s\n",fileName);
            }
         }
         memset(maxStringOn,0,sizeof(int)*ndim);
         memset(maxStringOff,0,sizeof(int)*ndim);
         memset(stringOn,0,sizeof(int)*ndim);
         memset(stringOff,0,sizeof(int)*ndim);
         memset(totalOn,0,sizeof(int)*ndim);
         run = runNumber;
      }

      if (ev < 0) {
         break;
      }
      else if ((lgdHits==NULL) || (lgdHits->nhits==0)) {
         continue;   
      }

      for (chan=0;chan<ndim;chan++) {
         int row=chan/LGD_MATRIX_ROWS;
         int col=chan%LGD_MATRIX_ROWS;
         if (lgd_image.yield[row][col] > HIT_THRESHOLD) {
            if (stringOff[chan] > maxStringOff[chan]) {
	       maxStringOff[chan] = stringOff[chan];
            }
            stringOff[chan] = 0;
            stringOn[chan]++;
            totalOn[chan]++;
         }
         else {
            if (stringOn[chan] > maxStringOn[chan]) {
	       maxStringOn[chan] = stringOn[chan];
            }
            stringOn[chan] = 0;
            stringOff[chan]++;
         }
      }
   }
}

void lgdstrings_(char *fmt, int len) {
   char *form = malloc(len+1);
   strncpy(form,fmt,len);
   form[len] = ' ';
   *(strchr(form,' ')) = 0;
   lgdStrings(form);
   free(form);
}
