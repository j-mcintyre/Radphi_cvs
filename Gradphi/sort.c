/* Following is the subroutine that sorts array hits and returns it 
 * in the array sorthits
*/
void sort(float *hits,float *sorthits,int dimension,int nhits) {

 float min;
 int n;
 int k;
 int l;
 
 min=100; /* some ridiculously huge time */
 if (nhits>0) {
  for (n=0;n<nhits;n++) { /* n is the index of the SORTED array */
   for(k=0;k<nhits;k++) { /* looking for the minimum */ 
    if(min>hits[6+dimension*k]) {
     min=hits[6+dimension*k];
     l=k;
    }
   } /* now min contains the shortest time, and l=number of this hit */
   for (k=0;k<dimension;k++) { /* fill the sorted hits array */
    sorthits[k+dimension*n]=hits[k+dimension*l];
   }        
   hits[6+dimension*l]=622; /* to ensure that this entry will not be minimum anymore */
   min=100; /* reset the minimum value so that next search can be done */
  }
 } 
}   

