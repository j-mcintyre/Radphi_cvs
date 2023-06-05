/* 
        mapfrascii.c
         Created           : 11-May-1993 by Thom Sulanke
*/

static char rcsid[] = "$Id: mapfrascii.c,v 1.2 1996/11/08 15:03:58 lfcrob Exp $";

/* $Log: mapfrascii.c,v $
/* Revision 1.2  1996/11/08 15:03:58  lfcrob
/* Added $Log$
/* */ 

#include <stdio.h>                   /* I/O definitions                      */
#include <string.h>
#include <map_manager.h>
#include <map_internal.h>
#include <stdlib.h>

int read_subsystem( char subsystemname[]);
int read_item( char itemname[], int *length, int *type);
int read_int_array( int length, int *atime, int array[]);
int read_float_array( int length, int *atime, float farray[]);
int read_char_array( int *length, int *atime, char carray[]);
char *progname;
 

/***********************************************************************/
/*                                                                     */
/*   mapfrascii                                                        */
/*   ----------                                                        */
/*                                                                     */
/*         Created     : 11-May-1993    Author : Thom Sulanke          */
/*         Purpose     : read information in ascii for transfer from   */
/*                       incompatible systems                          */
/*                                                                     */
/***********************************************************************/
 
main( argc, argv )

     int argc;
     char *argv[];
{
  char subsystemname[MAP_NAMELENGTH], itemname[MAP_NAMELENGTH];
  char nextsubsystemname[MAP_NAMELENGTH], nextitemname[MAP_NAMELENGTH];
  int i;                /* index */
  int *array;
  char *carray;
  float *farray;
  int atime, vtime;
  int length, olength;
  int type;
  int nvals;

  progname = argv[0];

  if ( argc != 2 )
    {
      fprintf(stderr,"Usage: %s themap\n", progname);
      exit(0);
    }

  while( read_subsystem(subsystemname) )
    {
      while( read_item(itemname,&length,&type) )
	{
	  map_add_item( argv[1], subsystemname, itemname, length, type );
	  switch ( type ) {
	  case 0:
	    array = (int *) calloc(length,sizeof(int));
	    if ( array == NULL )
	      {
		fprintf(stderr,"%s: error allocating memory for array.\n", progname);
		exit(0);
	      }
	    while ( read_int_array(length, &atime, array) )
	      map_put_int( argv[1], subsystemname, itemname, length,
			  array, atime);
	    free(array);      
	    break;
	  case 1:
	    farray = (float *) calloc(length,sizeof(float));
	    if ( farray == NULL )
	      {
		fprintf(stderr,"%s: error allocating memory for array.\n", progname);
		exit(0);
	      }
	    while ( read_float_array(length, &atime, farray) )
	      map_put_float( argv[1], subsystemname, itemname, length,
			    farray, atime);
	    free(farray);      
	    break;
	  case 2:
	    carray = (char *) calloc(length,sizeof(char));
	    if ( carray == NULL )
	      {
		fprintf(stderr,"%s: error allocating memory for array.\n", progname);
		exit(0);
	      }
	    while ( read_char_array(&length, &atime, carray) )
	      map_put_char( argv[1], subsystemname, itemname, length,
			    carray, atime);
	    free(carray);      
	    break;
	  } /* switch ( type ) */
	} /* while */
    } /* while */
  
  return 0;
}

int read_subsystem( char subsystemname[])

{
  char c[2];
  int i;
  int int_tmp;

  if ( scanf("%1s",c) != 1 )
    {
      fprintf(stderr,"error when expecting S.\n");
      exit(0);
    }

  if ( c[0] == 'E' )
    return 0;

  if ( c[0] != 'S' )
    {
      fprintf(stderr,"saw %s when expecting S.\n",c);
      exit(0);
    }

  for(i=0;i<MAP_NAMELENGTH;i++)
  {
  	scanf("%x",&int_tmp);
  	subsystemname[i] = int_tmp;
  }
  return 1;
}

int read_item( char itemname[], int *length, int *type)

{
  char c[2];
  int i;
  int int_tmp;

  if ( scanf("%1s",c) != 1 )
    {
      fprintf(stderr,"error when expecting I.\n");
      exit(0);
    }

  if ( c[0] == 'E' )
    return 0;

  if ( c[0] != 'I' )
    {
      fprintf(stderr,"saw %s when expecting I.\n",c);
      exit(0);
    }

  scanf("%d %d",length,type);
  for(i=0;i<MAP_NAMELENGTH;i++)
  {
  	scanf("%x",&int_tmp);
  	itemname[i] = int_tmp;
  }
  return 1;
}

int read_int_array( int length, int *atime, int array[])

{
  char c[2];
  int i;

  if ( scanf("%1s",c) != 1 )
    {
      fprintf(stderr,"error when expecting A.\n");
      exit(0);
    }

  if ( c[0] == 'E' )
    return 0;

  if ( c[0] != 'A' )
    {
      fprintf(stderr,"saw %s when expecting A.\n",c);
      exit(0);
    }

  scanf("%d",atime);
  for ( i=0; i<length; i++ )
  scanf("%d",&array[i]);
  return 1;
}

int read_float_array( int length, int *atime, float farray[])

{
  char c[2];
  int i;

  if ( scanf("%1s",c) != 1 )
    {
      fprintf(stderr,"error when expecting F.\n");
      exit(0);
    }

  if ( c[0] == 'E' )
    return 0;

  if ( c[0] != 'F' )
    {
      fprintf(stderr,"saw %s when expecting F.\n",c);
      exit(0);
    }

  scanf("%d",atime);
  for ( i=0; i<length; i++ )
    {
      char s[100];

      /*scanf("%f",&farray[i]);   <--- this cannot handle NaNs. */

      scanf("%s",s); farray[i] = strtod(s,NULL);

      /*printf("Value: %lf\n",farray[i]);*/
    }

  return 1;
}

int read_char_array( int *length, int *atime, char carray[])

{
  char c[2];
  int i;
  int int_tmp;

  if ( scanf("%1s",c) != 1 )
    {
      fprintf(stderr,"error when expecting C.\n");
      exit(0);
    }

  if ( c[0] == 'E' )
    return 0;

  if ( c[0] != 'C' )
    {
      fprintf(stderr,"saw %s when expecting C.\n",c);
      exit(0);
    }

  scanf("%d %d",atime,length);
  for ( i=0; i<*length; i++ )
  {
  	scanf("%x",&int_tmp);
  	carray[i] = int_tmp;
  	}
  return 1;
}

