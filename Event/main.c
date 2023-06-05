/* The event display here represented was designed, coded & documented by   */
/* Troels Petersen under the supervision of Elton Smith, with help from     */
/*   Dave  Heddle and James Muguira.					    */
/* Comments can be found around more complicated code leaving the basics    */
/* uncommented. Upon problems and/or questions, please contact Troels       */
/* Petersen at: <petersen@cebaf.gov> or <petersen@fys.ku.dk>		    */
/*                                                                          */
/*       modifications by R.T. Jones  July 2000                             */
/*         - split the code in setup.c into setup.c, draw.c and input.c     */
/*           to help partition the evolution of this code.                  */
/*         - changed Hv_COLORWIDTH because the rainbow has more colors now  */
/*           than it did before and the color key overflowed the window.    */
/*         - increased the width of the feedback area to 50 columns.        */
/*         - changed color indexing from absolute on scale [0,28] to use    */
/*           symbolic color names of Hv_red...Hv_white.                     */
/*                                                                          */
/*       June 2000   DSA                                                    */
/*          - modified to read .evt files                                   */
/*          - modified to use makeHits and setup_makeHits from library      */
/*          - still to do: use "hits" groups for other detectors?           */
/*       Other problems: - Postscript hardcopy is wonky, especially for     */
/*                         Landscape-style plots                            */
/*                       - resizing?                                        */
/*                       - lots of hard-wired geometry dependence           */
/*                                                                          */
/*       Feb 2000    DSA                                                    */
/*          - modify to use HITS groups for Monte Carlo data                */
/*               no times yet, just ADC                                     */
/*               BGV needs to be fixed up                                   */
/*          - added Monte Carlo Event dump feature (should make an option)  */
/*                                                                          */
/*       Dec. 1999   DSA and Adam Gurson                                    */
/*          - added unpacking routines calls if RAW data, so can analyze    */
/*            .itp files as well as .itape files                            */
/*          - added makeLGDHits, makePhotonPatch, calls to clusterizer      */
/*          - fixed seg. fault for `Go to Event' when file not on pwd       */
/*                                                                          */
/*       June/July 1999      Lisa Kaufman and Dan Steiner                   */
/*          - added BSD, BGV, tdc and adc groups                            */
/*          - added BSD, BGV to main view                                   */
/*          - checks run number and runs RPD or BSD/BGV setup               */
/*          - added view for BSD                                            */
/*          - BSD, BGV data                                                 */
/*                                                                          */
/*       June/July 1998          D.S.A.                                     */
/*          - added MAM mass**2 word display                                */
/*          - initialize mam_energy, mam_mass properly                      */
/*          - creating a new view for the photon tagger                     */
/*          - keyboard accelerators added - altN for next event             */
/*          - CPV data  UPV data    (Sandy Sligh  and Jennifer Knowles)     */ 
/*                                                                          */
/*       April 1998          D.S.A.                                         */
/*          - added UPV, CPV tdc groups...                                  */
/*          - added new menu (Options) and pedestal on/off option           */
/*          - pedestal subtraction on/off.                                  */
/*          - add MAM energy word display, and threshold option             */
/*                                                                          */
/*       May 1997            D.S. Armstrong                                 */
/*         - new routine GetData                                            */
/*         - skip Reconstruction routine                                    */
/*         - don't try to count the number of events in the file anymore    */
/*         - remove, for now, previous event possibility                    */
/*         - modifications to ReadForward, ReadGoto, NewFile, MakeNew to    */
/*            use processed itape data                                      */
/*         - modifed Recoil to use ADC data for RPD hits                    */
/*         - new routine NewInput for data coming from stdin, not a file    */
/*         - CPV text output (crude) in LeadGlass                           */
/*         - add warnings if only one group found                           */
/*         - add simple scaler handling - just write to stdout              */
/*         - pedestal subtraction for Pb/glass                              */
/*         - trigger type written out                                       */
/*         - trigger types can be selected                                  */
/*         - added online (Dispatcher) capability                           */

#include "Event.h"

int main(unsigned int argc, char **argv)
{
  Hv_VaInitialize(argc, argv,
		  Hv_WELCOMETEXT,    "Version 3.0  R.T. Jones",
		  Hv_USEVIRTUALVIEW, True,
		  Hv_USEWELCOMEVIEW, True,
		  Hv_WIDTH,          1000,
		  Hv_HEIGHT,         850,
		  NULL);
  Init();
  Hv_Go();
}
