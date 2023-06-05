/* public function members (c interface) */

int fileOpen(char *fname, char *pmode);
int fileGetEvent(int fd, itape_header_t *event);
int fileSkipEvents(int fd, int count);
void fileClose(int fd);

/* public function members (fortran interface) */
int filetest_(char *infile);
