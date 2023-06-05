
#ifndef POPUPS_H_INCLUDED
#define POPUPS_H_INCLUDED
/* $Log: popups.h,v $
/* Revision 1.3  1997/01/30 18:25:43  lfcrob
/* Added *correct* log line
/*
 */
#include <Xm/XmAll.h>

#define WARNING_OK 1
#define WARNING_CANCEL 2

void showError(Widget w,char *string);
int showWarning(Widget w,char *string,XtAppContext app);

#endif
