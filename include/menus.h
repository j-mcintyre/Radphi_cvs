/* menus.h 4.1 
   10 Dec 1995 18:48:34
   /auto/home/lemond/e852/source/include/SCCS/s.menus.h
   @(#)menus.h	4.1
*/
typedef struct _menu_item{
  char *label;
  Widget w;
  WidgetClass *class;
  char mnemonic;
  char *accelerator;
  char *accel_text;
  void (*callback)();
  XtPointer callback_data;
  struct _menu_item *subitems;
} MenuItem;

Widget BuildMenu(Widget,int,char *,char,MenuItem*);
