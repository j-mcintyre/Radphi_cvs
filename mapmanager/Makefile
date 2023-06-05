#! gnumake

OBJECTS = map_create.o map_add_item.o map_add_item1.o map_rem_item.o \
          map_rem_sub.o \
          map_next_sub.o map_next_item.o map_put_int.o map_put_float.o \
          map_put_char.o \
          map_get_int.o map_get_float.o map_get_char.o map_fnd_int.o \
          map_fnd_float.o \
          map_fillname.o map_find_pred.o \
          map_error.o map_syserror.o map_warn.o map_log_mess.o \
	  map_read.o map_write.o map_add_sub.o \
          map_overwrite.o map_find_same.o map_time_pred.o map_dump.o \
          map_find_item.o map_stat_item.o map_open_close.o \
          map_find_sub.o map_fix_item.o map_fix_sub.o map_fix_head.o \
	  map_stat_head.o map_rem_arr.o \
          map_add_item_.o map_put_int_.o map_put_flo_.o map_get_int_.o \
          map_get_flo_.o map_fnd_int_.o map_fnd_flo_.o map_stat_item_.o \
          map_rem_item_.o map_rem_sub_.o map_rem_arr_.o map_dump_.o \
          map_log_mess_.o \
	  mapf_add_item.o mapf_get_float.o mapf_get_int.o mapf_put_float.o \
	  mapf_put_int.o mapf_rem_item.o mapf_stat_item.o \
	  mapf_get_char.o mapf_put_char.o \
	  map_read_int.o map_read_float.o map_read_char.o \
	  map_read_pointer.o map_read_name.o map_read_table.o \
	  map_read_nametable.o map_read_subsystem.o map_read_arrayheader.o \
	  map_read_subtableheader.o map_read_item.o \
	  map_write_int.o map_write_float.o map_write_char.o \
	  map_write_pointer.o map_write_name.o map_write_table.o \
	  map_write_nametable.o map_write_subsystem.o map_write_arrayheader.o \
	  map_write_subtableheader.o map_write_item.o \
	  map_overwrite_int.o map_overwrite_float.o map_overwrite_char.o \
	  map_overwrite_pointer.o map_overwrite_name.o map_overwrite_table.o \
	  map_overwrite_nametable.o map_overwrite_subsystem.o \
	  map_overwrite_arrayheader.o \
	  map_overwrite_subtableheader.o map_overwrite_item.o

#
#	read only routines
#
OBJECTS_RO = map_next_sub.o map_next_item.o \
          map_get_int.o map_get_float.o map_get_char.o map_fnd_int.o \
          map_fnd_float.o \
          map_fillname.o map_find_pred.o \
          map_error.o map_syserror.o map_warn.o map_log_mess.o \
	  map_read.o \
          map_find_same.o map_time_pred.o map_dump.o \
          map_find_item.o map_stat_item.o map_open_close_ro.o \
          map_find_sub.o \
	  map_stat_head.o \
          map_get_int_.o \
          map_get_flo_.o map_fnd_int_.o map_fnd_flo_.o map_stat_item_.o \
          map_dump_.o map_log_mess_.o \
	  mapf_get_float.o mapf_get_int.o \
	  mapf_stat_item.o map_readonly_stub.o \
	  mapf_get_char.o mapf_put_char.o \
	  map_read_int.o map_read_float.o map_read_char.o \
	  map_read_pointer.o map_read_name.o map_read_table.o \
	  map_read_nametable.o map_read_subsystem.o map_read_arrayheader.o \
	  map_read_subtableheader.o map_read_item.o
#
#	No Fortran: read only routines
#
OBJECTS_NOF_RO = map_next_sub.o map_next_item.o \
          map_get_int.o map_get_float.o map_get_char.o map_fnd_int.o \
          map_fnd_float.o \
          map_fillname.o map_find_pred.o \
          map_error.o map_syserror.o map_warn.o map_log_mess.o \
	  map_read.o \
          map_find_same.o map_time_pred.o map_dump.o \
          map_find_item.o map_stat_item.o map_open_close_ro.o \
          map_find_sub.o \
	  map_stat_head.o \
          map_get_int_.o \
          map_get_flo_.o map_fnd_int_.o map_fnd_flo_.o map_stat_item_.o \
          map_dump_.o map_log_mess_.o  map_readonly_stub.o \
	  map_read_int.o map_read_float.o map_read_char.o \
	  map_read_pointer.o map_read_name.o map_read_table.o \
	  map_read_nametable.o map_read_subsystem.o map_read_arrayheader.o \
	  map_read_subtableheader.o map_read_item.o \


LIBFILE =  $(LIBDIR)/libmapmanager.a
LIBFILE_RO = $(LIBDIR)/libmapmanager_ro.a
LIBFILE_VX_RO = $(LIBDIR)/libmapmanager_vx_ro.a

#
# target definitions:
#

include $(RADPHI_HOME)/makefiles/rules.make


lib: $(LIBFILE) $(LIBFILE_RO)

$(LIBFILE_VX_RO):  $(OBJECTS_NOF_RO)
	ldppc -o vxmman.o -r $?	

#	$(LD) -o $@ $(OBJECTS_RO)

$(LIBFILE): $(OBJECTS)
	$(AR) rv $@ $?	

$(LIBFILE_RO): $(OBJECTS_RO)
	$(AR) rv $@ $?	


#
#  make production files
#

exe: mapcreate mapdump mapmerge maptoascii mapfrascii mapfixit mapPutItem map_create_float_item map_create_int_item map_remove_array map_create_char_item map_fill_char_item_with_zeros
exe-install: mapcreate.install mapdump.install mapmerge.install maptoascii.install mapfrascii.install mapfixit.install mapPutItem.install map_create_float_item.install map_create_int_item.install map_remove_array.install map_create_char_item.install map_fill_char_item_with_zeros.install

mapcreate : mapcreate.o -lmapmanager -ldata
	$(CC) $(ARCHFLAGS)  -o $@ $^

mapdump : mapdump.o -lmapmanager -ldata
	$(CC) $(ARCHFLAGS)  -o $@ $^

mapmerge : mapmerge.o -lmapmanager -ldata
	$(CC) $(ARCHFLAGS)  -o $@ $^

maptoascii : maptoascii.o -lmapmanager -ldata
	$(CC) $(ARCHFLAGS)  -o $@ $^

mapfrascii : mapfrascii.o -lmapmanager -ldata
	$(CC) $(ARCHFLAGS)  -o $@ $^

mapPutItem: mapPutItem.o -lmapmanager -ldata
	$(CC) $(ARCHFLAGS) -o $@ $^

mapfixit : mapfixit.o -lmapmanager -ldata
	$(CC) $(ARCHFLAGS)  -o $@ $^

map_create_subsystem: map_create_subsystem.o -lmapmanager -ldata
	$(CC) $(ARCHFLAGS)  -o $@ $^

map_create_float_item: map_create_float_item.o -lmapmanager -ldata
	$(CC) $(ARCHFLAGS)  -o $@ $^

map_create_int_item: map_create_int_item.o -lmapmanager -ldata
	$(CC) $(ARCHFLAGS)  -o $@ $^

map_create_char_item: map_create_char_item.o -lmapmanager -ldata
	$(CC) $(ARCHFLAGS)  -o $@ $^

map_remove_array: map_remove_array.o  -lmapmanager -ldata
	$(CC) $(ARCHFLAGS)  -o $@ $^

map_fill_char_item_with_zeros: map_fill_char_item_with_zeros.o -lmapmanager -ldata
	$(CC) $(ARCHFLAGS)  -o $@ $^

#
#  make test programs
#

test : test_sig test_map2 test1 test2 test_map_c test_map_f test_big_write \
       test_big_read test_read test_write test2_read test2_write

#
#  make test program and library in local directory
#

test2 : test2.o libmapmanager.a -ldata
	cc -o $@ $^

test_sig : test_sig.o libmapmanager.a -ldata
	cc -o $@ $^
#
#  make test program and library in local directory
#

test1 : test1.o libmapmanager.a -ldata
	cc -o $@ $^

#
#  Make test program and library in local directory
#

test_map_c : test_map_c.o libmapmanager.a -ldata
	cc -o $@ $^
test_map2 : test_map2.o libmapmanager.a -ldata
	cc -o $@ $^
test_big_write : test_big_write.o libmapmanager.a -ldata
	cc -o $@ $^
test_big_read : test_big_read.o libmapmanager.a -ldata
	cc -o $@ $^
test_map_f : test_map_f.o libmapmanager.a -ldata
	f77 -o $@ $^
test_write : test_write.o libmapmanager.a -ldata
	cc -o $@ $^
test_read : test_read.o libmapmanager.a -ldata
	cc -o $@ $^
test2_write : test2_write.o libmapmanager.a -ldata
	cc -o $@ $^
test2_read : test2_read.o libmapmanager.a -ldata
	cc -o $@ $^
#
#  make program in local directory
#
#
#	compile open and close routines without locking.
#

map_open_close_ro.o : map_open_close.c
	$(CC) $(CDBG) $(ARCHFLAGS) $(CFLAGS) map_open_close.c -c -o map_open_close_ro.o $(INCLUDE) -Dreadonly


-include *.depend
