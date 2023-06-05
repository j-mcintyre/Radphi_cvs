#! gnumake

include $(RADPHI_HOME)/makefiles/rules.make

LIBFILE = $(LIBDIR)/libHv.a

SRCS = \
          Hv_attributes.c \
          Hv_balloon.c \
          Hv_barplot.c \
          Hv_box.c \
          Hv_button.c \
          Hv_buttongroup.c \
          Hv_callback.c \
          Hv_choice.c \
          Hv_color.c \
          Hv_colorbutton.c \
          Hv_compound.c \
          Hv_configure.c \
          Hv_connect.c \
          Hv_container.c \
          Hv_dialogs_w32.c \
          Hv_dlogmenus.c \
          Hv_dlogs.c \
          Hv_dlogsupport.c \
          Hv_drag.c \
          Hv_draw.c \
          Hv_drawsupport.c \
          Hv_drawtools.c \
          Hv_feedback.c \
          Hv_fonts.c \
          Hv_fractdone.c \
          Hv_gateway.c \
          Hv_gateway_xm.c \
          Hv_gateway_w32.c \
          Hv_histo.c \
          Hv_init.c \
          Hv_io.c \
          Hv_items.c \
          Hv_led.c \
          Hv_linklist.c \
          Hv_mapitem.c \
          Hv_mapmath.c \
          Hv_mapreader.c \
          Hv_math.c \
          Hv_menubutton.c \
          Hv_menus.c \
          Hv_offscreen.c \
          Hv_opengl_w32.c \
          Hv_optionbutton.c \
          Hv_pic.c \
          Hv_plotsupport.c \
          Hv_pointer.c \
          Hv_postscript.c \
          Hv_print.c \
          Hv_privdraw.c \
          Hv_psdraw.c \
          Hv_radiodlogs.c \
          Hv_rainbowscale.c \
          Hv_rectsupport.c \
          Hv_region.c \
          Hv_scroll.c \
          Hv_separator.c \
          Hv_slider.c \
          Hv_standard.c \
          Hv_textentryitem.c \
          Hv_textitem.c \
          Hv_undo.c \
          Hv_utilities.c \
          Hv_views.c \
          Hv_virtual.c \
          Hv_wheel.c \
          Hv_worldimage.c \
          Hv_worldline.c \
          Hv_worldpolygon.c \
          Hv_worldrect.c \
          Hv_worldwedge.c \
          Hv_xyplot.c \
          Hv_xyplotedit.c \
          Hv_xyplotio.c

OBJS = ${SRCS:.c=.o}

lib: $(OBJS)
	$(AR) rv $(LIBFILE) $(OBJS)

exe:
	@echo "Hv only builds as a library."

