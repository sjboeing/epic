AM_FCFLAGS = 				\
	-I $(top_builddir)/src/3d/	\
	-I $(top_builddir)/src/utils	\
	-I $(top_builddir)/src/mpi	\
	-I $(top_builddir)/src/netcdf

AM_LDFLAGS = 						\
	$(top_builddir)/src/3d/utils/options.o		\
	$(top_builddir)/src/3d/utils/parameters.o


epic3d_modelsdir = $(prefix)/bin

epic3d_models_PROGRAMS = epic3d-models

epic3d_models_SOURCES =		\
	beltrami_3d.f90		\
	robert_3d.f90		\
	moist_3d.f90		\
	epic3d-models.f90

epic3d_models_LDADD = 					\
	$(top_builddir)/src/utils/libepic_utils.la	\
	$(top_builddir)/src/mpi/libepic_mpi.la		\
	$(top_builddir)/src/netcdf/libepic_netcdf.la

clean-local:
	rm -f *.mod
