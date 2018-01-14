#cmake_minimum_required(VERSION 2.8)

# Basic compiler/linker configuration
#file(GLOB SOURCES "*.c")
#file(GLOB SOURCES "equations/*.c")
#file(GLOB SOURCES "magnetic/*.c")
#file(GLOB SOURCES "tools/*.c")
#file(GLOB SOURCES "tools/sycamera/*.c")
#add_executable(soft ${SOURCES})

set(CMAKE_MODULE_PATH "${CMAKE_BINARY_DIR}/Modules/" CACHE PATH "CMake module path" FORCE)

include_directories("${PROJECT_SOURCE_DIR}/src/include")

# Set libraries to use
option(BUILD_TESTS "Build the test programs" OFF)
option(DEBUG "Compile with debug symbols and no optimzations" OFF)
option(USE_HDF5 "Link in support for HDF5" OFF)
option(USE_MATLAB "Link in support for MATLAB *.mat files" ON)
option(USE_MPI "Compile with MPI support" OFF)

set(mainprogram
	"${PROJECT_SOURCE_DIR}/src/main.c"
)
set(main
	"${PROJECT_SOURCE_DIR}/src/counter.c"
	"${PROJECT_SOURCE_DIR}/src/ctsv.c"
	"${PROJECT_SOURCE_DIR}/src/diag.c"
	"${PROJECT_SOURCE_DIR}/src/domain.c"
	"${PROJECT_SOURCE_DIR}/src/particles.c"
	"${PROJECT_SOURCE_DIR}/src/quantities.c"
	"${PROJECT_SOURCE_DIR}/src/readfile.c"
	"${PROJECT_SOURCE_DIR}/src/rkf45.c"
	"${PROJECT_SOURCE_DIR}/src/settings.c"
	"${PROJECT_SOURCE_DIR}/src/sfile.c"
	"${PROJECT_SOURCE_DIR}/src/ssdt.c"
	"${PROJECT_SOURCE_DIR}/src/util.c"
	"${PROJECT_SOURCE_DIR}/src/vector.c"
)

set(equations
	"${PROJECT_SOURCE_DIR}/src/equations/equations.c"
	"${PROJECT_SOURCE_DIR}/src/equations/GCM.c"
	"${PROJECT_SOURCE_DIR}/src/equations/GCM_rel.c"
	"${PROJECT_SOURCE_DIR}/src/equations/particle.c"
	"${PROJECT_SOURCE_DIR}/src/equations/particle_rel.c"
)

set(magnetic
	"${PROJECT_SOURCE_DIR}/src/magnetic/axis.c"
	"${PROJECT_SOURCE_DIR}/src/magnetic/interp2.c"
	"${PROJECT_SOURCE_DIR}/src/magnetic/magfield.c"
	"${PROJECT_SOURCE_DIR}/src/magnetic/magnetic_circ.c"
	"${PROJECT_SOURCE_DIR}/src/magnetic/magnetic_field.c"
	"${PROJECT_SOURCE_DIR}/src/magnetic/magnetic_num.c"
)

set(sycamera
	#"${PROJECT_SOURCE_DIR}/src/tools/sycamera/ascii.c"
	#"${PROJECT_SOURCE_DIR}/src/tools/sycamera/binary.c"
	"${PROJECT_SOURCE_DIR}/src/tools/sycamera/bsspectrum.c"
	"${PROJECT_SOURCE_DIR}/src/tools/sycamera/bremsdist.c"
	"${PROJECT_SOURCE_DIR}/src/tools/sycamera/cone_delta.c"
	"${PROJECT_SOURCE_DIR}/src/tools/sycamera/cone_dist.c"
	"${PROJECT_SOURCE_DIR}/src/tools/sycamera/df_interp.c"
	"${PROJECT_SOURCE_DIR}/src/tools/sycamera/df_readfile.c"
	"${PROJECT_SOURCE_DIR}/src/tools/sycamera/distfunc.c"
	"${PROJECT_SOURCE_DIR}/src/tools/sycamera/green.c"
	#"${PROJECT_SOURCE_DIR}/src/tools/sycamera/hdf5_sfile.c"
	"${PROJECT_SOURCE_DIR}/src/tools/sycamera/image.c"
	"${PROJECT_SOURCE_DIR}/src/tools/sycamera/isotropic.c"
	#"${PROJECT_SOURCE_DIR}/src/tools/sycamera/map.c"
	"${PROJECT_SOURCE_DIR}/src/tools/sycamera/pas2.c"
	"${PROJECT_SOURCE_DIR}/src/tools/sycamera/pcyl.c"
	"${PROJECT_SOURCE_DIR}/src/tools/sycamera/pcyl_lookup.c"
	"${PROJECT_SOURCE_DIR}/src/tools/sycamera/pcyl_lookup_pol.c"
	"${PROJECT_SOURCE_DIR}/src/tools/sycamera/pdist.c"
	"${PROJECT_SOURCE_DIR}/src/tools/sycamera/pdist_lookup.c"
	"${PROJECT_SOURCE_DIR}/src/tools/sycamera/pdist_spec_lookup.c"
	"${PROJECT_SOURCE_DIR}/src/tools/sycamera/polimage.c"
	"${PROJECT_SOURCE_DIR}/src/tools/sycamera/polspectrometer.c"
	"${PROJECT_SOURCE_DIR}/src/tools/sycamera/ptotal.c"
	"${PROJECT_SOURCE_DIR}/src/tools/sycamera/space3d.c"
	"${PROJECT_SOURCE_DIR}/src/tools/sycamera/spectrometer.c"
	"${PROJECT_SOURCE_DIR}/src/tools/sycamera/sphere.c"
	"${PROJECT_SOURCE_DIR}/src/tools/sycamera/sycamera.c"
	"${PROJECT_SOURCE_DIR}/src/tools/sycamera/sycout.c"
	"${PROJECT_SOURCE_DIR}/src/tools/sycamera/topview.c"
	"${PROJECT_SOURCE_DIR}/src/tools/sycamera/lookup/Zeta11_lookup.c"
	"${PROJECT_SOURCE_DIR}/src/tools/sycamera/lookup/Zeta12_lookup.c"
	"${PROJECT_SOURCE_DIR}/src/tools/sycamera/lookup/Zeta21_lookup.c"
)

set(tests
	"${PROJECT_SOURCE_DIR}/tests/main.c"
	"${PROJECT_SOURCE_DIR}/tests/bremsspec.c"
	"${PROJECT_SOURCE_DIR}/tests/bremsdist.c"
	"${PROJECT_SOURCE_DIR}/tests/distfunc.c"
	"${PROJECT_SOURCE_DIR}/tests/hyperbola.c"
)
set(tools
	"${PROJECT_SOURCE_DIR}/src/tools/orbit.c"
	"${PROJECT_SOURCE_DIR}/src/tools/tools.c"
)

if (DEBUG)
	message(STATUS "Compiling in DEBUG mode")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -Wno-unused-parameter -fopenmp -O0 -g -pg -no-pie -D_FILE_OFFSET_BITS=64")
else (DEBUG)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -Wno-unused-parameter -fopenmp -O3 -march=native -D_FILE_OFFSET_BITS=64")
endif (DEBUG)

if (USE_HDF5)
	set(main ${main} "shdf5.c")
endif (USE_HDF5)
if (USE_MATLAB)
	set(main ${main} "smat.c")
endif (USE_MATLAB)
if (USE_MPI)
	set(main ${main} "smpi.c")
endif (USE_MPI)

add_executable(soft ${mainprogram} ${main} ${equations} ${magnetic} ${tools} ${sycamera})
target_link_libraries(soft m)
install(TARGETS soft DESTINATION /usr/bin)

# Test programs
if (BUILD_TESTS)
	add_executable(runtest ${tests} ${main} ${equations} ${magnetic} ${tools} ${sycamera})
	target_link_libraries(runtest m)
endif (BUILD_TESTS)

# Find GSL!
find_package(GSL REQUIRED)
if (GSL_FOUND)
	if ((GSL_VERSION VERSION_GREATER "2.0") OR (GSL_VERSION VERSION_EQUAL "2.0"))
		include_directories(${GSL_INCLUDE_DIRS})
		target_link_libraries(soft ${GSL_LIBRARIES})
		target_link_libraries(soft ${GSL_CBLAS_LIBRARY})

		if (BUILD_TESTS)
			target_link_libraries(runtest ${GSL_LIBRARIES})
			target_link_libraries(runtest ${GSL_CBLAS_LIBRARY})
		endif (BUILD_TESTS)
	else ((GSL_VERSION VERSION_GREATER "2.0") OR (GSL_VERSION VERSION_EQUAL "2.0"))
		message(FATAL_ERROR "SOFT needs GSL >= 2.0")
	endif ((GSL_VERSION VERSION_GREATER "2.0") OR (GSL_VERSION VERSION_EQUAL "2.0"))
endif (GSL_FOUND)

# Compile with HDF5 support?
if (USE_HDF5)
	find_package(HDF5 REQUIRED)
	if (HDF5_FOUND)
		include_directories(${HDF5_INCLUDE_DIRS})
		message(STATUS "Path to HDF5: ${HDF5_INCLUDE_DIRS}")
		target_link_libraries(soft ${HDF5_LIBRARIES})

		if (BUILD_TESTS)
			target_link_libraries(runtest ${HDF5_LIBRARIES})
		endif (BUILD_TESTS)
	else (HDF5_FOUND)
		message(FATAL_ERROR "HDF5 support was requested, but no HDF5 library was found on this system")
	endif (HDF5_FOUND)
endif (USE_HDF5)

# Compile with MATLAB support?
if (USE_MATLAB)
	find_package(Matlab COMPONENTS MAT_LIBRARY MX_LIBRARY)
	if (Matlab_FOUND)
		include_directories(${Matlab_INCLUDE_DIRS})
		message(${Matlab_MAT_LIBRARY})
		message(${Matlab_MX_LIBRARY})
		target_link_libraries(soft ${Matlab_MAT_LIBRARY} ${Matlab_MX_LIBRARY})

		if (BUILD_TESTS)
			target_link_libraries(runtest ${Matlab_MAT_LIBRARY} ${Matlab_MX_LIBRARY})
		endif (BUILD_TESTS)
	else (Matlab_FOUND)
		message(FATAL_ERROR "MATLAB support was requested, but no MATLAB installation was found")
	endif (Matlab_FOUND)
endif (USE_MATLAB)

# Find OpenMP!
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OPENMP_C_FLAGS}")
endif (OPENMP_FOUND)

# Find MPI
if (USE_MPI)
	find_package(MPI)
	if (MPI_C_FOUND)
		include_directories(${MPI_C_INCLUDE_PATH})

		set_target_properties(soft PROPERTIES COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS}")
		set_target_properties(soft PROPERTIES LINK_FLAGS "${MPI_C_LINK_FLAGS}")
		target_link_libraries(soft "${MPI_C_LIBRARIES}")

		if (BUILD_TESTS)
			set_target_properties(runtest PROPERTIES COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS}")
			set_target_properties(runtest PROPERTIES LINK_FLAGS "${MPI_C_LINK_FLAGS}")
			target_link_libraries(runtest "${MPI_C_LIBRARIES}")
		endif (BUILD_TESTS)
	else (MPI_C_FOUND)
		message(FATAL_ERROR "MPI support was requested, but no MPI library was found")
	endif (MPI_C_FOUND)
endif (USE_MPI)

configure_file(
	"${PROJECT_SOURCE_DIR}/src/include/config.h.in"
	"${PROJECT_SOURCE_DIR}/src/include/config.h"
)

