# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/cheen/mc410

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/cheen/mc410

# Include any dependencies generated for this target.
include CMakeFiles/proj6B.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/proj6B.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/proj6B.dir/flags.make

CMakeFiles/proj6B.dir/proj6.cxx.o: CMakeFiles/proj6B.dir/flags.make
CMakeFiles/proj6B.dir/proj6.cxx.o: proj6.cxx
	$(CMAKE_COMMAND) -E cmake_progress_report /home/cheen/mc410/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/proj6B.dir/proj6.cxx.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/proj6B.dir/proj6.cxx.o -c /home/cheen/mc410/proj6.cxx

CMakeFiles/proj6B.dir/proj6.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/proj6B.dir/proj6.cxx.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/cheen/mc410/proj6.cxx > CMakeFiles/proj6B.dir/proj6.cxx.i

CMakeFiles/proj6B.dir/proj6.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/proj6B.dir/proj6.cxx.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/cheen/mc410/proj6.cxx -o CMakeFiles/proj6B.dir/proj6.cxx.s

CMakeFiles/proj6B.dir/proj6.cxx.o.requires:
.PHONY : CMakeFiles/proj6B.dir/proj6.cxx.o.requires

CMakeFiles/proj6B.dir/proj6.cxx.o.provides: CMakeFiles/proj6B.dir/proj6.cxx.o.requires
	$(MAKE) -f CMakeFiles/proj6B.dir/build.make CMakeFiles/proj6B.dir/proj6.cxx.o.provides.build
.PHONY : CMakeFiles/proj6B.dir/proj6.cxx.o.provides

CMakeFiles/proj6B.dir/proj6.cxx.o.provides.build: CMakeFiles/proj6B.dir/proj6.cxx.o

# Object files for target proj6B
proj6B_OBJECTS = \
"CMakeFiles/proj6B.dir/proj6.cxx.o"

# External object files for target proj6B
proj6B_EXTERNAL_OBJECTS =

proj6B: CMakeFiles/proj6B.dir/proj6.cxx.o
proj6B: CMakeFiles/proj6B.dir/build.make
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkFiltersFlowPaths-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkCommonExecutionModel-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkCommonDataModel-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkCommonMath-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkCommonCore-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtksys-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkCommonMisc-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkCommonSystem-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkCommonTransforms-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkFiltersGeneral-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkCommonComputationalGeometry-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkFiltersCore-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkFiltersSources-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkIOCore-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkzlib-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkverdict-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkFiltersHybrid-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkImagingSources-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkImagingCore-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkRenderingCore-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkCommonColor-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkFiltersExtraction-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkFiltersStatistics-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkImagingFourier-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkalglib-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkFiltersGeometry-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkIOAMR-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkFiltersAMR-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkParallelCore-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkIOLegacy-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkhdf5_hl-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkhdf5-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkglew-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkImagingStencil-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkRenderingContext2D-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkRenderingFreeType-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkfreetype-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkIOLSDyna-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkIOXML-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkIOGeometry-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkIOXMLParser-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkexpat-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtktiff-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkjpeg-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkInfovisCore-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkImagingColor-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkInfovisLayout-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkFiltersModeling-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkImagingHybrid-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkIOImage-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkDICOMParser-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkmetaio-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkpng-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkoggtheora-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkRenderingContextOpenGL2-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkRenderingOpenGL2-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkFiltersVerdict-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkjsoncpp-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtklibxml2-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkFiltersSelection-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkFiltersHyperTree-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkRenderingImage-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkDomainsChemistryOpenGL2-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkDomainsChemistry-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkViewsCore-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkInteractionWidgets-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkImagingGeneral-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkInteractionStyle-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkRenderingAnnotation-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkRenderingVolume-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkViewsContext2D-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkTestingIOSQL-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkIOSQL-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtksqlite-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkTestingGenericBridge-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkGeovisCore-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkproj4-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkIOExodus-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkexoIIc-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkNetCDF-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkNetCDF_cxx-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkRenderingLabel-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkIOMINC-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkIOExport-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkIOEnSight-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkFiltersProgrammable-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkFiltersTexture-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkIOMovie-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkFiltersGeneric-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkTestingRendering-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkIONetCDF-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkFiltersImaging-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkFiltersSMP-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkRenderingVolumeOpenGL2-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkIOParallel-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkFiltersParallel-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkImagingStatistics-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkImagingMorphological-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkImagingMath-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkFiltersParallelImaging-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkInteractionImage-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkIOImport-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkIOVideo-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkIOInfovis-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkChartsCore-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkIOParallelXML-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkRenderingLOD-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkIOPLY-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkViewsInfovis-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkverdict-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkoggtheora-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkRenderingOpenGL2-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkglew-6.3.so.1
proj6B: /usr/lib/x86_64-linux-gnu/libSM.so
proj6B: /usr/lib/x86_64-linux-gnu/libICE.so
proj6B: /usr/lib/x86_64-linux-gnu/libX11.so
proj6B: /usr/lib/x86_64-linux-gnu/libXext.so
proj6B: /usr/lib/x86_64-linux-gnu/libXt.so
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkjsoncpp-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkexoIIc-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkIONetCDF-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkNetCDF_cxx-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkNetCDF-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkhdf5_hl-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkhdf5-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkFiltersParallel-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtklibxml2-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkParallelCore-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkIOLegacy-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkIOXML-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkIOXMLParser-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkexpat-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkIOGeometry-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkInfovisLayout-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkViewsCore-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkInteractionWidgets-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkFiltersHybrid-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkFiltersModeling-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkImagingHybrid-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkIOImage-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkIOCore-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtktiff-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkjpeg-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkDICOMParser-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkmetaio-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkpng-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkRenderingAnnotation-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkImagingColor-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkRenderingVolume-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkInteractionStyle-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkRenderingLabel-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkFiltersImaging-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkImagingGeneral-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkImagingSources-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkChartsCore-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkRenderingContext2D-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkRenderingFreeType-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkRenderingCore-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkFiltersSources-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkCommonColor-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkFiltersGeometry-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkfreetype-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkzlib-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkInfovisCore-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkFiltersExtraction-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkFiltersGeneral-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkCommonComputationalGeometry-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkFiltersCore-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkFiltersStatistics-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkImagingFourier-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkImagingCore-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkCommonExecutionModel-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkCommonDataModel-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkCommonMisc-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkCommonSystem-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtksys-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkCommonTransforms-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkCommonMath-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkCommonCore-6.3.so.1
proj6B: /home/cheen/Desktop/410/VTK-build/lib/libvtkalglib-6.3.so.1
proj6B: CMakeFiles/proj6B.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable proj6B"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/proj6B.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/proj6B.dir/build: proj6B
.PHONY : CMakeFiles/proj6B.dir/build

CMakeFiles/proj6B.dir/requires: CMakeFiles/proj6B.dir/proj6.cxx.o.requires
.PHONY : CMakeFiles/proj6B.dir/requires

CMakeFiles/proj6B.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/proj6B.dir/cmake_clean.cmake
.PHONY : CMakeFiles/proj6B.dir/clean

CMakeFiles/proj6B.dir/depend:
	cd /home/cheen/mc410 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/cheen/mc410 /home/cheen/mc410 /home/cheen/mc410 /home/cheen/mc410 /home/cheen/mc410/CMakeFiles/proj6B.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/proj6B.dir/depend
