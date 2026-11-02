
################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2026 NCrystal developers                                   ##
##                                                                            ##
##  Licensed under the Apache License, Version 2.0 (the "License");           ##
##  you may not use this file except in compliance with the License.          ##
##  You may obtain a copy of the License at                                   ##
##                                                                            ##
##      http://www.apache.org/licenses/LICENSE-2.0                            ##
##                                                                            ##
##  Unless required by applicable law or agreed to in writing, software       ##
##  distributed under the License is distributed on an "AS IS" BASIS,         ##
##  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  ##
##  See the License for the specific language governing permissions and       ##
##  limitations under the License.                                            ##
##                                                                            ##
################################################################################

cmake_minimum_required(VERSION 3.16...3.31)

project( ExampleNCrystalProject VERSION 0.0.1 LANGUAGES CXX C )

#Find NCrystal:

#Find NCrystal directly if possible (e.g. via CMAKE_PREFIX_PATH or NCrystal_DIR,
#or when this project is itself built with scikit-build-core against
#ncrystal-core from Python wheels), and only otherwise ask "ncrystal-config
#--show cmakedir" (which is not always possible, e.g. when cross-compiling):
find_package( NCrystal 4.4.7 QUIET )
if ( NOT NCrystal_FOUND )
  execute_process(
    COMMAND ncrystal-config --show cmakedir
    OUTPUT_VARIABLE "NCrystal_DIR" OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  find_package( NCrystal 4.4.7 REQUIRED )
endif()

#Add our own targets, and link to NCrystal::NCrystal as appropriate:

add_executable( "testapp" "main.cc" )
target_link_libraries( "testapp" NCrystal::NCrystal )
install (TARGETS "testapp" RUNTIME DESTINATION bin)
