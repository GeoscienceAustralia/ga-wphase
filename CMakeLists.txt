cmake_minimum_required(VERSION 3.17.2...3.26)

project(_libtau
  VERSION 0.3.5
  DESCRIPTION "libtau travel time calculator"
  LANGUAGES C Fortran
)

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module NumPy)

# F2PY headers
execute_process(
  COMMAND "${Python_EXECUTABLE}"
  -c "import numpy.f2py; print(numpy.f2py.get_include())"
  OUTPUT_VARIABLE F2PY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

set(f2py_module_name "_libtau")
set(f2py_module_c "${f2py_module_name}module.c")

add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_c}"
  COMMAND ${Python_EXECUTABLE} -m "numpy.f2py"
                   "${CMAKE_SOURCE_DIR}/wphase/psi/_libtau.pyf"
                   --lower
  DEPENDS
    wphase/psi/libtau.f
    wphase/psi/ttimes_subrout.f
    wphase/psi/_libtau.pyf
)

python_add_library(
  ${CMAKE_PROJECT_NAME} MODULE
  "${f2py_module_c}"
  "${F2PY_INCLUDE_DIR}/fortranobject.c"
  "${CMAKE_SOURCE_DIR}/wphase/psi/ttimes_subrout.f"
  "${CMAKE_SOURCE_DIR}/wphase/psi/libtau.f"
  WITH_SOABI
)

target_include_directories(${CMAKE_PROJECT_NAME} PUBLIC
                           ${F2PY_INCLUDE_DIR}
                          )

target_link_libraries(${CMAKE_PROJECT_NAME} PUBLIC Python::NumPy)

install(TARGETS ${CMAKE_PROJECT_NAME} DESTINATION wphase/psi)
