cmake_minimum_required(VERSION 3.18)

set(CMAKE_C_COMPILER "riscv32-unknown-elf-gcc")

project(mpqc LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/libs/)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/libs/)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(mpqc_SRCS
  src/masked_poly_frommsg.c
  src/frodo_sampler.c
  src/haetae_sampler.c
  src/hawk_sampler.c
  src/gadgets.c
  src/masked_bin2tri.c
  src/masked_cmp.c
  src/masked_fisheryates.c
  src/masked_i2c.c
  src/masked_reject.c
  src/masked_repand.c
  src/masked_keccak.c
  src/masked_representations.c
  src/masked_sorting.c
  src/randombytes.c
  src/bike_rotr.c
)
include_directories(include)

add_compile_options(-march=rv32g -O3 -ffreestanding -static -ffunction-sections)
add_link_options(-ffreestanding -static -ffunction-sections -Tlink.ld)
#add_compile_options(-Wall -Wextra -Wpedantic -Wno-maybe-uninitialized)
#add_compile_options(-Wmissing-prototypes -Wredundant-decls -Wshadow -Wpointer-arith)
#add_compile_options(-O3 -fomit-frame-pointer)

foreach(X RANGE 1 6)
  add_library(base_set${X} ${mpqc_SRCS})
  target_compile_definitions(base_set${X} PUBLIC SET${X})
endforeach()

foreach(X RANGE 1 6)
  add_library(masked_core_set${X} ${mpqc_SRCS})
  target_compile_definitions(masked_core_set${X} PUBLIC SET${X} USE_MASKED_ISA)
endforeach()

foreach(X RANGE 1 6)
  add_library(masked_ext_set${X} ${mpqc_SRCS})
  target_compile_definitions(masked_ext_set${X} PUBLIC SET${X} USE_MASKED_ISA USE_MASKED_EXT)
endforeach()

foreach(X RANGE 1 6)
  add_library(masked_3_set${X} ${mpqc_SRCS})
  target_compile_definitions(masked_3_set${X} PUBLIC SET${X} USE_MASKED_EXT)
endforeach()


add_subdirectory(eval)
