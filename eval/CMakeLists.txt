#foreach(X RANGE 1 4)
#    add_executable(test_base_set${X} test.c)
#    target_link_libraries(test_base_set${X} base_set${X})
#    target_compile_definitions(test_base_set${X} PUBLIC SET${X})
#endforeach()

foreach(X RANGE 1 6)
    add_executable(measure_speed_base_set${X} measure_speed.c)
    set_target_properties(measure_speed_base_set${X} PROPERTIES LINK_DEPENDS ../link.ld)
    target_link_libraries(measure_speed_base_set${X} base_set${X})
    target_compile_definitions(measure_speed_base_set${X} PUBLIC SET${X})
endforeach()

foreach(X RANGE 1 6)
    add_executable(measure_speed_masked_core_set${X} measure_speed.c)
    set_target_properties(measure_speed_masked_core_set${X} PROPERTIES LINK_DEPENDS ../link.ld)
    target_link_libraries(measure_speed_masked_core_set${X} masked_core_set${X})
    target_compile_definitions(measure_speed_masked_core_set${X} PUBLIC SET${X})
endforeach()

foreach(X RANGE 1 6)
    add_executable(measure_speed_masked_ext_set${X} measure_speed.c)
    set_target_properties(measure_speed_masked_ext_set${X} PROPERTIES LINK_DEPENDS ../link.ld)
    target_link_libraries(measure_speed_masked_ext_set${X} masked_ext_set${X})
    target_compile_definitions(measure_speed_masked_ext_set${X} PUBLIC SET${X})
endforeach()

foreach(X RANGE 1 6)
    add_executable(measure_speed_masked_3_set${X} measure_speed.c)
    set_target_properties(measure_speed_masked_3_set${X} PROPERTIES LINK_DEPENDS ../link.ld)
    target_link_libraries(measure_speed_masked_3_set${X} masked_3_set${X})
    target_compile_definitions(measure_speed_masked_3_set${X} PUBLIC SET${X})
endforeach()

add_executable(ttest ttest.c)
set_target_properties(ttest PROPERTIES LINK_DEPENDS ../link.ld)
target_link_libraries(ttest masked_ext_set1)
target_compile_definitions(ttest PUBLIC SET1)

#add_executable(masked_core-test test.c)
#target_link_libraries(masked_core-test masked_core)



