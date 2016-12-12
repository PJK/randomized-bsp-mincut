file(GLOB KS_FILES ../karger-stein/*.cpp)
file(GLOB RC_FILES ../recursive-contract/*.cpp)
file(GLOB SOURCES ../*.cpp)


add_executable(square_root square_root.cpp ${SOURCES} ../input/InputIterator.cpp ${KS_FILES} ${RC_FILES})
target_link_libraries(square_root ${MPI_LIBRARIES})

add_executable(square_root_project_gigaedge square_root_project_gigaedge.cpp ${SOURCES} ../input/InputIterator.cpp ${KS_FILES} ${RC_FILES})
target_link_libraries(square_root_project_gigaedge ${MPI_LIBRARIES})

add_executable(transition transition.cpp ${SOURCES} ../input/InputIterator.cpp ${KS_FILES} ${RC_FILES})
target_link_libraries(transition ${MPI_LIBRARIES})

add_executable(boost_stoer_wagner boost_stoer_wagner.cpp ../input/InputIterator.cpp)
add_executable(karger_stein karger_stein.cpp ../input/InputIterator.cpp ${KS_FILES})

add_executable(sorting_test sorting_test.cpp ../utils.cpp)
target_link_libraries(sorting_test ${MPI_LIBRARIES})

add_executable(connectivity connectivity.cpp ../input/InputIterator.cpp)
add_executable(simple_generator simple_generator.cpp)
add_executable(levine_transform levine_transform.cpp)
add_executable(complete_graph_generator complete_graph_generator.cpp)