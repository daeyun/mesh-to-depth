add_library(mesh2depth SHARED mesh2depth.cc)
target_link_libraries(mesh2depth eigen3 pthread common depth_tracer ray_mesh_intersection camera)
