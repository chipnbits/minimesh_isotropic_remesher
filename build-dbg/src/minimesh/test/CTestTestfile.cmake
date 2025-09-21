# CMake generated Testfile for 
# Source directory: /home/sghys/projects/CPSC524-modeling/src/minimesh/test
# Build directory: /home/sghys/projects/CPSC524-modeling/build-dbg/src/minimesh/test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(minimesh.properties "/home/sghys/projects/CPSC524-modeling/build-dbg/src/minimesh/test/minimesh_tests" "--test-case=*mesh*")
set_tests_properties(minimesh.properties PROPERTIES  WORKING_DIRECTORY "/home/sghys/projects/CPSC524-modeling" _BACKTRACE_TRIPLES "/home/sghys/projects/CPSC524-modeling/src/minimesh/test/CMakeLists.txt;29;add_test;/home/sghys/projects/CPSC524-modeling/src/minimesh/test/CMakeLists.txt;0;")
add_test(minimesh.subdiv "/home/sghys/projects/CPSC524-modeling/build-dbg/src/minimesh/test/minimesh_tests" "--test-case=*subdiv*")
set_tests_properties(minimesh.subdiv PROPERTIES  WORKING_DIRECTORY "/home/sghys/projects/CPSC524-modeling" _BACKTRACE_TRIPLES "/home/sghys/projects/CPSC524-modeling/src/minimesh/test/CMakeLists.txt;32;add_test;/home/sghys/projects/CPSC524-modeling/src/minimesh/test/CMakeLists.txt;0;")
add_test(minimesh.utils "/home/sghys/projects/CPSC524-modeling/build-dbg/src/minimesh/test/minimesh_tests" "--test-case=*utils*" "--success")
set_tests_properties(minimesh.utils PROPERTIES  WORKING_DIRECTORY "/home/sghys/projects/CPSC524-modeling" _BACKTRACE_TRIPLES "/home/sghys/projects/CPSC524-modeling/src/minimesh/test/CMakeLists.txt;35;add_test;/home/sghys/projects/CPSC524-modeling/src/minimesh/test/CMakeLists.txt;0;")
