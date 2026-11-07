include_guard(GLOBAL)

function(collie_apply_compiler_options target)
  target_compile_features("${target}" PUBLIC cxx_std_17)

  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options("${target}" PRIVATE -Wall -Wextra)
  endif()

  if(COLLIE_ENABLE_SANITIZERS)
    if(NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
      message(FATAL_ERROR
        "COLLIE_ENABLE_SANITIZERS requires a GNU- or Clang-compatible compiler"
      )
    endif()

    target_compile_options("${target}" PRIVATE
      -fsanitize=address,undefined
      -fno-omit-frame-pointer
    )
    target_link_options("${target}" PRIVATE
      -fsanitize=address,undefined
      -fno-omit-frame-pointer
    )
  endif()
endfunction()
