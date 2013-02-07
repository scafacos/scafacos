
#define DEBUG_HEADER(file)  \
    write(file,'("[PE ", I6.6, ", file: ",a, ", line ", I0, "] ")', advance='no') debug_my_rank, __FILE__, __LINE__


#define DEBUG_INFO(format, ...)                    \
        call debug_ipefile_open();                  \
          DEBUG_HEADER(debug_ipefile);              \
          write(debug_ipefile, format) __VA_ARGS__; \
        call debug_ipefile_close();

#define DEBUG_WARNING_ALL(format, ...)             \
        call debug_ipefile_open();                  \
          DEBUG_HEADER(debug_ipefile);              \
          write(debug_ipefile, format) __VA_ARGS__; \
        call debug_ipefile_close();                 \
        DEBUG_HEADER(debug_stdout);                 \
        write(debug_stdout, format) __VA_ARGS__;

#define DEBUG_WARNING(format, ...)                 \
        call debug_ipefile_open();                  \
          DEBUG_HEADER(debug_ipefile);              \
          write(debug_ipefile, format) __VA_ARGS__; \
        call debug_ipefile_close();                 \
        if (debug_my_rank == 0) DEBUG_HEADER(debug_stdout);              \
        if (debug_my_rank == 0) write(debug_stdout, format) __VA_ARGS__;

#define DEBUG_ERROR(format, ...)                   \
        call debug_ipefile_open();                  \
          DEBUG_HEADER(debug_ipefile);              \
          write(debug_ipefile, format) __VA_ARGS__; \
        call debug_ipefile_close();                 \
        DEBUG_HEADER(debug_stdout);                 \
        write(debug_stdout, format) __VA_ARGS__;    \
        call debug_mpi_abort();

#define DEBUG_ERROR_NO_HEADER(format, ...)         \
        call debug_ipefile_open();                  \
          write(debug_ipefile, format) __VA_ARGS__; \
        call debug_ipefile_close();                 \
        write(debug_stdout, format) __VA_ARGS__;    \
        call debug_mpi_abort();

#define DEBUG_ERROR_NO_DIAGFILE(format, ...)         \
        write(debug_stdout, format) __VA_ARGS__;    \
        call debug_mpi_abort();

#define DEBUG_DATA(format, ...)                    \
        call debug_ipefile_open();                  \
          write(debug_ipefile, format) __VA_ARGS__; \
        call debug_ipefile_close();


