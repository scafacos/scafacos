
#define MSEG_BORDER_UPDATE_REDUCTION

#define MSEG_INFO

/*#define MSS_ROOT*/


/* do reduce+bcast instead of allreduce on jugene */
#ifdef JUGENE
# define GLOBAL_REDUCEBCAST_THRESHOLD  0
#endif
