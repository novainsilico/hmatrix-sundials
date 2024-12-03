/* Check function return value...
    opt == 0 means SUNDIALS function allocates memory so check if
             returned NULL pointer
    opt == 1 means SUNDIALS function returns a flag so check if
             flag >= 0
    opt == 2 means function allocates memory so check if returned
             NULL pointer
*/
int check_flag(void *flagvalue, const char *funcname, int opt, void (*report_error)(int,const char*, const char*, char*, void*));

/* Give void a name so that we can map it to the Haskell type UserData */
typedef void UserData;

/* The use of this macro requires that a function

   void debug(char*)

   is in scope */
#ifdef ENABLE_DEBUG
#define DEBUG(...) do { \
  char *formatted_str; \
  if (asprintf(&formatted_str, __VA_ARGS__)) { \
    debug(formatted_str); \
    free(formatted_str); \
  } } while (0)
#else
#define DEBUG(...) do {} while (0)
#endif
