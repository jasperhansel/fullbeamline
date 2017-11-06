/* pargs.h */

typedef struct _pa_entry {
  char opt;
  void *dst;
  void (*cnv)(char *arg,void *dst);
  char *arg;
  char *help;
} pa_entry;

int pargs(int *argc,char **argv[],pa_entry *pargtab);
void usage(pa_entry *tp,int mode,void (*print)(char *s));
void pa_bool(char *arg,void *dst);
void pa_int(char *arg,void *dst);
void pa_long(char *arg,void *dst);
void pa_dbl(char *arg,void *dst);
void pa_str(char *arg,void *dst);
void pa_strcpy(char *arg,void *dst);
void printvers(char *v,void (*print)(char *s));

void pa_error(char *arg,void *dst);
void pa_print(char *s);
