/* GPT: Memory manager */

/* Passed as parameter to notifyfunc's in case of a synchronization event */
typedef enum gptsyncstatus
{
  gptcreatearray,
  gptremovearray,
  gptremoveset,

  gptadditem,
  gptremoveitem
} gptsyncstatus ;

/* Encapsulation of an array to be synchronized, must be in a set */
typedef struct gptsyncarray /* sa */
{
  int itemcount ;
  size_t blocksize ;
  void *memory ;

  void (*notifyfunc)(struct gptsyncarray *psa, gptsyncstatus stat, int n, void *info) ;
  void *info ;
  struct gptsyncarray *next ;
} gptsyncarray ;

/* Encapsulation of a set of arrays to be synchronized */
typedef struct gptsyncset /* sp */
{
  int itemcount ;
  gptsyncarray *gptsyncarrays ;

  int chunkcount ;
  struct gptsyncset *next ;
} gptsyncset ;

/* Standard memory management with automatic GPT error checking */
void *gptmalloc(size_t size) ;
void *gptrealloc(void *mem, size_t size) ;
void *gptcalloc(size_t num, size_t size) ;
void gptfree(void *mem) ;

/* item array functions */
void *gptmallocitemarray(size_t blocksize, int chunkcount, int len) ;
void *gptmallocitem(void *mem, size_t blocksize, int chunkcount, int len) ;
void *gptfreeitem(void *mem, size_t blocksize, int chunkcount, int len, int n) ;
void gptfreeitemarray(void *mem, size_t blocksize, int chunkcount, int len) ;

/* Synchronized set and array functions */
gptsyncset *gptcreatesyncset(int itemcount, int chunkcount) ;
int gptremovesyncset(gptsyncset *psp) ;
gptsyncarray *gptcreatesyncarray(gptsyncset *psp, size_t blocksize,
        void (*notifyfunc)(gptsyncarray *psa, gptsyncstatus stat, int n, void *info),
        void *info) ;
int gptremovesyncarray(gptsyncset *psp, gptsyncarray *psa) ;
void gptcreatesyncitem(gptsyncset *psp) ;
void gptremovesyncitem(gptsyncset *psp, int n) ;

/* Block malloc functions */
struct gptmemblock
{
  char *memstart ; // Start of BIG memory chunk
  char *memend ;   // End of memory chunk, value indicates first unavailable byte.
  char *memnext ;  // Start of first available byte/End of smallest possible BIG block
  int overshoot ;
} ;

void gptmemblockinit(struct gptmemblock *pmem, size_t size, int overshoot) ;
void gptmemblockresize(struct gptmemblock *pmem) ;
void gptmemblockdone(struct gptmemblock *pmem) ;

void *gptblockmalloc(struct gptmemblock *pmem, size_t size) ;
void *gptblockcalloc(struct gptmemblock *pmem, size_t num, size_t size) ;
void gptblockfree(struct gptmemblock *pmem,void *buf) ;
