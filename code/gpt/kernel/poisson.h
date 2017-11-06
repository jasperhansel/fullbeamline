/* poisson.h: Interface of the Rostock Multigrid poisson solver
 *
 * Authors: Gisela Pöplau, Bas van der Geer, Marieke de Loos
 */

enum MGreturn
{
  NoError=0,			// OK
  NoMemory,			// Out of memory
  NoConvergence,		// Requested accuracy can not be achieved
  InvalidBoundaries,		// Invalid combination of boundary conditions
  NotImplementedYet,		// Requested features are not implemented yet
  NoGrids                       // No coarser grid have been constructed
} ;

enum MGboundary
{
  DirichletZero = 'D',		// V=0 at boundary edges
  Dirichlet     = 'A',		// V=approximated at boundary edge (Current=Brick approximation)
  Open          = 'O',		// V=0 at infinity
  Periodic      = 'P',		// Vlow=Vupp
  Ellipse       = 'E',		// Two out of three boundaries can be set to 'Ellipse' to simulate a beam pipe
  Neumann       = 'N'           // normal derivative at boundary edges=0
} ;

enum MGsolverType
{
  SolverMG,			// multigrid-algorithm is used
  SolverMGCG,			// Default:multigrid-preconditioned cg-algorithm is used
  SolverCG,			// preconditioned cg-algorithm is used if MG or MG_CG fails
  SolverSOR			// use SOR procedure, for comparisions only
} ;

struct MGinfo
{
  enum MGboundary low[3] ;	// Lower x,y and z boundary conditions
  enum MGboundary upp[3] ;	// Upper x,y and z boundary conditions

  enum MGsolverType solvertype ;
  struct gptmemblock mem ;
} ;

//enum MGreturn !!
int multigridpoisson283(int N[3], double *meshlines[3], double *F, double *U, double *error, struct MGinfo *bc ) ;
