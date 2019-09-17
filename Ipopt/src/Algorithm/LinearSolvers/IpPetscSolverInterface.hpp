#ifndef __IPPETSCSOLVERINTERFACE_HPP__
#define __IPPETSCSOLVERINTERFACE_HPP__

#include "IpSparseSymLinearSolverInterface.hpp"
// extern "C" {
#include "petscksp.h"
#include "petscmat.h"
// }

namespace Ipopt
{
typedef struct {
  Vec x, b;             /* solution vector, right-hand-side vector */
  Mat A;                /* sparse matrix */
  KSP ksp;              /* linear solver context */
  PC pc;
} UserCtx;

/** Interface to the linear solver Petsc, derived from
 *  SparseSymLinearSolverInterface and MumpsSolverInterface.
 */
class PetscSolverInterface: public SparseSymLinearSolverInterface
{
public:
   /** @name Constructor/Destructor */
   //@{
   /** Constructor */
   UserCtx ctx;
   double* values = NULL;
   int dim;        /* grid dimensions */
   int usedirect = 0;
   PetscSolverInterface();

//    /** Destructor */
//    virtual ~MumpsSolverInterface();
//    //@}

   bool InitializeImpl(
      const OptionsList& options,
      const std::string& prefix
   );

   /** @name Methods for requesting solution of the linear system. */
   //@{
   virtual ESymSolverStatus InitializeStructure(
      Index        dim,
      Index        nonzeros,
      const Index* airn,
      const Index* ajcn
   );

   virtual double* GetValuesArrayPtr();

   virtual ESymSolverStatus MultiSolve(
      bool         new_matrix,
      const Index* airn,
      const Index* ajcn,
      Index        nrhs,
      double*      rhs_vals,
      bool         check_NegEVals,
      Index        numberOfNegEVals
   );

   virtual Index NumberOfNegEVals() const;
//    //@}

   //* @name Options of Linear solver */
   //@{
   virtual bool IncreaseQuality();

   virtual bool ProvidesInertia() const
   {
      return false;
   }

   EMatrixFormat MatrixFormat() const
   {
      return CSR_Format_0_Offset;
   }
   //@}

//    static void RegisterOptions(
//       SmartPtr<RegisteredOptions> roptions
//    );

//    virtual bool ProvidesDegeneracyDetection() const;

//    virtual ESymSolverStatus DetermineDependentRows(
//       const Index*      ia,
//       const Index*      ja,
//       std::list<Index>& c_deps
//    );

private:
   /**@name Default Compiler Generated Methods
    * (Hidden to avoid implicit creation/calling).
    * These methods are not implemented and
    * we do not want the compiler to implement
    * them for us, so we declare them private
    * and do not define them. This ensures that
    * they will not be implicitly created/called. */
   //@{
   /** Copy Constructor */
   PetscSolverInterface(
      const PetscSolverInterface&
   );

   /** Default Assignment Operator */
   void operator=(
      const PetscSolverInterface&
   );
   //@}


   /** @name Information about most recent factorization/solve */
   //@{
   /** Number of negative eigenvalues */
   Index negevals_;
   //@}

   /** @name Initialization flags */
   //@{
   /** Flag indicating if internal data is initialized.
    *  For initialization, this object needs to have seen a matrix.
    */
   bool initialized_;
   /** Flag indicating if the matrix has to be refactorized because
    *  the pivot tolerance has been changed.
    */
   bool pivtol_changed_;
   /** Flag that is true if we just requested the values of the
    *  matrix again (SYMSOLVER_CALL_AGAIN) and have to factorize
    *  again.
    */
   bool refactorize_;

//    /** @name Solver specific data/options */
//    //@{
//    /** Pivot tolerance */
//    Number pivtol_;

//    /** Maximal pivot tolerance */
//    Number pivtolmax_;

//    /** Percent increase in memory */
//    Index mem_percent_;

//    /** Permutation and scaling method in MUMPS */
//    Index mumps_permuting_scaling_;

//    /** Pivot order in MUMPS. */
//    Index mumps_pivot_order_;

//    /** Scaling in MUMPS */
//    Index mumps_scaling_;

//    /** Threshold in MUMPS to state that a constraint is linearly dependent */
//    Number mumps_dep_tol_;

   /** Flag indicating whether the TNLP with identical structure has
    *  already been solved before.
    */
   bool warm_start_same_structure_;
   //@}

   /** Flag indicating if symbolic factorization has already been called */
   bool have_symbolic_factorization_;

//    /** @name Internal functions */
//    //@{
//    /** Call MUMPS (job=1) to perform symbolic manipulations, and reserve
//     *  memory.
//     */
//    ESymSolverStatus SymbolicFactorization();

//    /** Call MUMPS (job=2) to factorize the Matrix.
//     *  It is assumed that the first nonzeros_ element of a_ contain the values
//     *  of the matrix to be factorized.
//     */
//    ESymSolverStatus Factorization(
//       bool  check_NegEVals,
//       Index numberOfNegEVals
//    );

//    /** Call MUMPS (job=3) to do the solve. */
//    ESymSolverStatus Solve(
//       Index   nrhs,
//       double* rhs_vals
//    );
//    //@}
};

} // namespace Ipopt
#endif
