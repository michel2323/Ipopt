#include "IpPetscSolverInterface.hpp"

namespace Ipopt
{

PetscSolverInterface::PetscSolverInterface()
{
   DBG_START_METH("PetscSolverInterface::PetscSolverInterface()",
                  dbg_verbosity);
   int ierr = PetscInitialize(NULL, NULL, (char*)0, NULL);
   ierr = KSPCreate(PETSC_COMM_SELF,&ctx.ksp);

}

bool PetscSolverInterface::InitializeImpl(
   const OptionsList& options,
   const std::string& prefix
)
{
   // Reset all private data
   initialized_ = false;
   pivtol_changed_ = false;
   refactorize_ = false;
   have_symbolic_factorization_ = false;

   return true;
}

ESymSolverStatus PetscSolverInterface::InitializeStructure(
   Index        dim,
   Index        nonzeros,
   const Index* ia,
   const Index* ja
)
{
   DBG_START_METH("PetscSolverInterface::InitializeStructure", dbg_verbosity);

   ESymSolverStatus retval = SYMSOLVER_SUCCESS;

   if (values == NULL) delete [] values;
   values = new double[nonzeros];
   MatCreateSeqSBAIJWithArrays(PETSC_COMM_SELF, 1, dim, dim, const_cast<PetscInt*>(ia), const_cast<PetscInt*>(ja), values, &ctx.A);
   VecCreateSeqWithArray(PETSC_COMM_SELF,1,dim,NULL,&ctx.b);
   VecCreateSeqWithArray(PETSC_COMM_SELF,1,dim,NULL,&ctx.x);
   this->dim = dim;
   have_symbolic_factorization_ = false;
   initialized_ = true;
   return retval;
}

double* PetscSolverInterface::GetValuesArrayPtr()
{
   DBG_START_METH("PetscSolverInterface::GetValuesArrayPtr", dbg_verbosity)
   DBG_ASSERT(initialized_);
   return values;
}

ESymSolverStatus PetscSolverInterface::MultiSolve(
   bool         new_matrix,
   const Index* ia,
   const Index* ja,
   Index        nrhs,
   double*      rhs_vals,
   bool         check_NegEVals,
   Index        numberOfNegEVals
)
{
   DBG_START_METH("PetscSolverInterface::MultiSolve", dbg_verbosity);
   DBG_ASSERT(!check_NegEVals || ProvidesInertia());
   DBG_ASSERT(initialized_);
   double *sol = new double[dim];
   KSPSetOperators(ctx.ksp, ctx.A, ctx.A);
   for (Index i = 0; i < nrhs; i++) {
     VecPlaceArray(ctx.b, rhs_vals + i*dim);
     VecPlaceArray(ctx.x, sol);
     KSPSolve(ctx.ksp, ctx.b, ctx.x);
     VecGetArray(ctx.x, &sol);
     for (int j = 0; j < dim ; j++) {
         rhs_vals[i*dim + j] = sol[j];
     }
     VecRestoreArray(ctx.x, &sol);
   }

   int nzero, npos;
   if (check_NegEVals && ProvidesInertia()) {
     MatGetInertia(ctx.A, &numberOfNegEVals, &nzero, &npos);
     negevals_ = numberOfNegEVals;
   }
   VecResetArray(ctx.x);
   VecResetArray(ctx.b);
   return SYMSOLVER_SUCCESS;
}

Index PetscSolverInterface::NumberOfNegEVals() const
{
   DBG_START_METH("PetscSolverInterface::NumberOfNegEVals", dbg_verbosity);
   DBG_ASSERT(negevals_ >= 0);
   return negevals_;
}

bool PetscSolverInterface::IncreaseQuality()
{
   DBG_START_METH("MumpsTSolverInterface::IncreaseQuality", dbg_verbosity);
   return false;
}
}