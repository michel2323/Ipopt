#include "IpPetscSolverInterface.hpp"
//#include "cholmod.h"

namespace Ipopt
{

PetscSolverInterface::PetscSolverInterface()
{
   DBG_START_METH("PetscSolverInterface::PetscSolverInterface()",
                  dbg_verbosity);
   int ierr = PetscInitialize(NULL, NULL, (char*)0, NULL);
   ierr = KSPCreate(PETSC_COMM_SELF,&ctx.ksp);
   //KSPSetType(ctx.ksp, KSPCG);

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
   PetscInt *pia;
   PetscInt *pja;
   if (values == NULL) delete [] values;
   values = new double[nonzeros];
   pia = new PetscInt[nonzeros];
   pja = new PetscInt[nonzeros];
   for (int i = 0; i < nonzeros; i++) pia[i] = ia[i];
   for (int i = 0; i < nonzeros; i++) pja[i] = ja[i];
   MatCreateSeqSBAIJWithArrays(PETSC_COMM_SELF, 1, dim, dim, pia, pja, values, &ctx.A);
   //MatCreateSeqSBAIJWithArrays(PETSC_COMM_SELF, 1, dim, dim, const_cast<PetscInt*>(ia), const_cast<PetscInt*>(ja), values, &ctx.A);
   //MatCreateSeqAIJWithArrays(PETSC_COMM_SELF, dim, dim, const_cast<PetscInt*>(ia), const_cast<PetscInt*>(ja), values, &ctx.A);

   VecCreateSeqWithArray(PETSC_COMM_SELF,1,dim,NULL,&ctx.b);
   VecCreateSeqWithArray(PETSC_COMM_SELF,1,dim,NULL,&ctx.x);
   this->dim = dim;
   //std::cout << "dim: " << dim << std::endl;
   //have_symbolic_factorization_ = false;
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
   Mat mat;
   //if (usedirect > 0) {
   //MatConvert(ctx.A, MATSEQAIJCUSPARSE, MAT_INITIAL_MATRIX, &mat); 
     MatConvert(ctx.A, MATSEQAIJ, MAT_INITIAL_MATRIX, &mat); 
     KSPSetOperators(ctx.ksp, mat, mat);
     KSPGetPC(ctx.ksp,&ctx.pc);
     PCSetType(ctx.pc,PCCHOLESKY);
     PCFactorSetMatSolverType(ctx.pc, MATSOLVERCHOLMOD);
     //PCFactorSetMatSolverType(ctx.pc, MATSOLVERCUSPARSE);
     usedirect = 0;
     //}
   //else {
     //KSPSetOperators(ctx.ksp, ctx.A, ctx.A);
     //KSPGetPC(ctx.ksp,&ctx.pc);
     //PCSetType(ctx.pc,PCNONE);
   //}
   for (Index i = 0; i < nrhs; i++) {
     VecPlaceArray(ctx.b, rhs_vals + i*dim);
     VecPlaceArray(ctx.x, sol);
     KSPSolve(ctx.ksp, ctx.b, ctx.x);
     PetscInt its;
     KSPGetIterationNumber(ctx.ksp, &its);
     std::cout << "Iterations: " << its << std::endl;
     VecGetArray(ctx.x, &sol);
     for (int j = 0; j < dim ; j++) {
         rhs_vals[i*dim + j] = sol[j];
     }
     VecRestoreArray(ctx.x, &sol);
   }

   //PetscInt nzero, npos;
   PetscInt nzero, npos;
   if (check_NegEVals && ProvidesInertia()) {
     PetscInt pnumberOfNegEVals = numberOfNegEVals;
     MatGetInertia(ctx.A, &pnumberOfNegEVals, &nzero, &npos);
     //MatGetInertia(ctx.A, &numberOfNegEVals, &nzero, &npos);
     negevals_ = pnumberOfNegEVals;
     //negevals_ = numberOfNegEVals;
   }
   VecResetArray(ctx.x);
   VecResetArray(ctx.b);
   MatDestroy(&mat);
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
   std::cout << "asking for better quality\n" << std::endl;
   PetscReal rtol;
   PetscReal abstol;
   PetscReal dtol;
   PetscInt maxits;
   //if (usedirect == 0) {
     //usedirect++;
     //return true;
   //}

   //KSPGetTolerances(ctx.ksp, &rtol, &abstol, &dtol, &maxits);
   //std::cout << "rtol: " << rtol << " abstol: " << abstol << " dtol: " << dtol << " maxits: " << maxits << std::endl;
   //maxits*=2;
   //KSPSetTolerances(ctx.ksp, rtol, abstol, dtol, maxits);
   return false;
}
}
