// Copyright (C) 2005, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Andreas Waechter             IBM    2005-12-24

#ifndef __IPDENSEGENMATRIX_HPP__
#define __IPDENSEGENMATRIX_HPP__

#include "IpUtils.hpp"
#include "IpMatrix.hpp"
#include "IpDenseVector.hpp"
#include "IpDenseSymMatrix.hpp"

namespace Ipopt
{

/** forward declarations */
class DenseGenMatrixSpace;

/** Class for dense general matrices.
 *
 * Matrix elements are stored in one array in "Fortran" format.
 */
class DenseGenMatrix: public Matrix
{
public:

   /**@name Constructors / Destructors */
   //@{
   /** Constructor, taking the owner_space. */
   DenseGenMatrix(
      const DenseGenMatrixSpace* owner_space
   );

   /** Destructor */
   ~DenseGenMatrix();
   //@}

   /** Create a new DenseGenMatrix from same MatrixSpace */
   SmartPtr<DenseGenMatrix> MakeNewDenseGenMatrix() const;

   /** Retrieve the array for storing the matrix elements.
    *
    *  This is the non-const version, and it is assume that afterwards
    *  the calling method will set all matrix elements.  The matrix
    *  elements are stored one column after each other.
    */
   Number* Values()
   {
      initialized_ = true;
      ObjectChanged();
      return values_;
   }

   /** Retrieve the array that stores the matrix elements.
    *
    *  This is the const version, i.e., read-only.  The matrix
    *  elements are stored one column after each other.
    */
   const Number* Values() const
   {
      DBG_ASSERT(initialized_);
      return values_;
   }

   /** Method for copying the content of another matrix into this matrix */
   void Copy(
      const DenseGenMatrix& M
   );

   /** Set this matrix to be a multiple of the identity matrix.
    *
    *  This assumes that this matrix is square.
    */
   void FillIdentity(
      Number factor = 1.
   );

   /** Method for scaling the columns of the matrix.
    *
    *  The scaling factors are given in form of a DenseVector.
    */
   void ScaleColumns(
      const DenseVector& scal_vec
   );

   /** Method for adding the product of two matrices to this matrix. */
   void AddMatrixProduct(
      Number                alpha,
      const DenseGenMatrix& A,
      bool                  transA,
      const DenseGenMatrix& B,
      bool                  transB,
      Number                beta
   );

   /** Method for adding a high-rank update to this matrix.
    *
    *  It computes M = alpha*V1^T V2 + beta*M, where V1 and V2 are
    *  MultiVectorMatrices.
    */
   void HighRankUpdateTranspose(
      Number                   alpha,
      const MultiVectorMatrix& V1,
      const MultiVectorMatrix& V2,
      Number                   beta
   );

   /** Method for computing the Cholesky factorization of a positive
    *  definite matrix.
    *
    *  The factor is stored in this matrix, as
    *  lower-triangular matrix, i.e., M = J * J^T.
    *
    *  @return false if the factorization could not be done, e.g., when
    *  the matrix is not positive definite.
    */
   bool ComputeCholeskyFactor(
      const DenseSymMatrix& M
   );

   /** Method for computing an eigenvalue decomposition of the given
    *  symmetrix matrix M.
    *
    *  On return, this matrix contains the
    *  eigenvalues in its columns, and Evalues contains the
    *  eigenvalues.
    *
    *  @return false, if there were problems during the computation
    */
   bool ComputeEigenVectors(
      const DenseSymMatrix& M,
      DenseVector&          Evalues
   );

   /** Method for performing one backsolve with an entire matrix on
    *  the right hand side, assuming that the this matrix is square
    *  and contains a lower triangular matrix.
    *
    *  The incoming right hand side B is overwritten with the solution
    *  X of op(A)*X = alpha*B. op(A) = A or op(A) = A^T.
    */
   void CholeskyBackSolveMatrix(
      bool            trans,
      Number          alpha,
      DenseGenMatrix& B
   ) const;

   /** Method for performing a solve of a linear system for one
    *  vector, assuming that this matrix contains the Cholesky factor
    *  for the linear system.
    *
    *  The vector b contains the right hand
    *  side on input, and contains the solution on output.
    */
   void CholeskySolveVector(
      DenseVector& b
   ) const;

   /** Method for performing a solve of a linear system for one
    *  right-hand-side matrix, assuming that this matrix contains the
    *  Cholesky factor for the linear system.
    *
    *  The matrix B contains the right hand sides on input, and contains
    *  the solution on output.
    */
   void CholeskySolveMatrix(
      DenseGenMatrix& B
   ) const;

   /** Method for computing the LU factorization of an unsymmetric matrix.
    *
    *  The factorization is done in place. */
   bool ComputeLUFactorInPlace();

   /** Method for using a previously computed LU factorization for a
    *  backsolve with a matrix on the rhs.
    */
   void LUSolveMatrix(
      DenseGenMatrix& B
   ) const;

   /** Method for using a previously computed LU fatorization for a
    *  backsolve with a single vector.
    */
   void LUSolveVector(
      DenseVector& b
   ) const;

protected:
   /**@name Overloaded methods from Matrix base class*/
   //@{
   virtual void MultVectorImpl(
      Number        alpha,
      const Vector& x,
      Number        beta,
      Vector&       y
   ) const;

   virtual void TransMultVectorImpl(
      Number        alpha,
      const Vector& x,
      Number        beta,
      Vector&       y
   ) const;

   virtual bool HasValidNumbersImpl() const;

   virtual void ComputeRowAMaxImpl(
      Vector& rows_norms,
      bool    init
   ) const;

   virtual void ComputeColAMaxImpl(
      Vector& cols_norms,
      bool    init
   ) const;

   virtual void PrintImpl(
      const Journalist&  jnlst,
      EJournalLevel      level,
      EJournalCategory   category,
      const std::string& name,
      Index              indent,
      const std::string& prefix
   ) const;
   //@}

private:
   /**@name Default Compiler Generated Methods
    * (Hidden to avoid implicit creation/calling).
    * These methods are not implemented and
    * we do not want the compiler to implement
    * them for us, so we declare them private
    * and do not define them. This ensures that
    * they will not be implicitly created/called. */
   //@{
   /** Default Constructor */
   DenseGenMatrix();

   /** Copy Constructor */
   DenseGenMatrix(
      const DenseGenMatrix&
   );

   /** Default Assignment Operator */
   void operator=(
      const DenseGenMatrix&
   );
   //@}

   const DenseGenMatrixSpace* owner_space_;

   /** Array for storing the matrix elements (one columns after each other) */
   Number* values_;

   /** Flag indicating whether the values_ array has been initialized */
   bool initialized_;

   /** Enum for factorization type */
   enum Factorization
   {
      NONE,
      LU,
      CHOL
   };

   /** Flag indicating if and which factorization has been applied */
   Factorization factorization_;

   /** Array for storing the pivot sequences if the matrix has been LU-factorized */
   Index* pivot_;
};

/** This is the matrix space for DenseGenMatrix.
 */
class DenseGenMatrixSpace: public MatrixSpace
{
public:
   /** @name Constructors / Destructors */
   //@{
   /** Constructor for matrix space for DenseGenMatrices.
    *
    * Takes in dimension of the matrices.
    */
   DenseGenMatrixSpace(
      Index nRows,
      Index nCols
   );

   /** Destructor */
   ~DenseGenMatrixSpace()
   { }
   //@}

   /** Method for creating a new matrix of this specific type. */
   DenseGenMatrix* MakeNewDenseGenMatrix() const
   {
      return new DenseGenMatrix(this);
   }

   virtual Matrix* MakeNew() const
   {
      return MakeNewDenseGenMatrix();
   }

};

inline SmartPtr<DenseGenMatrix> DenseGenMatrix::MakeNewDenseGenMatrix() const
{
   return owner_space_->MakeNewDenseGenMatrix();
}

} // namespace Ipopt
#endif
