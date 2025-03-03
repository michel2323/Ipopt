// Copyright (C) 2004, 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPZEROMATRIX_HPP__
#define __IPZEROMATRIX_HPP__

#include "IpUtils.hpp"
#include "IpMatrix.hpp"

namespace Ipopt
{

/** Class for Matrices with only zero entries.
 */
class ZeroMatrix: public Matrix
{
public:
   /**@name Constructors / Destructors */
   //@{
   /** Constructor, taking the corresponding matrix space. */
   ZeroMatrix(
      const MatrixSpace* owner_space
   );

   /** Destructor */
   ~ZeroMatrix();
   //@}

protected:
   /**@name Methods overloaded from matrix */
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

   virtual void ComputeRowAMaxImpl(
      Vector& rows_norms,
      bool    init
   ) const
   { }

   virtual void ComputeColAMaxImpl(
      Vector& cols_norms,
      bool    init
   ) const
   { }

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
    * they will not be implicitly created/called.
    */
   //@{
   /** Default Constructor */
   ZeroMatrix();

   /** Copy Constructor */
   ZeroMatrix(
      const ZeroMatrix&
   );

   /** Default Assignment Operator */
   void operator=(
      const ZeroMatrix&
   );
   //@}
};

/** Class for matrix space for ZeroMatrix. */
class ZeroMatrixSpace: public MatrixSpace
{
public:
   /** @name Constructors / Destructors */
   //@{
   /** Constructor, given the number of row and columns. */
   ZeroMatrixSpace(
      Index nrows,
      Index ncols
   )
      : MatrixSpace(nrows, ncols)
   { }

   /** Destructor */
   virtual ~ZeroMatrixSpace()
   { }
   //@}

   virtual Matrix* MakeNew() const
   {
      return MakeNewZeroMatrix();
   }

   /** Method for creating a new matrix of this specific type. */
   ZeroMatrix* MakeNewZeroMatrix() const
   {
      return new ZeroMatrix(this);
   }

private:
   /**@name Default Compiler Generated Methods
    * (Hidden to avoid implicit creation/calling).
    * These methods are not implemented and
    * we do not want the compiler to implement
    * them for us, so we declare them private
    * and do not define them. This ensures that
    * they will not be implicitly created/called.
    */
   //@{
   /** Default Constructor */
   ZeroMatrixSpace();

   /** Copy Constructor */
   ZeroMatrixSpace(
      const ZeroMatrixSpace&
   );

   /** Default Assignment Operator */
   void operator=(
      const ZeroMatrixSpace&
   );
   //@}
};
} // namespace Ipopt
#endif
