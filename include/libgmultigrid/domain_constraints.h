#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace LWS {
    
    /* 
    An interface for any constraint that could be part of a saddle problem
    to be solved with multigrid. Note that there isn't any functionality in
    this library for automatically incorporating these constraints; this
    needs to be implemented by the user in their MultigridDomain.

    Five functions need to be implemented:
    - int NumConstraintRows() -- returns how many rows the constraints
        should occupy in the matrix.
    - int NumExpectedCols() -- returns how many columns (i.e. degrees of freedom)
        the constraints have.
    - void AddTriplets(std::vector<Eigen::Triplet<double>>&) -- adds the (row, col, val)
        triplets corresponding to the constraint matrix to the given list.
    - void SetTargetValues(Eigen::VectorXd&) -- given a vector with NumConstraintRows()
        entries, sets each entry to be the corresponding value of that constraint function.
    - void NegativeConstraintValues(Eigen::VectorXd &b, Eigen::VectorXd &targets) --
        given a vector b with NumConstraintRows() entries, and given the target
        values of the constraint function, fills b with the corresponding negative
        values of the function.
     */

    template<typename T>
    class DomainConstraints {
        public:
        // Should return how many rows will be occupied by the matrix block for
        // this constraint set.
        int NumConstraintRows() const {
            return static_cast<const T&>(*this).NumConstraintRows();
        }

        // Should return now many columns the matrix block for this constraint set 
        // will span, i.e. the total number of degrees of freedom.
        int NumExpectedCols() const {
            return static_cast<const T&>(*this).NumExpectedCols();
        }

        // Returns the dimensions of the full (square) saddle matrix, which looks like:
        // [ A    B^T ]
        // [ B    0   ]
        // where A is the kernel of the problem being solved, and B is the block
        // representing all of the constraints.
        int SaddleNumRows() const {
            return NumConstraintRows() + NumExpectedCols();
        }

        // Fills the given matrix with entries encoding this constraint set,
        // resizing it to be NumConstraintRows() x NumExpectedCols() if it wasn't already. 
        // No offsets are added; the sparse matrix contains only this constraint block
        // and nothing more.
        void FillConstraintMatrix(Eigen::SparseMatrix<double> &B) const {
            std::vector<Eigen::Triplet<double>> triplets;
            int nRows = NumConstraintRows();
            int nCols = NumExpectedCols();
            static_cast<const T&>(*this).AddTriplets(triplets);

            B.resize(nRows, nCols);
            B.setFromTriplets(triplets.begin(), triplets.end());
        }

        // Given a dense matrix of size SaddleNumRows() x SaddleNumRows() 
        // (which is the size of a full saddle matrix), fills the blocks
        // for B and B^T at the bottom and right edges of the matrix.
        // The input matrix must alerady have this size.
        void FillDenseBlock(Eigen::MatrixXd &A) const {
            std::vector<Eigen::Triplet<double>> triplets;
            static_cast<const T&>(*this).AddTriplets(triplets);
            int offset = NumExpectedCols();

            for (auto &t : triplets) {
                // Copy into lower-left block
                A(offset + t.row(), t.col()) = t.value();
                // Copy transpose into upper-right block
                A(t.col(), offset + t.row()) = t.value();
            }
        }
        
        // Given a vector containing the target values for the constraint function,
        // updates each target value with the current corresponding target value
        // for the constraint.
        void UpdateTargetValues(Eigen::VectorXd &targets) const {
            int nConstrs = NumConstraintRows();
            if (targets.rows() != nConstrs) {
                targets.setZero(nConstrs);
            }
            static_cast<const T&>(*this).SetTargetValues(targets);
        }

        // Fills the given vector b with the negated current values of this
        // constraint function. Entries are written starting from the index of
        // the argument "offset."
        double FillConstraintValues(Eigen::VectorXd &b, Eigen::VectorXd &targets, int offset) {
            // Fill values of constraint function in a generic way
            int nConstrs = NumConstraintRows();
            Eigen::VectorXd b_constrs(nConstrs);
            static_cast<const T&>(*this).NegativeConstraintValues(b_constrs, targets);
            b.block(offset, 0, nConstrs, 1) = b_constrs;

            return b_constrs.lpNorm<Eigen::Infinity>();
        }
    };
}