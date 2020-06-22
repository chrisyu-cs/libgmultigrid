#pragma once

#include "multigrid_operator.h"
#include "nullspace_projector.h"
#include <Eigen/Core>

namespace LWS {

    template<typename Mult, typename Operator>
    class MultigridDomain {
        public:
        typedef Mult MultType;
        typedef Operator OperatorType;
        virtual ~MultigridDomain() {}

        virtual MultigridDomain<Mult, Operator>* Coarsen(Operator* prolongOp) const = 0;
        virtual Mult* GetMultiplier() const = 0;
        virtual Eigen::MatrixXd GetFullMatrix() const = 0;
        virtual Eigen::VectorXd DirectSolve(Eigen::VectorXd &b) const = 0;
        virtual int NumVertices() const = 0;
        virtual int NumRows() const  = 0;
        virtual Operator* MakeNewOperator() const = 0;
        virtual NullSpaceProjector* GetConstraintProjector() const = 0;
    };
}