# libgmultigrid
A header-only library for geometric multigrid.

## Installation

The recommended way to use this library is by adding it as a submodule to the project that you would like to use it in. Navigate to the directory you would like to clone the submodule into, and then run:
```
git submodule add https://github.com/icethrush/libgmultigrid.git
```
After this, add `libgmultigrid/include` to the include path of your project, and you should be able to use it in your code.

## Dependencies

This library uses Eigen, and thus requires the Eigen directories to be in the include path of your project or machine.

This library has only been included in projects that were built using Clang. Usage with gcc/g++ is possible, but the compiler may emit a different set of warnings that may need to be suppressed (if, for instance, you are compiling with -Wall).

## How to use

Let's suppose you want to solve a system `Ax = b`, and you have some way of performing matrix-vector products with the matrix `A`. To use this library, you need to implement a few classes with specific signatures.

### 1. Implement a vector multiplier

First, you must define a class `Mult` that extends `VectorMultiplier<Mult>` (i.e. a recursive template class) and implements a template method `Multiply` with the following signature (in `vector_multiplier.h`):

```
        template<typename V, typename Dest>
        void Multiply(V &v, Dest &b) const {
            ...
        }
```

The specific signatures of the template arguments `V` and `Dest` here are defined internally by Eigen, but generally, one can treat `V` as a read-only `VectorXd`, and `Dest` as a write-only `VectorXd`. As a simple example, here is a class that wraps a dense matrix from Eigen:

```
    class DenseMatrixMult : public VectorMultiplier<DenseMatrixMult> {
        public:
        DenseMatrixMult(Eigen::MatrixXd mat) {
            A = mat;
        }

        template<typename V, typename Dest>
        void Multiply(V &v, Dest &b) const;

        private:
        Eigen::MatrixXd A;

    };

    template<typename V, typename Dest>
    void DenseMatrixMult::Multiply(V &v, Dest &b) const {
        Eigen::VectorXd x;
        int nRows = v.size();
        x.setZero(nRows);

        for (int i = 0; i < nRows; i++) {
            x(i) = v[i];
        }

        Eigen::VectorXd Ax = A * x;

        for (int i = 0; i < nRows; i++) {
            b[i] = Ax(i);
        }
    }
```

Of course, for efficiency reasons, one would want to write a more complex multiplier that performs some sort of approximation. The `VectorMultiplier` can be arbitrarily complex, as long as it has a method called `Multiply` that matches the above signature.

### 2. Implement a prolongation operator

Next, you need to implement a class that extends `MultigridOperator` (in `multigrid_operator.h`). This will serve as the prolongation operator, which is what maps data between the different levels of the multigrid hierarchy. `MultigridOperator` has three virtual methods:

```
        virtual Eigen::VectorXd prolong(Eigen::VectorXd v) = 0;
        virtual Eigen::VectorXd restrictWithTranspose(Eigen::VectorXd v) = 0;
        virtual Eigen::VectorXd restrictWithPinv(Eigen::VectorXd v) = 0;
```




### 3. Implement the multigrid domain

The class `MultigridDomain` defines an interface with several operations relating to the geometric domain (e.g. a square grid, a space curve, a triangle mesh, etc.).


