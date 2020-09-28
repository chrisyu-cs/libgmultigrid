# libgmultigrid
A header-only library for performing V-cycles as a part of geometric multigrid. 

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

Let's suppose you want to solve a system `Ax = b`, where `A` can be produced by some routine on a geometric domain (e.g. a curve, a mesh); a common example is the mesh Laplacian. Suppose also that you have some way of performing matrix-vector products with the matrix `A`. Then, the system can be solved using the multigrid method, for which this library provides a framework. To use this library, you will need to implement a few classes with specific signatures.

### 1. Implement a vector multiplier

First, you must define a class `Mult` that extends `VectorMultiplier<Mult>` (i.e. a recursive template class) and implements a template method `Multiply` with the following signature (in `vector_multiplier.h`):

```c++
    template<typename V, typename Dest>
    void Multiply(V &v, Dest &b) const {
        ...
    }
```

The specific signatures of the template arguments `V` and `Dest` here are defined internally by Eigen, but generally, one can treat `V` as a read-only `VectorXd`, and `Dest` as a write-only `VectorXd`. As a simple example, here is a class that wraps a dense matrix from Eigen:

```c++
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

```c++
virtual Eigen::VectorXd prolong(Eigen::VectorXd v) = 0;
virtual Eigen::VectorXd restrictWithTranspose(Eigen::VectorXd v) = 0;
virtual Eigen::VectorXd restrictWithPinv(Eigen::VectorXd v) = 0;
```

Intuitively, the prolongation operator is a sparse matrix `J` of size N x M, where N is the size of the finer level, and M is the size of the coarser level; multiplying by the matrix then maps data from the coarse level to the fine level. The entries of this matrix represent some sort of interpolation of the data from the coarse grid onto the cells of the fine grid. For instance, if the multigrid hierarchy was constructed by subdivision of triangle meshes, then `J` should be the subdivision matrix at each level.

The method `prolong` should do the multiplication with `P` and return the result. Similarly, `restrictWithTranspose` should multiply by the transpose of `J`, and `restrictWithPinv` should multiply by the pseudoinverse of `J`; note that a convenience function `ApplyPinv` is provided that can do this multiplication, assuming `J` is an Eigen `SparseMatrix`.

Note, however, if you have constraints in your system, then things become more complicated. Explicit constraint values often cannot be prolonged just by using `J`, and must be handled separately. Methods such as that of [Braess and Sarazin](https://homepage.ruhr-uni-bochum.de/dietrich.braess/smooth.pdf) also exist that do not require explicit constraint rows and instead preserve constraint values by other means, but these also require additional logic. In general, your implementations of `prolong`, `restrictWithTranspose`, and `restrictWithPinv` will have to account for whatever data is included in your problem.

It's advisable to give your implementation of this class a constructor with no arguments, and a function to later initialize it with specific data such as `J`.

### 3. Implement the multigrid domain

The class `MultigridDomain` defines an interface with several operations relating to the geometric domain (e.g. a square grid, a space curve, a triangle mesh, etc.). You will have to implement a class that extends `MultigridDomain<Mult, Operator>`, where you use the classes that you implemented in the previous 2 steps as the arguments `Mult` and `Operator`. So, for instance, a class declaration might look like:
```c++
class ConstraintProjectorDomain : public MultigridDomain<BlockClusterTree, MatrixProjectorOperator> {
    ...
}
```
Here, `BlockClusterTree` and `MatrixProjectorOperator` are two classes corresponding to the two above steps, and we extend the specific template instantiation using these two classes.

#### Interface

Your class will have seven virtual functions to implement:

```c++
virtual MultigridDomain<Mult, Operator>* Coarsen(Operator* prolongOp) const = 0;
```
This applies one level of coarsening to the current domain and returns the result. For instance, on a space curve, this might create a coarser curve by deleting every other vertex. When this function finishes, the argument `prolongOp` should contain all the data necessary to perform prolongation and restriction between the fine domain and the coarse domain.


```c++
virtual Mult* GetMultiplier() const = 0;
```

This returns the matrix-vector product associated with the current domain; it should just be a one-liner that returns the appropriate field.

```c++
virtual Eigen::MatrixXd GetFullMatrix() const = 0;
```
This should construct and return the full left-hand side matrix `A` (in the system `Ax = b` that you want to solve) from the current geometric domain.


```c++
virtual Eigen::VectorXd DirectSolve(Eigen::VectorXd &b) const = 0;
```
This should get the matrix `A` using `GetFullMatrix()` and then solve `Ax = b` using a direct method.

```c++
virtual int NumVertices() const = 0;
```
This should return the number of vertices (or in general the number of elements) in the current domain.

```c++
virtual int NumRows() const  = 0;
```
This should return the number of rows in the system `Ax = b` on the current domain. This is not necessarily the same as `NumVertices()`. For instance, each vertex might have an associated 3-vector instead of just one value, making a `3V x 3V` system. Or, if there are explicit constraint rows, then these also need to be accounted for.

```c++
virtual Operator* MakeNewOperator() const = 0;
```
This should just construct and return a new `Operator`, which most likely is just `return new Foo()` with the operator class you defined previously.

### 4. Create a hierarchy and solve

All that's left to do is write the code that uses the above components to perform the actual linear solve. This code will look something like:

```c++
using Domain = YourDomain;
using Solver = MultigridHierarchy<Domain>;

// Get the right-hand side of the system.
Eigen::VectorXd b = getYourData( ... );

// Set up your initial domain -- the solver will automatically
// create the rest of the hierarchy.
Domain* domain = new Domain( ... );

// Create a solver using your initial domain.
Solver* multigrid = new Solver(domain);

// Pick a tolerance for the multigrid solve.
double tol = 0.001;

// Solve the system.
Eigen::VectorXd x = solver->solver->template VCycleSolve<Smoother>(b, tol);
```



