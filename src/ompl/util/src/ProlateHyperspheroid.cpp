/*********************************************************************
* Software License Agreement (BSD License)
*
*  Copyright (c) 2014, University of Toronto
*  All rights reserved.
*
*  Redistribution and use in source and binary forms, with or without
*  modification, are permitted provided that the following conditions
*  are met:
*
*   * Redistributions of source code must retain the above copyright
*     notice, this list of conditions and the following disclaimer.
*   * Redistributions in binary form must reproduce the above
*     copyright notice, this list of conditions and the following
*     disclaimer in the documentation and/or other materials provided
*     with the distribution.
*   * Neither the name of the University of Toronto nor the names of its
*     contributors may be used to endorse or promote products derived
*     from this software without specific prior written permission.
*
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
*  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
*  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
*  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
*  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
*  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
*  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
*  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
*  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
*  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
*  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
*  POSSIBILITY OF SUCH DAMAGE.
*********************************************************************/

/* Author: Jonathan Gammell*/

// The class's header
#include "ompl/util/ProlateHyperspheroid.h"
// For OMPL exceptions
#include "ompl/util/Exception.h"
// For OMPL information
#include "ompl/util/Console.h"
// For geometric equations like prolateHyperspheroidMeasure
#include "ompl/util/GeometricEquations.h"

// For std::make_shared
#include <algorithm>
#include <cmath>
#include <memory>

// Eigen core:
#include <eigen3/Eigen/Core>
// Inversion and determinants
#include <Eigen/LU>
#include <Eigen/Eigenvalues>
// SVD decomposition
#include <Eigen/SVD>

#include <boost/math/special_functions/beta.hpp>

struct ompl::ProlateHyperspheroid::PhsData
{
    /** \brief The dimension of the prolate hyperspheroid.*/
    unsigned int dim_;
    /** \brief Whether the transformation is up to date */
    bool isTransformUpToDate_;
    /** \brief The minimum possible transverse diameter of the PHS. Defined as the distance between the two foci*/
    double minTransverseDiameter_;
    /** \brief The transverse diameter of the PHS. */
    double transverseDiameter_;
    /** \brief The measure of the PHS. */
    double phsMeasure_;
    /** \brief The first focus of the PHS (i.e., the start state of the planning problem). Unlike other parts of Eigen,
     * variably sized matrices do not require special allocators. */
    Eigen::VectorXd xFocus1_;
    /** \brief The second focus of the PHS (i.e., the goal state of the planning problem). Unlike other parts of Eigen,
     * variably sized matrices do not require special allocators.  */
    Eigen::VectorXd xFocus2_;
    /** \brief The centre of the PHS. Defined as the average of the foci. Unlike other parts of Eigen, variably sized
     * matrices do not require special allocators. */
    Eigen::VectorXd xCentre_;
    /** \brief The rotation from PHS-frame to world frame. Is only calculated on construction. Unlike other parts of
     * Eigen, variably sized matrices do not require special allocators. */
    Eigen::MatrixXd rotationWorldFromEllipse_;
    /** \brief The transformation from PHS-frame to world frame. Is calculated every time the transverse diameter
     * changes. Unlike other parts of Eigen, variably sized matrices do not require special allocators. */
    Eigen::MatrixXd transformationWorldFromEllipse_;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

ompl::ProlateHyperspheroid::ProlateHyperspheroid(unsigned int n, const double focus1[], const double focus2[])
  : dataPtr_(std::make_shared<PhsData>())
{
    // Initialize the data:
    dataPtr_->dim_ = n;
    dataPtr_->transverseDiameter_ = 0.0;  // Initialize to something.
    dataPtr_->isTransformUpToDate_ = false;

    // Copy the arrays into their Eigen containers via the Eigen::Map "view"
    dataPtr_->xFocus1_ = Eigen::Map<const Eigen::VectorXd>(focus1, dataPtr_->dim_);
    dataPtr_->xFocus2_ = Eigen::Map<const Eigen::VectorXd>(focus2, dataPtr_->dim_);

    // Calculate the minimum transverse diameter
    dataPtr_->minTransverseDiameter_ = (dataPtr_->xFocus1_ - dataPtr_->xFocus2_).norm();

    // Calculate the centre:
    dataPtr_->xCentre_ = 0.5 * (dataPtr_->xFocus1_ + dataPtr_->xFocus2_);

    // Calculate the rotation
    updateRotation();
}

void ompl::ProlateHyperspheroid::setTransverseDiameter(double transverseDiameter)
{
    if (transverseDiameter < dataPtr_->minTransverseDiameter_)
    {
        OMPL_ERROR("%g < %g", transverseDiameter, dataPtr_->minTransverseDiameter_);
        throw Exception("Transverse diameter cannot be less than the distance between the foci.");
    }

    // Store and update if changed
    if (dataPtr_->transverseDiameter_ != transverseDiameter)
    {
        // Mark as out of date
        dataPtr_->isTransformUpToDate_ = false;

        // Store
        dataPtr_->transverseDiameter_ = transverseDiameter;

        // Update the transform
        updateTransformation();
    }
    // No else, the diameter didn't change
}

void ompl::ProlateHyperspheroid::transform(const double sphere[], double phs[]) const
{
    if (!dataPtr_->isTransformUpToDate_)
    {
        throw Exception("The transformation is not up to date in the PHS class. Has the transverse diameter been set?");
    }

    // Calculate the tranformation and offset, using Eigen::Map views of the data
    Eigen::Map<Eigen::VectorXd>(phs, dataPtr_->dim_) =
        dataPtr_->transformationWorldFromEllipse_ * Eigen::Map<const Eigen::VectorXd>(sphere, dataPtr_->dim_);
    Eigen::Map<Eigen::VectorXd>(phs, dataPtr_->dim_) += dataPtr_->xCentre_;
}

bool ompl::ProlateHyperspheroid::isInPhs(const double point[]) const
{
    if (!dataPtr_->isTransformUpToDate_)
    {
        // The transform is not up to date until the transverse diameter has been set
        throw Exception("The transverse diameter has not been set");
    }

    return (getPathLength(point) < dataPtr_->transverseDiameter_);
}

bool ompl::ProlateHyperspheroid::isOnPhs(const double point[]) const
{
    if (!dataPtr_->isTransformUpToDate_)
    {
        // The transform is not up to date until the transverse diameter has been set
        throw Exception("The transverse diameter has not been set");
    }

    return (getPathLength(point) == dataPtr_->transverseDiameter_);
}

unsigned int ompl::ProlateHyperspheroid::getPhsDimension() const
{
    return dataPtr_->dim_;
}

double ompl::ProlateHyperspheroid::getPhsMeasure() const
{
    if (!dataPtr_->isTransformUpToDate_)
    {
        // The transform is not up to date until the transverse diameter has been set, therefore we have no transverse
        // diameter and we have infinite measure
        return std::numeric_limits<double>::infinity();
    }

    // Calculate and return:
    return dataPtr_->phsMeasure_;
}

double ompl::ProlateHyperspheroid::getPhsMeasure(double tranDiam) const
{
    return prolateHyperspheroidMeasure(dataPtr_->dim_, dataPtr_->minTransverseDiameter_, tranDiam);
}

void ompl::ProlateHyperspheroid::transformSlice(const double sphere[], unsigned int fixedDimension,
                                                double fixedValue, double phs[]) const
{
    if (!dataPtr_->isTransformUpToDate_)
        throw Exception("The transformation is not up to date in the PHS class. Has the transverse diameter been set?");
    if (fixedDimension >= dataPtr_->dim_)
        throw Exception("The fixed PHS slice dimension is out of range.");

    const Eigen::MatrixXd shape = dataPtr_->transformationWorldFromEllipse_ *
                                  dataPtr_->transformationWorldFromEllipse_.transpose();
    const double fixedVariance = shape(fixedDimension, fixedDimension);
    const double fixedOffset = fixedValue - dataPtr_->xCentre_[fixedDimension];
    const double normalizedOffsetSquared = fixedOffset * fixedOffset / fixedVariance;
    if (normalizedOffsetSquared >= 1.0)
        throw Exception("The requested fixed-coordinate slice does not intersect the PHS.");

    const unsigned int sliceDimension = dataPtr_->dim_ - 1u;
    Eigen::VectorXd centre(sliceDimension);
    Eigen::MatrixXd conditionalShape(sliceDimension, sliceDimension);
    unsigned int row = 0u;
    for (unsigned int sourceRow = 0u; sourceRow < dataPtr_->dim_; ++sourceRow)
    {
        if (sourceRow == fixedDimension)
            continue;
        centre[row] = dataPtr_->xCentre_[sourceRow] +
                      shape(sourceRow, fixedDimension) * fixedOffset / fixedVariance;
        unsigned int column = 0u;
        for (unsigned int sourceColumn = 0u; sourceColumn < dataPtr_->dim_; ++sourceColumn)
        {
            if (sourceColumn == fixedDimension)
                continue;
            conditionalShape(row, column) =
                (1.0 - normalizedOffsetSquared) *
                (shape(sourceRow, sourceColumn) - shape(sourceRow, fixedDimension) *
                                                     shape(fixedDimension, sourceColumn) / fixedVariance);
            ++column;
        }
        ++row;
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(conditionalShape);
    if (solver.info() != Eigen::Success)
        throw Exception("Failed to factor the fixed-coordinate PHS slice.");
    const Eigen::VectorXd radii = solver.eigenvalues().cwiseMax(0.0).cwiseSqrt();
    const Eigen::VectorXd slice = centre + solver.eigenvectors() * radii.asDiagonal() *
                                              Eigen::Map<const Eigen::VectorXd>(sphere, sliceDimension);

    unsigned int sliceIndex = 0u;
    for (unsigned int dimension = 0u; dimension < dataPtr_->dim_; ++dimension)
        phs[dimension] = dimension == fixedDimension ? fixedValue : slice[sliceIndex++];
}

double ompl::ProlateHyperspheroid::getPhsMeasure(double tranDiam, unsigned int boundedDimension, double lowerBound,
                                                 double upperBound, int fixedDimension, double fixedValue) const
{
    const unsigned int dimension = dataPtr_->dim_;
    const double transverseRadius = 0.5 * tranDiam;
    const double conjugateRadius =
        0.5 * std::sqrt(tranDiam * tranDiam -
                        dataPtr_->minTransverseDiameter_ * dataPtr_->minTransverseDiameter_);
    Eigen::VectorXd radii = Eigen::VectorXd::Constant(dimension, conjugateRadius);
    radii[0] = transverseRadius;
    const Eigen::MatrixXd shape = dataPtr_->rotationWorldFromEllipse_ * radii.array().square().matrix().asDiagonal() *
                                  dataPtr_->rotationWorldFromEllipse_.transpose();

    auto coordinateFraction = [](double centre, double radius, double lower, double upper, unsigned int dim) {
        auto cdf = [dim](double value) {
            const double clamped = std::max(-1.0, std::min(1.0, value));
            if (clamped == 0.0)
                return 0.5;
            const double beta = boost::math::ibeta(0.5, 0.5 * (static_cast<double>(dim) + 1.0),
                                                   clamped * clamped);
            return 0.5 + std::copysign(0.5 * beta, clamped);
        };
        if (radius <= 0.0 || upper < centre - radius || lower > centre + radius)
            return 0.0;
        return std::max(0.0, cdf((upper - centre) / radius) - cdf((lower - centre) / radius));
    };

    if (fixedDimension < 0)
    {
        const double radius = std::sqrt(shape(boundedDimension, boundedDimension));
        return getPhsMeasure(tranDiam) *
               coordinateFraction(dataPtr_->xCentre_[boundedDimension], radius, lowerBound, upperBound, dimension);
    }

    const unsigned int fixed = static_cast<unsigned int>(fixedDimension);
    const double fixedVariance = shape(fixed, fixed);
    const double fixedOffset = fixedValue - dataPtr_->xCentre_[fixed];
    const double normalizedOffsetSquared = fixedOffset * fixedOffset / fixedVariance;
    if (normalizedOffsetSquared >= 1.0)
        return 0.0;

    Eigen::VectorXd conditionalCentre(dimension - 1u);
    Eigen::MatrixXd conditionalShape(dimension - 1u, dimension - 1u);
    unsigned int row = 0u;
    unsigned int boundedSliceDimension = 0u;
    for (unsigned int sourceRow = 0u; sourceRow < dimension; ++sourceRow)
    {
        if (sourceRow == fixed)
            continue;
        if (sourceRow == boundedDimension)
            boundedSliceDimension = row;
        conditionalCentre[row] = dataPtr_->xCentre_[sourceRow] +
                                 shape(sourceRow, fixed) * fixedOffset / fixedVariance;
        unsigned int column = 0u;
        for (unsigned int sourceColumn = 0u; sourceColumn < dimension; ++sourceColumn)
        {
            if (sourceColumn == fixed)
                continue;
            conditionalShape(row, column) =
                (1.0 - normalizedOffsetSquared) *
                (shape(sourceRow, sourceColumn) -
                 shape(sourceRow, fixed) * shape(fixed, sourceColumn) / fixedVariance);
            ++column;
        }
        ++row;
    }

    const unsigned int sliceDimension = dimension - 1u;
    const double sliceMeasure = unitNBallMeasure(sliceDimension) *
                                std::sqrt(std::max(0.0, conditionalShape.determinant()));
    const double boundedRadius =
        std::sqrt(std::max(0.0, conditionalShape(boundedSliceDimension, boundedSliceDimension)));
    return sliceMeasure * coordinateFraction(conditionalCentre[boundedSliceDimension], boundedRadius, lowerBound,
                                             upperBound, sliceDimension);
}

double ompl::ProlateHyperspheroid::getMinTransverseDiameter() const
{
    return dataPtr_->minTransverseDiameter_;
}

double ompl::ProlateHyperspheroid::getPathLength(const double point[]) const
{
    return (dataPtr_->xFocus1_ - Eigen::Map<const Eigen::VectorXd>(point, dataPtr_->dim_)).norm() +
           (Eigen::Map<const Eigen::VectorXd>(point, dataPtr_->dim_) - dataPtr_->xFocus2_).norm();
}

unsigned int ompl::ProlateHyperspheroid::getDimension() const
{
    return dataPtr_->dim_;
}

void ompl::ProlateHyperspheroid::updateRotation()
{
    // Mark the transform as out of date
    dataPtr_->isTransformUpToDate_ = false;

    // If the minTransverseDiameter_ is too close to 0, we treat this as a circle.
    double circleTol = 1E-9;
    if (dataPtr_->minTransverseDiameter_ < circleTol)
    {
        dataPtr_->rotationWorldFromEllipse_.setIdentity(dataPtr_->dim_, dataPtr_->dim_);
    }
    else
    {
        // Variables
        // The transverse axis of the PHS expressed in the world frame.
        Eigen::VectorXd transverseAxis(dataPtr_->dim_);
        // The matrix representation of the Wahba problem
        Eigen::MatrixXd wahbaProb(dataPtr_->dim_, dataPtr_->dim_);
        // The middle diagonal matrix in the SVD solution to the Wahba problem
        Eigen::VectorXd middleM(dataPtr_->dim_);

        // Calculate the major axis, storing as the first eigenvector
        transverseAxis = (dataPtr_->xFocus2_ - dataPtr_->xFocus1_) / dataPtr_->minTransverseDiameter_;

        // Calculate the rotation that will allow us to generate the remaining eigenvectors
        // Formulate as a Wahba problem, first forming the matrix a_j*a_i' where a_j is the transverse axis if the
        // ellipse in the world frame, and a_i is the first basis vector of the world frame (i.e., [1 0 .... 0])
        wahbaProb = transverseAxis * Eigen::MatrixXd::Identity(dataPtr_->dim_, dataPtr_->dim_).col(0).transpose();

        // Then run it through the  SVD solver
        Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::NoQRPreconditioner> svd(wahbaProb,
                                                                         Eigen::ComputeFullV | Eigen::ComputeFullU);

        // Then calculate the rotation matrix from the U and V components of SVD
        // Calculate the middle diagonal matrix
        middleM = Eigen::VectorXd::Ones(dataPtr_->dim_);
        // Make the last value equal to det(U)*det(V) (zero-based indexing remember)
        middleM(dataPtr_->dim_ - 1) = svd.matrixU().determinant() * svd.matrixV().determinant();

        // Calculate the rotation
        dataPtr_->rotationWorldFromEllipse_ = svd.matrixU() * middleM.asDiagonal() * svd.matrixV().transpose();
    }
}

void ompl::ProlateHyperspheroid::updateTransformation()
{
    // Variables
    // The radii of the ellipse
    Eigen::VectorXd diagAsVector(dataPtr_->dim_);
    // The conjugate diameters:
    double conjugateDiamater;

    // Calculate the conjugate radius
    conjugateDiamater = std::sqrt(dataPtr_->transverseDiameter_ * dataPtr_->transverseDiameter_ -
                                  dataPtr_->minTransverseDiameter_ * dataPtr_->minTransverseDiameter_);

    // Store into the diagonal matrix
    // All the elements but one are the conjugate radius
    diagAsVector.fill(conjugateDiamater / 2.0);

    // The first element in diagonal is the transverse radius
    diagAsVector(0) = 0.5 * dataPtr_->transverseDiameter_;

    // Calculate the transformation matrix
    dataPtr_->transformationWorldFromEllipse_ = dataPtr_->rotationWorldFromEllipse_ * diagAsVector.asDiagonal();

    // Calculate the measure:
    dataPtr_->phsMeasure_ =
        prolateHyperspheroidMeasure(dataPtr_->dim_, dataPtr_->minTransverseDiameter_, dataPtr_->transverseDiameter_);

    // Mark as up to date
    dataPtr_->isTransformUpToDate_ = true;
}
