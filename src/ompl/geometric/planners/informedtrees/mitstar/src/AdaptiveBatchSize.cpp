/*********************************************************************
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2023, Technical University of Munich
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
 *   * Neither the name of the University of Munich nor the names of its
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

// Authors: Liding Zhang

#include "ompl/geometric/planners/informedtrees/mitstar/AdaptiveBatchSize.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

namespace ompl
{
    namespace geometric
    {
        namespace mitstar
        {
            AdaptiveBatchSize::AdaptiveBatchSize(const DecayMethod &decay_method, ompl::base::Cost &solutionCost,
                                                 double informedMeasure, unsigned int &batchSize,
                                                 double &S_max_initial,
                                                 const unsigned int &maxSamples, const unsigned int &minSamples,
                                                 std::size_t dim)
              : decay_method_(decay_method)
              , _solutionCost_(solutionCost)
              , informedMeasure_(informedMeasure)
              , _batchSize_(batchSize)
              , _S_max_initial_(S_max_initial)
              , _maxSamples_(maxSamples)
              , _minSamples_(minSamples)
              , dim_(dim)
            {
            }
            unsigned int AdaptiveBatchSize::adjustBatchSize(DecayMethod decay_method)
            {
                switch (decay_method)
                {
                    case ITERATION_TIMES:
                        return AdaptiveBatchSize::adjustBatchSizeIterationTimes();
                    case LINEAR:
                        return AdaptiveBatchSize::adjustBatchSizeLinear();
                    case PARABOLA:
                        return AdaptiveBatchSize::adjustBatchSizeParabola();
                    case LOG:
                        return AdaptiveBatchSize::adjustBatchSizeLog();
                    case BRACHISTOCHRONE:
                        return AdaptiveBatchSize::adjustBatchSizeBrachistochrone();
                    default:
                        return _batchSize_;
                }
            }

            unsigned int AdaptiveBatchSize::adjustBatchSizeIterationTimes()
            {
                // Implementation for ITERATION_TIMES method
                if (_iteration_ >= 1000 && _iteration_ < 5000)
                {
                    return _batchSize_ = 150;
                }

                else if (_iteration_ >= 5000 && _iteration_ < 8000)
                {
                    return _batchSize_ = 100;
                }
                else if (_iteration_ >= 8000 && _iteration_ < 10000)
                {
                    return _batchSize_ = 80;
                }
                else
                {
                    return _batchSize_;
                }
            }

            unsigned int AdaptiveBatchSize::adjustBatchSizeLinear()
            {
                const double ratio = measureRatio();

                _batchSize_ = _minSamples_ + (_maxSamples_ - _minSamples_) * ratio;

                // Clamp batchSize_ to be within range
                _batchSize_ = std::max(_minSamples_, std::min(_maxSamples_, _batchSize_));
                return _batchSize_;
            }

            unsigned int AdaptiveBatchSize::adjustBatchSizeParabola()
            {
                const double ratio = measureRatio();

                // parabola curve
                double parabola_factor = std::sqrt(ratio);

                _batchSize_ = _minSamples_ + (_maxSamples_ - _minSamples_) * parabola_factor;
                // Clamp batchSize_ to be within [minSamples_, maxSamples_]
                _batchSize_ = std::max(_minSamples_, std::min(_maxSamples_, _batchSize_));
                return _batchSize_;
            }

            unsigned int AdaptiveBatchSize::adjustBatchSizeLog()
            {
                // Implementation for LOG & SIGMOID method
                if (std::isinf(_solutionCost_.value()))
                {
                    return _batchSize_;
                }

                const double ratio = measureRatio();

                // Logarithmic decay
                const double lambda = static_cast<double>(_minSamples_ + _maxSamples_) /
                                      static_cast<double>(dim_);  // Adjust this value to make the decay faster or slower
                // Sigmoid function to smooth ratio
                const double smoothedValue = 1.0 / (1.0 + std::exp(-10.0 * (ratio - 0.5)));
                const double decay_factor =
                    std::log(1.0 + lambda * smoothedValue) / std::log(1.0 + lambda);
                _batchSize_ = _minSamples_ + (_maxSamples_ - _minSamples_) * decay_factor;

                // Clamp batchSize_ to be within [minSamples_, maxSamples_]
                _batchSize_ = std::max(_minSamples_, std::min(_maxSamples_, _batchSize_));
                // areaList_[ratio]++;

                std::cout << "MIT* adaptive batch: cost=" << _solutionCost_.value()
                          << ", informed_measure=" << informedMeasure_
                          << ", initial_measure=" << _S_max_initial_ << ", ratio=" << ratio
                          << ", batch_size=" << _batchSize_ << std::endl;
                return _batchSize_;
            }

            unsigned int AdaptiveBatchSize::adjustBatchSizeBrachistochrone()
            {
                const double ratio = measureRatio();

                // Sigmoid function to smooth ratio
                double smoothedValue = 1 / (1 + exp(-10 * (ratio - 0.5)));
                // Brachistochrone curve decay
                double theta = M_PI / 2 * smoothedValue;
                double brachistochrone_factor = std::sqrt(sin(theta));
                _batchSize_ = _minSamples_ + (_maxSamples_ - _minSamples_) * brachistochrone_factor;

                // Clamp batchSize_ to be within [minSamples_, maxSamples_]
                _batchSize_ = std::max(_minSamples_, std::min(_maxSamples_, _batchSize_));
                return _batchSize_;
            }

            double AdaptiveBatchSize::measureRatio()
            {
                if (!std::isfinite(informedMeasure_) || informedMeasure_ < 0.0)
                {
                    return 1.0;
                }
                if (!std::isfinite(_S_max_initial_) || _S_max_initial_ <= 0.0)
                {
                    _S_max_initial_ = informedMeasure_;
                }
                if (_S_max_initial_ <= 0.0)
                {
                    return 0.0;
                }
                return std::max(
                    0.0,
                    std::min(1.0, informedMeasure_ / _S_max_initial_));
            }

        }  // namespace mitstar

    }  // namespace geometric

}  // namespace ompl
