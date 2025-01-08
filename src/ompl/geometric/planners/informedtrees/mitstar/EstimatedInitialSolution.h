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

// Authors: Liding Zhang, Ziqian Xu

#ifndef OMPL_GEOMETRIC_PLANNERS_INFORMEDTREES_MITSTAR_ESTIMATEDINITIALSOLUTION_
#define OMPL_GEOMETRIC_PLANNERS_INFORMEDTREES_MITSTAR_ESTIMATEDINITIALSOLUTION_

#include "ompl/geometric/planners/informedtrees/mitstar/RandomGeometricGraph.h"
#include "ompl/geometric/planners/informedtrees/MITstar.h"
#include "ompl/geometric/planners/informedtrees/mitstar/Edge.h"


namespace ompl
{
    namespace geometric
    {
        namespace mitstar
        {
            class EstimatedInitialSolution
            {
            public:
                EstimatedInitialSolution();
                ~EstimatedInitialSolution() = default;

                /** \brief Get the estimated heuristic initial solution */
                ompl::base::Cost EstimatedHeuristicInitialSolution(const mitstar::Edge &edge);

                /** \brief If we cannot find a true solution cost regarding current estimated heuristic initial solution, 
                 * we need to adjust it according to the expand factor */
                ompl::base::Cost adjustedEstimatedHeuristicInitialSolution(ompl::base::Cost EstimatedCost_, double expandfactor);

                /** \brief The expand factor of estimated heuristic initial solution*/
                double expand_factor = 0.0;

            private:

                /** \brief The part of the true solution cost in the forward tree */
                ompl::base::Cost trueCostPart{std::numeric_limits<double>::infinity()};

                /** \brief The part of the estimated heuristic initial solution cost in the reverse tree */
                ompl::base::Cost admissibleCostPart{std::numeric_limits<double>::infinity()};

                /** \brief The estimated heuristic initial cost */
                ompl::base::Cost EstimatedHeuristicInitialCost{std::numeric_limits<double>::infinity()};

                /** \brief The reliability of the estimated heuristic initial solution cost */
                double reliability = 0.0;
            };

        }  // namespace mitstar

    }  // namespace geometric

}  // namespace ompl

#endif  // OMPL_GEOMETRIC_PLANNERS_INFORMEDTREES_MITSTAR_ESTIMATEDINITIALSOLUTION_