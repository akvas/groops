/***********************************************/
/**
* @file gnssParametrizationConstraintMultivariateRandomWalk.h
*
* @brief Constrain parameters with a multivariate random walk.
* @see GnssParametrization
*
* @author Andreas Kvas
* @date 2026-10-04
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPARAMETRIZATIONCONSTRAINTMULTIVARIATERANDOMWALK__
#define __GROOPS_GNSSPARAMETRIZATIONCONSTRAINTMULTIVARIATERANDOMWALK__

// Latex documentation
#ifdef DOCSTRING_GnssParametrization
static const char *docstringGnssParametrizationConstraintMultivariateRandomWalk = R"(
\subsection{ConstraintMultivariateRandomWalk}\label{gnssParametrizationType:constraintMultivariateRandomWalk}
Add a pseudo observation equation (constraint) for each consecutive epoch of the selected parameters,
with one parameter selector for each dimension:
\begin{equation}
  0 = \mathbf{x}_{k} - \mathbf{x}_{k-1} + \epsilon, \hspace{15pt} \epsilon \sim \mathcal{N}(0, \mathbf{C}),
\end{equation}
where the configured covariance matrix $\mathbf{C}$ couples the parameter dimensions.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "gnss/gnss.h"
#include "classes/matrixGenerator/matrixGenerator.h"
#include "classes/parameterSelector/parameterSelector.h"
#include "gnss/gnssParametrization/gnssParametrization.h"

/***** CLASS ***********************************/

/** @brief Parameter constraints with a multivariate random walk.
* @ingroup gnssParametrizationGroup
* @see GnssParametrization */
class GnssParametrizationConstraintMultivariateRandomWalk : public GnssParametrizationBase
{
  std::string name;
  Matrix covarianceInverse;
  Gnss *gnss;
  Bool relativeToApriori;

public:
  class ParameterPerDimension
  {
  public:
    ParameterSelectorPtr parameterSelector;
  };

  std::vector<ParameterPerDimension> parameterPerDimension;

  GnssParametrizationConstraintMultivariateRandomWalk(Config &config);

  void init(Gnss *gnss, Parallel::CommunicatorPtr /*comm*/) override {this->gnss = gnss;}
  void constraints(const GnssNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const override;
};

/***********************************************/

template<> Bool readConfig(Config &config, const std::string &name, GnssParametrizationConstraintMultivariateRandomWalk::ParameterPerDimension &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  if(!readConfigSequence(config, name, mustSet, defaultValue, annotation))
    return FALSE;
  readConfig(config, "parameters", var.parameterSelector, Config::MUSTSET, "", "parameters to constrain");
  endSequence(config);
  return TRUE;
}

/***********************************************/

inline GnssParametrizationConstraintMultivariateRandomWalk::GnssParametrizationConstraintMultivariateRandomWalk(Config &config)
{
  try
  {
    MatrixGeneratorPtr covarianceGenerator;
    readConfig(config, "name",                  name,                  Config::OPTIONAL, "constraint.name", "");
    readConfig(config, "covariance",            covarianceGenerator,   Config::MUSTSET,  "", "covariance matrix of the random walk increments");
    readConfig(config, "parameterPerDimension", parameterPerDimension, Config::MUSTSET,  "", "parameter selector for each dimension");
    readConfig(config, "relativeToApriori",     relativeToApriori,     Config::DEFAULT,  "0", "constrain only dx and not full x=dx+x0");
    if(isCreateSchema(config)) return;

    if(parameterPerDimension.empty())
      throw(Exception("At least one parameter dimension must be selected."));

    Matrix covariance = covarianceGenerator->compute();
    if((covariance.getType() != Matrix::SYMMETRIC) ||
       (covariance.rows() != parameterPerDimension.size()) ||
       (covariance.columns() != parameterPerDimension.size()))
      throw(Exception("Covariance matrix must be symmetric and match the number of parameter dimensions."));

    fillSymmetric(covariance);
    inverse(covariance);
    fillSymmetric(covariance);
    covarianceInverse = covariance;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void GnssParametrizationConstraintMultivariateRandomWalk::constraints(const GnssNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const
{
  try
  {
    if(!isEnabled(normalEquationInfo, name))
      return;

    Vector x0 = Vector(normalEquationInfo.parameterCount());
    if(!relativeToApriori)
      x0 = gnss->aprioriParameter(normalEquationInfo);

    std::vector<std::vector<UInt>> indices;
    for(const auto &parameterPerDim : parameterPerDimension)
      indices.push_back(parameterPerDim.parameterSelector->indexVector(normalEquationInfo.parameterNames()));

    const UInt parameterCount = indices.front().size();
    for(const auto &dimensionIndices : indices)
      if(dimensionIndices.size() != parameterCount)
        throw(Exception("mismatch in selected parameter count"));

    const UInt dimension = indices.size();
    UInt count = 0;
    auto addNormal = [&](UInt i, UInt j, Double value)
    {
      if(i > j)
        std::swap(i, j);

      const UInt idBlock1 = normals.index2block(i);
      const UInt idBlock2 = normals.index2block(j);
      const UInt blockIndex1 = normals.blockIndex(idBlock1);
      const UInt blockIndex2 = normals.blockIndex(idBlock2);
      normals.setBlock(idBlock1, idBlock2);
      if(normals.isMyRank(idBlock1, idBlock2))
        normals.N(idBlock1, idBlock2)(i-blockIndex1, j-blockIndex2) += value;
    };
    auto addRightHandSide = [&](UInt index, Double value)
    {
      const UInt idBlock = normals.index2block(index);
      if(Parallel::isMaster(normalEquationInfo.comm))
        n.at(idBlock)(index - normals.blockIndex(idBlock), 0) += value;
    };

    for(UInt k = 0; k + 1 < parameterCount; k++)
    {
      Bool completeIncrement = TRUE;
      for(UInt dim = 0; dim < dimension; dim++)
        if((indices.at(dim).at(k) == NULLINDEX) || (indices.at(dim).at(k + 1) == NULLINDEX))
          completeIncrement = FALSE;
      if(!completeIncrement)
        continue;

      Vector difference(dimension);
      for(UInt dim = 0; dim < dimension; dim++)
        difference(dim) = x0.at(indices.at(dim).at(k + 1)) - x0.at(indices.at(dim).at(k));

      for(UInt dim1 = 0; dim1 < dimension; dim1++)
        for(UInt dim2 = 0; dim2 < dimension; dim2++)
        {
          const Double weight = covarianceInverse(dim1, dim2);
          const UInt previous1 = indices.at(dim1).at(k);
          const UInt current1 = indices.at(dim1).at(k + 1);
          const UInt previous2 = indices.at(dim2).at(k);
          const UInt current2 = indices.at(dim2).at(k + 1);

          addRightHandSide(previous1,  weight * difference(dim2));
          addRightHandSide(current1,  -weight * difference(dim2));
          if(Parallel::isMaster(normalEquationInfo.comm))
            lPl += difference(dim1) * weight * difference(dim2);

          if(dim1 <= dim2)
          {
            addNormal(previous1, previous2, weight);
            addNormal(current1, current2, weight);
          }
          addNormal(previous1, current2, -weight);
        }
      count += dimension;
    }

    if(Parallel::isMaster(normalEquationInfo.comm))
      obsCount += count;

    if(count)
      logStatus<<"constrain "<<name<<" ("<<count<<" parameters)"<<Log::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif