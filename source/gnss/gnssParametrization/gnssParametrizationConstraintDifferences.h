/***********************************************/
/**
* @file gnssParametrizationConstraintDifferences.h
*
* @brief Constraint on parameter differences.
* @see GnssParametrization
*
* @author Andreas Kvas
* @date 2025-05-05
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPARAMETRIZATIONCONSTRAINTDIFFERENCES__
#define __GROOPS_GNSSPARAMETRIZATIONCONSTRAINTDIFFERENCES__

// Latex documentation
#ifdef DOCSTRING_GnssParametrization
static const char *docstringGnssParametrizationConstraintDifferences = R"(
\subsection{ConstraintDifferences}\label{gnssParametrizationType:constraintDifferences}
Add a pseudo observation equation (constraint)
for each selected \configClass{parameters}{parameterSelectorType} in the form
\begin{equation}
  0 = -1 \cdot \Delta x_{k-1} + 1 \cdot \Delta x_k  + \epsilon,
\end{equation}.
This constrains the differences of consecutive parameters, with \config{sigma} used to weight the observation equations.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "gnss/gnss.h"
#include "classes/parameterSelector/parameterSelector.h"
#include "gnss/gnssParametrization/gnssParametrization.h"

/***** CLASS ***********************************/

/** @brief Parameter constraints.
* @ingroup gnssParametrizationGroup
* @see GnssParametrization */
class GnssParametrizationConstraintDifferences : public GnssParametrizationBase
{
  std::string          name;
  ParameterSelectorPtr parameterSelector;
  Double               sigma;
  Bool                 relativeToApriori;
  Gnss                *gnss;

public:
GnssParametrizationConstraintDifferences(Config &config);

  void init(Gnss *gnss, Parallel::CommunicatorPtr /*comm*/) override {this->gnss = gnss;}
  void constraints(const GnssNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const override;
};

/***********************************************/

inline GnssParametrizationConstraintDifferences::GnssParametrizationConstraintDifferences(Config &config)
{
  try
  {
    readConfig(config, "name",              name,              Config::OPTIONAL, "constraint.name",  "");
    readConfig(config, "parameters",        parameterSelector, Config::MUSTSET,  "",  "parameters to constrain");
    readConfig(config, "sigma",             sigma,             Config::MUSTSET,  "",  "sigma of the constraint (same unit as parameter)");
    readConfig(config, "relativeToApriori", relativeToApriori, Config::DEFAULT,  "0", "constrain only dx and not full x=dx+x0");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void GnssParametrizationConstraintDifferences::constraints(const GnssNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const
{
  try
  {
    if(!isEnabled(normalEquationInfo, name))
      return;

    Vector x0 = Vector(normalEquationInfo.parameterCount());
    if(!relativeToApriori)
      x0 = gnss->aprioriParameter(normalEquationInfo);

    const Double weight = 1./std::pow(sigma, 2);
    const std::vector<UInt> indices = parameterSelector->indexVector(normalEquationInfo.parameterNames());
    UInt count = 0;
    for(UInt k = 0; k + 1 < indices.size(); k++)
    {
      UInt i1 = indices.at(k);
      UInt i2 = indices.at(k + 1);
      if(i1 != NULLINDEX && i2 != NULLINDEX)
      {
        const UInt idBlock1    = normals.index2block(i1);
        const UInt blockIndex1 = normals.blockIndex(idBlock1);
        const UInt idBlock2    = normals.index2block(i2);
        const UInt blockIndex2 = normals.blockIndex(idBlock2);

        normals.setBlock(idBlock1, idBlock1);
        normals.setBlock(idBlock1, idBlock2);
        normals.setBlock(idBlock2, idBlock2);

        if(normals.isMyRank(idBlock1, idBlock1))
          normals.N(idBlock1, idBlock1)(i1-blockIndex1, i1-blockIndex1) += weight;
        if(normals.isMyRank(idBlock1, idBlock2))
          normals.N(idBlock1, idBlock2)(i1-blockIndex1, i2-blockIndex2) -= weight;
        if(normals.isMyRank(idBlock2, idBlock2))
          normals.N(idBlock2, idBlock2)(i2-blockIndex2, i2-blockIndex2) += weight;

        if(Parallel::isMaster(normalEquationInfo.comm))
        {
          n.at(idBlock1)(i1 - blockIndex1, 0) -= weight * (x0.at(i1) - x0.at(i2));
          n.at(idBlock2)(i2 - blockIndex2, 0) += weight * (x0.at(i1) - x0.at(i2));
          lPl += weight * std::pow(x0.at(i1) - x0.at(i2), 2);
          obsCount++;
        }
        count++;
      }
    }

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
