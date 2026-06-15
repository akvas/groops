/***********************************************/
/**
* @file gnssParametrizationConstraintVectorAutoregressiveModel.h
*
* @brief Constrain parameters with a vector autoregressive model.
* @see GnssParametrization
*
* @author Andreas Kvas
* @date 2026-06-12
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPARAMETRIZATIONCONSTRAINTVARMODEL__
#define __GROOPS_GNSSPARAMETRIZATIONCONSTRAINTVARMODEL__

// Latex documentation
#ifdef DOCSTRING_GnssParametrization
static const char *docstringGnssParametrizationConstraintVectorAutoregressiveModel = R"(
\subsection{ConstraintVectorAutoregressiveModel}\label{gnssParametrizationType:vectorAutoregressiveModel}
Add a pseudo observation equation (constraint)
for each selected \configClass{parameters}{parameterSelectorType}
based on autoregressive (AR) models up to order $p$ in the form
\begin{equation}
  0 = x_i - \sum_{k=1}^{p} b_k x_{i-k}  + \epsilon, \hspace{15pt} \epsilon \sim \mathcal{N}(0, \sigma^2)
\end{equation}.
The first epochs $p - 1$ epochs are constrained using AR models of order zero to $p$.
This is equivalent to applying a zero constraint with Toeplitz covariance matrix to the parameter time series.
See \configClass{autoregressiveModelSequence}{autoregressiveModelSequenceType} for the detailed theoretical background.

The standard deviation $\sigma$ (\config{sigma}) is used to weight the observation equations.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "gnss/gnss.h"
#include "misc/kalmanProcessing.h"
#include "classes/parameterSelector/parameterSelector.h"
#include "gnss/gnssParametrization/gnssParametrization.h"

/***** CLASS ***********************************/

/** @brief Parameter constraints with autoregressive model.
* @ingroup gnssParametrizationGroup
* @see GnssParametrization */
class GnssParametrizationConstraintVectorAutoregressiveModel : public GnssParametrizationBase
{
  std::string                      name;
  AutoregressiveModelSequencePtr   arSequence;
  Gnss                             *gnss;

public:
  class ParameterPerDimension
  {
    public:
      ParameterSelectorPtr parameterSelector;
  };

  std::vector<ParameterPerDimension> parameterPerDimension;

  GnssParametrizationConstraintVectorAutoregressiveModel(Config &config);

  void init(Gnss *gnss, Parallel::CommunicatorPtr /*comm*/) override {this->gnss = gnss;}
  void constraints(const GnssNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const override;
};

/***********************************************/

template<> Bool readConfig(Config &config, const std::string &name, GnssParametrizationConstraintVectorAutoregressiveModel::ParameterPerDimension &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  if(!readConfigSequence(config, name, mustSet, defaultValue, annotation))
    return FALSE;
  readConfig(config, "parameters", var.parameterSelector, Config::MUSTSET,  "",  "parameters to constrain");
  endSequence(config);
  return TRUE;
}


/***********************************************/

inline GnssParametrizationConstraintVectorAutoregressiveModel::GnssParametrizationConstraintVectorAutoregressiveModel(Config &config)
{
  try
  {
    readConfig(config, "name",                             name,                   Config::OPTIONAL,      "constraint.name",  "");
    readConfig(config, "autoregressiveModelSequence",      arSequence,             Config::MUSTSET,  "",  "autoregressive model sequence");
    readConfig(config, "parameterPerDimension",            parameterPerDimension,  Config::MUSTSET,  "",  "parameter selector for each dimension");
    if(isCreateSchema(config)) return;

    if(arSequence->dimension() != 1)
      throw("Only one dimensional AR models are supported.");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void GnssParametrizationConstraintVectorAutoregressiveModel::constraints(const GnssNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const
{
  try
  {
    if(!isEnabled(normalEquationInfo, name))
      return;

    // const std::vector<UInt> indices = parameterSelector->indexVector(normalEquationInfo.parameterNames());

    // auto index = arSequence->distributedNormalsBlockIndex(indices.size());
    // Single::forEach(index.size(), [&](UInt k)
    // {
    //   UInt i1 = indices.at(index[k].first);
    //   UInt i2 = indices.at(index[k].second);
    //   if(i1 != NULLINDEX && i2 != NULLINDEX)
    //   {
    //     const UInt idBlock1    = normals.index2block(i1);
    //     const UInt blockIndex1 = normals.blockIndex(idBlock1);
    //     const UInt idBlock2    = normals.index2block(i2);
    //     const UInt blockIndex2 = normals.blockIndex(idBlock2);

    //     normals.setBlock(idBlock1, idBlock1);
    //     normals.setBlock(idBlock1, idBlock2);
    //     normals.setBlock(idBlock2, idBlock2);

    //     Matrix N_ik = arSequence->distributedNormalsBlock(normals.blockCount(), index[k].first, index[k].second);
    //     if(normals.isMyRank(idBlock1, idBlock2))
    //       normals.N(idBlock1, idBlock2)(i1-blockIndex1, i2-blockIndex2) += N_ik(0, 0);

    //   }
    // });

    // if(Parallel::isMaster(normalEquationInfo.comm))
    //     obsCount += indices.size();
    // const UInt count = indices.size();

    // if(count)
    //   logStatus<<"constrain "<<name<<" ("<<count<<" parameters)"<<Log::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
