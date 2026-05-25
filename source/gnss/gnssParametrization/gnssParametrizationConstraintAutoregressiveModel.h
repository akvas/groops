/***********************************************/
/**
* @file gnssParametrizationConstraintAutoregressiveModel.h
*
* @brief Constrain parameters with an autoregressive model.
* @see GnssParametrization
*
* @author Andreas Kvas
* @date 2026-05-06
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPARAMETRIZATIONCONSTRAINTARMODEL__
#define __GROOPS_GNSSPARAMETRIZATIONCONSTRAINTARMODEL__

// Latex documentation
#ifdef DOCSTRING_GnssParametrization
static const char *docstringGnssParametrizationConstraintAutoregressiveModel = R"(
\subsection{Constraints}\label{gnssParametrizationType:autoregressiveModel}
Add a pseudo observation equation (constraint)
for each selected \configClass{parameters}{parameterSelectorType}
\begin{equation}
  0 = x_i - \sum_{k=1}^{p} b_k x_{i-k}  + \epsilon, \hspace{15pt} \epsilon \sim \mathcal{N}(0, \sigma^2)
\end{equation}.
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
class GnssParametrizationConstraintAutoregressiveModel : public GnssParametrizationBase
{
  std::string                      name;
  ParameterSelectorPtr             parameterSelector;
  AutoregressiveModelSequencePtr   arSequence;
  Gnss                             *gnss;

public:
GnssParametrizationConstraintAutoregressiveModel(Config &config);

  void init(Gnss *gnss, Parallel::CommunicatorPtr /*comm*/) override {this->gnss = gnss;}
  void constraints(const GnssNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const override;
};

/***********************************************/

inline GnssParametrizationConstraintAutoregressiveModel::GnssParametrizationConstraintAutoregressiveModel(Config &config)
{
  try
  {
    readConfig(config, "name",                             name,              Config::OPTIONAL,      "constraint.name",  "");
    readConfig(config, "parameters",                       parameterSelector, Config::MUSTSET,  "",  "parameters to constrain");
    readConfig(config, "autoregressiveModelSequence",      arSequence,        Config::MUSTSET,  "",  "autoregressive model sequence");
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

inline void GnssParametrizationConstraintAutoregressiveModel::constraints(const GnssNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const
{
  try
  {
    if(!isEnabled(normalEquationInfo, name))
      return;

    const std::vector<UInt> indices = parameterSelector->indexVector(normalEquationInfo.parameterNames());

    auto index = arSequence->distributedNormalsBlockIndex(indices.size());
    Single::forEach(index.size(), [&](UInt k)
    {
      UInt i1 = indices.at(index[k].first);
      UInt i2 = indices.at(index[k].second);
      if(i1 != NULLINDEX && i2 != NULLINDEX)
      {
        const UInt idBlock1    = normals.index2block(i1);
        const UInt blockIndex1 = normals.blockIndex(idBlock1);
        const UInt idBlock2    = normals.index2block(i2);
        const UInt blockIndex2 = normals.blockIndex(idBlock2);

        normals.setBlock(idBlock1, idBlock1);
        normals.setBlock(idBlock1, idBlock2);
        normals.setBlock(idBlock2, idBlock2);

        Matrix N_ik = arSequence->distributedNormalsBlock(normals.blockCount(), index[k].first, index[k].second);
        std::cerr<<N_ik.rows()<<"|"<<N_ik.columns()<<std::endl;
        if(normals.isMyRank(idBlock1, idBlock2))
          normals.N(idBlock1, idBlock2)(i1-blockIndex1, i2-blockIndex2) += N_ik(0, 0);

      }
    });

    if(Parallel::isMaster(normalEquationInfo.comm))
        obsCount += indices.size();
    const UInt count = indices.size();

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
