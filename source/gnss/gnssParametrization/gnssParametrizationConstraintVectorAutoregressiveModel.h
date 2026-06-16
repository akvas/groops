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
based on vector autoregressive (AR) models up to order $p$ in the form
\begin{equation}
  0 = \mathbf{x}_i - \sum_{k=1}^{p} \mathbf{\Phi}_k \mathbf{x}_{i-k}  + \epsilon, \hspace{15pt} \epsilon \sim \mathcal{N}(0, \mathbf{\Sigma})
\end{equation}.
The first epochs $p - 1$ epochs are constrained using AR models of order zero to $p$.
This is equivalent to applying a zero constraint with block-Toeplitz covariance matrix to the parameter sequence.
See \configClass{autoregressiveModelSequence}{autoregressiveModelSequenceType} for the detailed theoretical background.
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
  Bool                             relativeToApriori;

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
    readConfig(config, "relativeToApriori",                relativeToApriori,      Config::DEFAULT,  "0", "constrain only dx and not full x=dx+x0");
    if(isCreateSchema(config)) return;

    if(arSequence->dimension() != parameterPerDimension.size())
      throw(Exception("Dimension mismatch between AR models and selected parameters."));
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

    Vector x0 = Vector(normalEquationInfo.parameterCount());
    if(!relativeToApriori)
      x0 = gnss->aprioriParameter(normalEquationInfo);

    const auto paranames = normalEquationInfo.parameterNames();
    std::vector< std::vector<UInt> > indices;
    for(auto parameterPerDim : parameterPerDimension)
      indices.push_back(parameterPerDim.parameterSelector->indexVector(normalEquationInfo.parameterNames()));

    UInt parameterCount = 0;
    for(auto idx : indices)
    {
      if((parameterCount != 0) && (idx.size() != parameterCount))
        throw(Exception("mismatch in selected parameter count"));

      parameterCount = idx.size();
    }
    const UInt ndim = indices.size();

    auto index = arSequence->distributedNormalsBlockIndex(parameterCount);
    Single::forEach(index.size(), [&](UInt k)
    {
      Matrix N_ik = arSequence->distributedNormalsBlock(parameterCount, index[k].first, index[k].second);
      if(N_ik.getType() == Matrix::SYMMETRIC)
        fillSymmetric(N_ik);

      for(UInt dim1 = 0; dim1 < ndim; dim1++)
        for(UInt dim2 = 0; dim2 < ndim; dim2++)
        {
          UInt i1 = indices.at(dim1).at(index[k].first);
          UInt i2 = indices.at(dim2).at(index[k].second);
          if(i1 != NULLINDEX && i2 != NULLINDEX)
          {
            const UInt idBlock1    = normals.index2block(i1);
            const UInt blockIndex1 = normals.blockIndex(idBlock1);
            const UInt idBlock2    = normals.index2block(i2);
            const UInt blockIndex2 = normals.blockIndex(idBlock2);

            normals.setBlock(idBlock1, idBlock1);
            normals.setBlock(idBlock1, idBlock2);
            normals.setBlock(idBlock2, idBlock2);

            if(normals.isMyRank(idBlock1, idBlock2))
              normals.N(idBlock1, idBlock2)(i1-blockIndex1, i2-blockIndex2) += N_ik(dim1, dim2);

            if(Parallel::isMaster(normalEquationInfo.comm))
            {
              n.at(idBlock1)(i1 - blockIndex1, 0) -= N_ik(dim1, dim2) * x0.at(i2);
              lPl += N_ik(dim1, dim2) * x0.at(i1) * x0.at(i2);
              if(index[k].first != index[k].second)  // consider lower triangle
              {
                n.at(idBlock1)(i2 - blockIndex2, 0) -= N_ik(dim2, dim1) * x0.at(i1);
                lPl += N_ik(dim2, dim1) * x0.at(i1) * x0.at(i2);
              }
            }

          }
        }
    }, FALSE/*timing*/);

    const UInt count = parameterCount * ndim;
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
