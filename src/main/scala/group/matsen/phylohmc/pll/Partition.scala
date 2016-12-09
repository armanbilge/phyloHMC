package group.matsen.phylohmc.pll

import com.sun.jna.Native

class Partition(tips: Int, clvBuffers: Int, states: Int, sites: Int, rateMatrices: Int, probMatrices: Int, rateCats: Int, scaleBuffers: Int, sse: Boolean = false, avx: Boolean = false, avx2: Boolean = false, avx512: Boolean = false, tipPatternCompression: Boolean = false) {

  private[this] val attributes = {
    var attributes = 0
    if (sse) attributes |= PLL_ATTRIB_ARCH_SSE
    if (avx) attributes |= PLL_ATTRIB_ARCH_AVX
    if (avx2) attributes |= PLL_ATTRIB_ARCH_AVX2
    if (avx512) attributes |= PLL_ATTRIB_ARCH_AVX512
    if (tipPatternCompression) attributes |= PLL_ATTRIB_PATTERN_TIP
    attributes
  }

  private[this] val self = pll_partition_create(tips, clvBuffers, states, sites, rateMatrices, probMatrices, rateCats, scaleBuffers, attributes)
  self.read()

  private[this] val sumtable = pll_aligned_alloc(new size_t(self.sites * self.rate_cats * self.states_padded * Native.getNativeSize(classOf[Double])), if (sse) PLL_ALIGNMENT_SSE else if (avx || avx2) PLL_ALIGNMENT_AVX else PLL_ALIGNMENT_CPU)

  override def finalize(): Unit = {
    pll_aligned_free(sumtable)
    pll_partition_destroy(self)
  }

  def setTipStates(tipIndex: Int, dataType: DataType, sequence: String) = pll_set_tip_states(self, tipIndex, dataType match {
    case Bin => pll_map_bin
    case Nt => pll_map_nt
    case AA => pll_map_aa
  }, sequence)

  def setTipCLV(tipIndex: Int, clv: Array[Double]) = pll_set_tip_clv(self, tipIndex, clv)

  def setPatternWeights(patternWeights: Array[Int]) = pll_set_pattern_weights(self, patternWeights)

  def setSubstParams(paramsIndex: Int, params: Array[Double]) = pll_set_subst_params(self, paramsIndex, params)

  def setFrequencies(paramsIndex: Int, params: Array[Double]) = pll_set_frequencies(self, paramsIndex, params)

  def setCategoryRates(rates: Array[Double]) = pll_set_category_rates(self, rates)

  def setCategoryWeights(rateWeights: Array[Double]) = pll_set_category_weights(self, rateWeights)

  def updateEigen(paramsIndex: Int) = pll_update_eigen(self, paramsIndex)

  def updateProbMatrices(paramsIndices: Array[Int], matrixIndices: Array[Int], branchLengths: Array[Double], count: Int) = pll_update_prob_matrices(self, paramsIndices, matrixIndices, branchLengths, count)

  def updateInvariantSites() = pll_update_invariant_sites(self)

  def updateInvariantSitesProportion(paramsIndex: Int, propInvar: Double) = pll_update_invariant_sites_proportion(self, paramsIndex, propInvar)

  def updatePartials(operations: Operations) = pll_update_partials(self, operations.pll_operation, operations.count)

  def computeRootLogLikelihood(clvIndex: Int, scalerIndex: Int, freqsIndex: Array[Int], persiteLnL: Array[Double]) = pll_compute_root_loglikelihood(self, clvIndex, scalerIndex, freqsIndex, persiteLnL)

  def computeEdgeLogLikelihood(parentCLVIndex: Int, parentScalerIndex: Int, childCLVIndex: Int, childScalerIndex: Int, matrixIndex: Int, freqsIndex: Array[Int], persiteLnL: Array[Double]) = pll_compute_edge_loglikelihood(self, parentCLVIndex, parentScalerIndex, childCLVIndex, childScalerIndex, matrixIndex, freqsIndex, persiteLnL)

  def updateSumtable(parentCLVIndex: Int, childCLVIndex: Int, paramsIndices: Array[Int]) = pll_update_sumtable(self, parentCLVIndex, childCLVIndex, paramsIndices, sumtable)

  def computeLikelihoodDerivatives(parentScalerIndex: Int, childScalerIndex: Int, branchLength: Double, paramsIndices: Array[Int], df: Array[Double], ddf: Array[Double]) = pll_compute_likelihood_derivatives(self, parentScalerIndex, childScalerIndex, branchLength, paramsIndices, sumtable, df, ddf)

}
