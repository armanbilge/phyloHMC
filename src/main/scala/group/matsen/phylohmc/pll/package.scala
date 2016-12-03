package group.matsen.phylohmc

import com.sun.jna.{Native, NativeLibrary, Pointer}

package object pll {

  sealed trait DataType
  object Bin extends DataType
  object Nt extends DataType
  object AA extends DataType

  def computeGammaCats(alpha: Double, categories: Int): Array[Double] = {
    val rates = Array.ofDim[Double](categories)
    pll_compute_gamma_cats(alpha, categories, rates)
    rates
  }

  Native.register("pll")

  private[pll] val PLL_ATTRIB_ARCH_SSE = 1 << 0
  private[pll] val PLL_ATTRIB_ARCH_AVX = 1 << 1
  private[pll] val PLL_ATTRIB_ARCH_AVX2 = 1 << 2
  private[pll] val PLL_ATTRIB_ARCH_AVX512 = 1 << 3
  private[pll] val PLL_ATTRIB_PATTERN_TIP = 1 << 4

  private[pll] val PLL_ALIGNMENT_CPU = new size_t(8)
  private[pll] val PLL_ALIGNMENT_SSE = new size_t(16)
  private[pll] val PLL_ALIGNMENT_AVX = new size_t(32)

  private[pll] val pll_map_bin = NativeLibrary.getInstance("pll").getGlobalVariableAddress("pll_map_bin")

  private[pll] val pll_map_nt = NativeLibrary.getInstance("pll").getGlobalVariableAddress("pll_map_nt")

  private[pll] val pll_map_aa = NativeLibrary.getInstance("pll").getGlobalVariableAddress("pll_map_aa")

  @native private[pll] def pll_partition_create(tips: Int, clv_buffers: Int, states: Int, sites: Int, rate_matrices: Int, prob_matrices: Int, rate_cats: Int, scale_buffers: Int, attributes: Int): pll_partition_t

  @native private[pll] def pll_partition_destroy(partition: pll_partition_t): Unit

  @native private[pll] def pll_set_tip_states(partition: pll_partition_t, tip_index: Int, map: Pointer, sequence: String): Int

  @native private[pll] def pll_set_tip_clv(partition: pll_partition_t, tip_index: Int, clv: Array[Double]): Unit

  @native private[pll] def pll_set_pattern_weights(partition: pll_partition_t, pattern_weights: Array[Int]): Unit

  @native private[pll] def pll_set_subst_params(partition: pll_partition_t, params_index: Int, params: Array[Double]): Unit

  @native private[pll] def pll_set_frequencies(partition: pll_partition_t, params_index: Int, frequencies: Array[Double]): Unit

  @native private[pll] def pll_set_category_rates(partition: pll_partition_t, rates: Array[Double]): Unit

  @native private[pll] def pll_set_category_weights(partition: pll_partition_t, rate_weights: Array[Double]): Unit

  @native private[pll] def pll_update_eigen(partition: pll_partition_t, params_index: Int): Unit

  @native private[pll] def pll_update_prob_matrices(partition: pll_partition_t, params_indices: Array[Int], matrix_indices: Array[Int], branch_lengths: Array[Double], count: Int): Unit

  @native private[pll] def pll_update_invariant_sites(partition: pll_partition_t): Int

  @native private[pll] def pll_update_invariant_sites_proportion(partition: pll_partition_t, params_index: Int, prop_invar: Double): Int

  @native private[pll] def pll_aligned_alloc(size: size_t, alignment: size_t): Pointer

  @native private[pll] def pll_aligned_free(ptr: Pointer): Unit

  @native private[pll] def pll_update_partials(partition: pll_partition_t, operations: pll_operation_t, count: Int): Double

  @native private[pll] def pll_compute_root_loglikelihood(partition: pll_partition_t, clv_index: Int, scaler_index: Int, freqs_index: Array[Int], persite_lnl: Array[Double]): Double

  @native private[pll] def pll_compute_edge_loglikelihood(partition: pll_partition_t, parent_clv_index: Int, parent_scaler_index: Int, child_clv_index: Int, child_scaler_index: Int, matrix_index: Int, freqs_index: Array[Int], persite_lnl: Array[Double]): Double

  @native private[pll] def pll_update_sumtable(partition: pll_partition_t, parent_clv_index: Int, child_clv_index: Int, params_indices: Array[Int], sumtable: Pointer): Int

  @native private[pll] def pll_compute_likelihood_derivatives(partition: pll_partition_t, parent_scaler_index: Int, child_scaler_index: Int, branch_length: Double, params_indices: Array[Int], sumtable: Pointer, d_f: Array[Double], dd_f: Array[Double]): Double

  @native private[pll] def pll_compute_gamma_cats(alpha: Double, categories: Int, output_rates: Array[Double]): Int

}
