package org.fredhutch.matsen.phylohmc

import java.util

import com.sun.jna.{Native, Structure}

package object pll {

  Native.register("pll")

  private[pll] trait PartitionPointer extends Structure.ByReference

  private[pll] class OperationPointer(parent_clv_index: Int, parent_scaler_index: Int, child1_clv_index: Int, child1_matrix_index: Int, child1_scaler_index: Int, child2_clv_index: Int, child2_matrix_index: Int, child2_scaler_index: Int) extends Structure with Structure.ByReference {
    override def getFieldOrder: util.List[_] = util.Arrays.asList("parent_clv_index", "parent_scaler_index", "child1_clv_index", "child1_matrix_index", "child1_scaler_index", "child2_clv_index", "child2_matrix_index", "child2_scaler_index")
  }

  @native private[pll] def pll_partition_create(tips: Int, clv_buffers: Int, states: Int, sites: Int, rate_matrices: Int, prob_matrices: Int, rate_cats: Int, scale_buffers: Int, attributes: Int): PartitionPointer

  @native private[pll] def pll_partition_destroy(partition: PartitionPointer): Unit

  @native private[pll] def pll_set_tip_states(partition: PartitionPointer, tip_index: Int, map: Array[Int], sequence: String): Int

  @native private[pll] def pll_set_tip_clv(partition: PartitionPointer, tip_index: Int, clv: Array[Double]): Unit

  @native private[pll] def pll_set_pattern_weights(partition: PartitionPointer, pattern_weights: Array[Int]): Unit

  @native private[pll] def pll_set_subst_params(partition: PartitionPointer, params_index: Int, params: Array[Double]): Unit

  @native private[pll] def pll_set_frequencies(partition: PartitionPointer, params_index: Int, frequencies: Array[Double]): Unit

  @native private[pll] def pll_update_prob_matrices(partition: PartitionPointer, params_index: Int, matrix_indices: Array[Int], branch_lengths: Array[Double], count: Int): Unit

  @native private[pll] def pll_update_invariant_sites(partition: PartitionPointer): Int

  @native private[pll] def pll_update_invariant_sites_proportion(partition: PartitionPointer, params_index: Int, prop_invar: Double): Int

  @native private[pll] def pll_update_partials(partition: PartitionPointer, operations: OperationPointer, count: Int): Double

  @native private[pll] def pll_compute_root_loglikelihood(partition: PartitionPointer, clv_index: Int, scaler_index: Int, freqs_index: Array[Int], persite_lnl: Array[Double]): Double

  @native private[pll] def pll_compute_edge_loglikelihood(partition: PartitionPointer, parent_clv_index: Int, parent_scaler_index: Int, child_clv_index: Int, child_scaler_index: Int, matrix_index: Int, freqs_index: Array[Int], persite_lnl: Array[Double]): Double

  @native private[pll] def pll_update_sumtable(partition: PartitionPointer, parent_clv_index: Int, child_clv_index: Int, params_indices: Array[Int], sumtable: Array[Double]): Int

  @native private[pll] def pll_compute_likelihood_derivatives(partition: PartitionPointer, parent_scaler_index: Int, child_scaler_index: Int, branch_length: Double, params_indices: Array[Int], sumtable: Array[Double], d_f: Array[Double], dd_f: Array[Double]): Double

}
