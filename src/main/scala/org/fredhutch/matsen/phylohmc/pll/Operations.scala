package org.fredhutch.matsen.phylohmc.pll

class Operations(val count: Int) {

  private[pll] val pll_operation = new pll_operation()
  private[this] val operations = pll_operation.toArray(count).asInstanceOf[Array[pll_operation]]

  def update(index: Int)(parentCLVIndex: Int = operations(index).parent_clv_index, parentScalerIndex: Int = operations(index).parent_scaler_index, child1CLVIndex: Int = operations(index).child1_clv_index, child1MatrixIndex: Int = operations(index).child1_matrix_index, child1ScalerIndex: Int = operations(index).child1_scaler_index, child2CLVIndex: Int = operations(index).child2_clv_index, child2MatrixIndex: Int = operations(index).child2_matrix_index, child2ScalerIndex: Int = operations(index).child2_scaler_index): Unit = {
    val operation = operations(index)
    operation.parent_clv_index = parentCLVIndex
    operation.parent_scaler_index = parentScalerIndex
    operation.child1_clv_index = child1CLVIndex
    operation.child1_matrix_index = child1MatrixIndex
    operation.child1_scaler_index = child1ScalerIndex
    operation.child2_clv_index = child2CLVIndex
    operation.child2_matrix_index = child2MatrixIndex
    operation.child2_scaler_index = child2ScalerIndex
  }

}
