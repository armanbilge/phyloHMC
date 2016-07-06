package org.fredhutch.matsen.phylohmc.pll

import org.fredhutch.matsen.phylohmc.{Branch, Taxon, Tree}

class TreeLikelihood[N](patterns: Map[Map[Taxon, Char], Int], gtr: GTR, mu: Double, sse: Boolean = false, avx: Boolean = false, avx2: Boolean = false, avx512: Boolean = false, tipPatternCompression: Boolean = false) extends (Tree[Double, N] => (Double, IndexedSeq[Double])) {

  val taxa = patterns.keys.head.keySet
  val taxaToInt = taxa.zipWithIndex.toMap

  val leafCount = taxa.size
  val internalNodeCount = taxa.size - 2
  val nodeCount = leafCount + internalNodeCount
  val branchCount = nodeCount - 1
  val patternCount = patterns.size

  private[this] val partition = {
    val partition = new Partition(leafCount, internalNodeCount, 4, patternCount, 1, branchCount, 1, 0, sse = sse, avx = avx, avx2 = avx2, avx512 = avx512, tipPatternCompression = tipPatternCompression)
    val patternSeq = patterns.toSeq
    taxa.foreach { t =>
      partition.setTipStates(taxaToInt(t), Nt, patternSeq.map(_._1(t)).mkString)
    }
    partition.setPatternWeights(patternSeq.map(_._2).toArray)
    partition.setCategoryRates(Array(mu))
    partition.setCategoryWeights(Array(1.0))
    partition.setFrequencies(0, gtr.pi.toArray)
    partition.setSubstParams(0, gtr.rates.toArray)
    partition.updateEigen(0)
    partition
  }

  private[this] val operations = new Operations(internalNodeCount)

  override def apply(t: Tree[Double, N]): (Double, IndexedSeq[Double]) = {

    val branchCount = t.branches.size

    partition.updateProbMatrices(Array(0), (0 until branchCount).toArray, t.lengths.toArray, branchCount)

    def recurse(parent: N, child: N, i: Int = 0, m: Map[N, Int] = Map()): (Int, Int, Map[N, Int]) = if (t.isLeaf(child)) {
      val clv = taxaToInt(t.taxa(child))
      (i, clv, m + (child -> clv))
    } else {
      val children = t.children(child, parent)
      val (left, right) = (children.head, children.tail.head)
      val (j, clv1, mp) = recurse(child, left, i, m)
      val (k, clv2, mpp) = recurse(child, right, j, mp)
      val clv = k + leafCount
      operations.update(k)(parentCLVIndex = clv, child1CLVIndex = clv1, child2CLVIndex = clv2, child1MatrixIndex = t.branchesToIndex(Branch(child, left)), child2MatrixIndex = t.branchesToIndex(Branch(child, right)))
      (k + 1, clv, mpp + (child -> clv))
    }

    val rho = t.nodes.filter(t.isLeaf).head
    val root = t.neighbors(rho).head
    val (_, rootCLV, m) = recurse(rho, root)
    val rhoCLV = taxaToInt(t.taxa(rho))
    val nodesToCLVs = m + (rho -> rhoCLV)
    partition.updatePartials(operations)
    val logL = partition.computeEdgeLogLikelihood(rhoCLV, -1, rootCLV, -1, t.branchesToIndex(Branch(rho, root)), Array(0), null)
    val df = new Array[Double](1)
    val ddf = new Array[Double](1)
    val dfs = new Array[Double](branchCount)
    t.branches.foreach { b =>
      val i = t.branchesToIndex(b)
      partition.updateSumtable(nodesToCLVs(b.head), nodesToCLVs(b.tail), Array(0))
      partition.computeLikelihoodDerivatives(-1, -1, t.lengths(i), Array(0), df, ddf)
      dfs(i) = df(0)
    }
    (logL, dfs.toIndexedSeq)
  }

}
