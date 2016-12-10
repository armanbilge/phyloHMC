package group.matsen.phylohmc.pll

import group.matsen.phylohmc.{Branch, Patterns, Tree}

import scala.collection.mutable

class TreeLikelihood[N](val patterns: Patterns, val gtr: GTR, val mu: Double, sse: Boolean = false, avx: Boolean = false, avx2: Boolean = false, avx512: Boolean = false, tipPatternCompression: Boolean = false) extends (Tree[Double, N] => (Double, IndexedSeq[Double])) {

  val taxa = patterns.head.keySet
  val taxaToInt = taxa.zipWithIndex.toMap

  val leafCount = taxa.size
  val internalNodeCount = taxa.size - 2
  val nodeCount = leafCount + internalNodeCount
  val branchCount = nodeCount - 1
  val patternCount = patterns.multiplicities.size

  private[this] val partition = {
    val partition = new Partition(leafCount, 3 * internalNodeCount, 4, patternCount, 1, branchCount, 1, 0, sse = sse, avx = avx, avx2 = avx2, avx512 = avx512, tipPatternCompression = tipPatternCompression)

    partition.setPatternWeights(patterns.multiplicities.map(_._2).toArray)
    taxa.foreach { taxon =>
      val i = taxaToInt(taxon)
      partition.setTipStates(i, Nt, patterns.multiplicities.map(_._1).map(_(taxon)).map(Vector('A', 'C', 'G', 'T')).mkString)
    }
    partition.setCategoryRates(Array(mu))
    partition.setCategoryWeights(Array(1.0))
    partition.setFrequencies(0, gtr.pi.toArray)
    partition.setSubstParams(0, gtr.rates.toArray)
    partition.updateEigen(0)
    partition
  }

  private[this] val operations = new Operations(3 * internalNodeCount)

  override def apply(t: Tree[Double, N]): (Double, IndexedSeq[Double]) = {

    val branchCount = t.branches.size

    partition.updateProbMatrices(Array(0), (0 until branchCount).toArray, t.lengths.toArray, branchCount)

    val indices = mutable.Map[(N, N), Int]()

    var i = leafCount
    t.nodes.filter(t.isLeaf).foreach { n =>
      indices((n, t.neighbors(n).head)) = taxaToInt(t.taxa(n))
    }

    def postorder(n: N, p: N): Unit = {
      if (t.isInternal(n)) {
        t.children(n, p).foreach(postorder(_, n))
        val c = t.children(n, p)
        val l = c.head
        val r = c.tail.head
        operations.update(i - leafCount)(parentCLVIndex = i, child1CLVIndex = indices((l, n)), child1MatrixIndex = t.branchesToIndex(Branch(l, n)), child2CLVIndex = indices((r, n)), child2MatrixIndex = t.branchesToIndex(Branch(r, n)))
        indices((n, p)) = i
        i += 1
      }
    }

    def preorder(n: N, p: N): Unit = {
      if (t.isInternal(p)) {
        val c = t.children(p, n)
        val l = c.head
        val r = c.tail.head
        operations.update(i - leafCount)(parentCLVIndex = i, child1CLVIndex = indices((l, p)), child1MatrixIndex = t.branchesToIndex(Branch(l, p)), child2CLVIndex = indices((r, p)), child2MatrixIndex = t.branchesToIndex(Branch(r, p)))
        indices((p, n)) = i
        i += 1
      }
      t.children(n, p).foreach(preorder(_, n))
    }

    val rootEdge = t.branches.head
    postorder(rootEdge.head, rootEdge.tail)
    postorder(rootEdge.tail, rootEdge.head)
    t.children(rootEdge.head, rootEdge.tail).foreach(preorder(_, rootEdge.head))
    t.children(rootEdge.tail, rootEdge.head).foreach(preorder(_, rootEdge.tail))

    partition.updatePartials(operations)
    val b = t.branches.head
    val logL = partition.computeEdgeLogLikelihood(indices((b.head, b.tail)), -1, indices((b.tail, b.head)), -1, t.branchesToIndex(b), Array(0), null)
    val df = new Array[Double](1)
    val ddf = new Array[Double](1)
    val dfs = new Array[Double](branchCount)
    t.branches.foreach { b =>
      val i = t.branchesToIndex(b)
      partition.updateSumtable(indices((b.head, b.tail)), indices((b.tail, b.head)), Array(0))
      partition.computeLikelihoodDerivatives(-1, -1, t.lengths(i), Array(0), df, ddf)
      dfs(i) = -df(0)
    }
    (logL, dfs.toIndexedSeq)
  }

}
