package group.matsen.phylohmc.pll

import java.util.concurrent.Executors

import group.matsen.phylohmc.{Pattern, Patterns, Tree}
import spire.algebra.VectorSpace

import scala.collection.parallel.ExecutionContextTaskSupport
import scala.concurrent.ExecutionContext

class ParallelTreeLikelihood[N](val threads: Int, val patterns: Patterns, val gtr: GTR, val mu: Double, sse: Boolean = false, avx: Boolean = false, avx2: Boolean = false, avx512: Boolean = false, tipPatternCompression: Boolean = false) extends (Tree[Double, N] => (Double, IndexedSeq[Double])) {

  private[this] val likelihoods = patterns.multiplicities.grouped(patterns.multiplicities.length / threads)
    .map(x => Patterns.from(x: _*)(Patterns.configuration.compact[Pattern]))
    .map(x => new TreeLikelihood[N](x, gtr, mu, sse, avx, avx2, avx512, tipPatternCompression)).toParArray
  likelihoods.tasksupport = new ExecutionContextTaskSupport(new ExecutionContext {
    val threadPool = Executors.newFixedThreadPool(threads)
    def execute(runnable: Runnable) = threadPool.submit(runnable)
    def reportFailure(t: Throwable) = throw t
  })

  override def apply(t: Tree[Double, N]): (Double, IndexedSeq[Double]) = {
    import spire.std.double._
    import spire.std.seq._
    val (ls, dLs) = likelihoods.map(l => l(t)).unzip
    (ls.sum, VectorSpace[IndexedSeq[Double], Double].sum(dLs.arrayseq))
  }


}
