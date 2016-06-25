package org.fredhutch.matsen.phylohmc

import org.apache.commons.math3.linear.{Array2DRowRealMatrix, CholeskyDecomposition, LUDecomposition}
import spire.algebra._
import spire.random.{Dist, Gaussian, Generator, Uniform}
import spire.std.seq._
import spire.syntax.innerProductSpace._
import spire.syntax.order._

abstract class PhyloHMC[R : Trig : Uniform : Gaussian, N](val posterior: Tree[R, N] => (R, IndexedSeq[R]), val M: Matrix[R], val alpha: R, val eps: R, val L: Int, val RToDouble: R => Double)(implicit val rng: Generator, implicit val f: Field[R], implicit val n: NRoot[R], implicit val s: Signed[R], implicit val o: Order[R]) extends (Z[R, N] => Z[R, N]) {

  val (invM, choleskyL) = {
    val apacheM = new Array2DRowRealMatrix(M.size, M.size)
    M.indices.foreach(Function.tupled((i, j) => apacheM.setEntry(i, j, RToDouble(M(i, j)))))
    val apacheInvM = new LUDecomposition(apacheM).getSolver.getInverse
    val apacheCholeskyL = new CholeskyDecomposition(apacheM).getL
    (Matrix(M.size)((i, j) => Field[R].fromDouble(apacheInvM.getEntry(i, j))), Matrix(M.size)((i, j) => Field[R].fromDouble(apacheCholeskyL.getEntry(i, j))))
  }

  val uniform = Dist.uniform(Field[R].zero, Field[R].one)
  val gaussian = Dist.gaussian(Field[R].zero, Field[R].one)
  val sqrtalpha = NRoot[R].sqrt(alpha)
  val sqrt1malpha = NRoot[R].sqrt(1 - alpha)

  def U(q: Tree[R, N]): (R, IndexedSeq[R]) = {
    val (p, dP) = posterior(q)
    (-p, -dP)
  }

  def K(p: IndexedSeq[R]): (R, IndexedSeq[R]) = {
    val invMp = invM * p
    ((p dot invMp) / 2, invMp)
  }

  def leapprog(eps: R)(z: Z[R, N]): Z[R, N]

  def flipMomentum(z: Z[R, N]): Z[R, N] = z.copy(p = -z.p)(_K = (z.k, -z.dK))

  def corruptMomentum(z: Z[R, N]): Z[R, N] = {
    val r = IndexedSeq.fill(z.p.size)(rng.next(gaussian))
    val pp = sqrt1malpha *: z.p + sqrtalpha *: (choleskyL * r)
    z.copy(p = pp)(_K = K(pp))
  }

  def simulateDynamics(z: Z[R, N]): Z[R, N] = (0 until L).foldLeft(z)((z, _) => leapprog(eps)(z))

  override def apply(z: Z[R, N]): Z[R, N] = {
    val zp = flipMomentum(simulateDynamics(z))
    val a = Trig[R].exp(z.H - zp.H) min 1
    corruptMomentum(flipMomentum(if (rng.next(uniform) < a) zp else z))
  }

}
