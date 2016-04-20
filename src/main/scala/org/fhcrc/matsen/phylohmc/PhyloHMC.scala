package org.fhcrc.matsen.phylohmc

import spire.algebra._
import spire.random.{Dist, Gaussian, Generator, Uniform}
import spire.syntax.innerProductSpace._
import spire.syntax.order._

case class Z[R : Field, N](q: Tree[R, N], p: Tree[R, N], Minv: Tree[Tree[R, N], N], L: Tree[Tree[R, N], N])(_U: => (R, Tree[R, N]), _K: => (R, Tree[R, N])) {
  lazy val (u, dU) = _U
  lazy val (k, dK) = _K
  lazy val H = u + k
  def copy(q: Tree[R, N] = this.q, p: Tree[R, N] = this.p, Minv: Tree[Tree[R, N], N] = this.Minv, L: Tree[Tree[R, N], N] = this.L)(_U: => (R, Tree[R, N]) = (u, dU), _K: => (R, Tree[R, N]) = (k, dK)): Z[R, N] =
    Z(q, p, Minv, L)(_U, _K)
  def nni(b: Branch[N], which: Boolean, _U: => (R, Tree[R, N]), _K: => (R, Tree[R, N])): Z[R, N] =
    Z(q.nni(b, which), p.nni(b, which), Minv.mapLengths(_.nni(b, which)).nni(b, which), L.mapLengths(_.nni(b, which)).nni(b, which))(_U, _K)
}

abstract class PhyloHMC[R : Field : NRoot : Trig : Order, N](val U: Tree[R, N] => (R, Tree[R, N]), val eps: R, val alpha: R, val L: Int)(implicit val rng: Generator, implicit val uniformImpl: Uniform[R], implicit val gaussianImpl: Gaussian[R]) extends (Z[R, N] => Z[R, N]) {

  val uniform = Dist.uniform(Field[R].fromInt(0), Field[R].fromInt(1))
  val gaussian = Dist.gaussian[R](Field[R].fromInt(0), Field[R].fromInt(1))
  val sqrtalpha = NRoot[R].sqrt(alpha)
  val sqrt1malpha = NRoot[R].sqrt(1 - alpha)

  def K(p: Tree[R, N], Minv: Tree[Tree[R, N], N]): (R, Tree[R, N])

  def leapfrog(eps: R)(z: Z[R, N]): Z[R, N] = {
    val halfEps = eps / 2
    val pp = z.p - halfEps *: z.dU
    val qp = z.q + eps *: z.dK
    val dUp = U(pp)
    val ppp = qp - halfEps *: dUp._2
    Z(qp, ppp, z.Minv, z.L)(dUp, K(ppp, z.Minv))
  }

  def leapprog(eps: R)(z: Z[R, N]): Z[R, N]

  def solveForEps(br: Branch[N], z: Z[R, N]): R = {
    val twoa = z.Minv(br) dot z.dU
    val b = z.Minv(br) dot z.p
    val c = z.q(br)
    val mbdiv2a = - b / twoa
    val discriminantdiv2a = NRoot[R].sqrt(b * b - 2 * twoa * c) / twoa
    val eps1 = mbdiv2a + discriminantdiv2a
    val eps2 = mbdiv2a - discriminantdiv2a
    if (eps1 >= 0 && eps2 >= 0) eps1 min eps2 else eps1 max eps2
  }

  def flipMomentum(z: Z[R, N]): Z[R, N] = z.copy(p = -z.p)(_K = (z.k, -z.dK))

  def corruptMomentum(z: Z[R, N]): Z[R, N] = {
    val r = z.q.mapLengths(_ => rng.next(gaussian))
    val pp = z.p.branches.foldLeft(z.p) { (p, b) =>
      val x = z.L(b) dot r
      p.updated(b, p(b) * sqrt1malpha + x * sqrtalpha)
    }
    z.copy(p = pp)(_K = K(pp, z.Minv))
  }

  def simulateDynamics(z: Z[R, N]): Z[R, N] = (0 until L).foldLeft(z)((z, _) => leapprog(eps)(z))

  override def apply(z: Z[R, N]): Z[R, N] = {
    val zp = flipMomentum(simulateDynamics(z))
    val a = Trig[R].exp(z.H - zp.H) min 1
    if (rng.next(uniform) < a) zp else z
  }

}

trait ReflectivePhyloHMC[R, N] extends PhyloHMC[R, N]

trait VuPhyloHMC[R, N] extends PhyloHMC[R, N]
