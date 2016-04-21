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
  def nni(b: Branch[N], which: Boolean, U: Tree[R, N] => (R, Tree[R, N]), K: Tree[Tree[R, N], N] => Tree[R, N] => (R, Tree[R, N])): Z[R, N] = {
    val qp = q.nni(b, which)
    val pp = p.updated(b, -p(b)).nni(b, which)
    val Minvp = Minv.mapLengths(_.nni(b, which)).nni(b, which)
    Z(qp, pp, Minvp, L.mapLengths(_.nni(b, which)).nni(b, which))(U(qp), K(Minvp)(pp))
  }
}

abstract class PhyloHMC[R : NRoot : Trig : Uniform : Gaussian, N](val U: Tree[R, N] => (R, Tree[R, N]), val eps: R, val alpha: R, val L: Int)(implicit val rng: Generator, implicit val f: Field[R], implicit val s: Signed[R], implicit val o: Order[R]) extends (Z[R, N] => Z[R, N]) {

  val uniform = Dist.uniform(Field[R].zero, Field[R].one)
  val gaussian = Dist.gaussian(Field[R].zero, Field[R].one)
  val sqrtalpha = NRoot[R].sqrt(alpha)
  val sqrt1malpha = NRoot[R].sqrt(1 - alpha)

  def K(Minv: Tree[Tree[R, N], N])(p: Tree[R, N]): (R, Tree[R, N])

  def leapfrog(eps: R)(z: Z[R, N]): Z[R, N] = {
    val halfEps = eps / 2
    val pp = z.p - halfEps *: z.dU
    val qp = z.q + eps *: z.dK
    val dUp = U(pp)
    val ppp = qp - halfEps *: dUp._2
    Z(qp, ppp, z.Minv, z.L)(dUp, K(ppp, z.Minv))
  }

  def leapprog(eps: R)(z: Z[R, N]): Z[R, N]

  def solveForEps(z: Z[R, N])(br: Branch[N]): R = {
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
    z.copy(p = pp)(_K = K(z.Minv)(pp))
  }

  def simulateDynamics(z: Z[R, N]): Z[R, N] = (0 until L).foldLeft(z)((z, _) => leapprog(eps)(z))

  override def apply(z: Z[R, N]): Z[R, N] = {
    val zp = flipMomentum(simulateDynamics(z))
    val a = Trig[R].exp(z.H - zp.H) min 1
    if (rng.next(uniform) < a) zp else z
  }

}

trait ReflectivePhyloHMC[R, N] extends PhyloHMC[R, N] {

  def leapprog(eps: R)(z: Z[R, N]): Z[R, N] = {
    z.q.branches.map(b => (b, solveForEps(z)(b))).filter(x => x._2 <= eps).toList.sortWith(_._2 < _._2).view.map(_._1).foldLeft({
      val zp = leapfrog(eps)(z)
      zp.copy(q = zp.q.mapLengths(Signed[R].abs))()
    }) { (z, b) =>
      rng.nextInt(3) match {
        case 0 => z
        case 1 => z.nni(b, false, U, K)
        case 2 => z.nni(b, true, U, K)
      }
    }
  }

}

trait VuPhyloHMC[R, N] extends PhyloHMC[R, N] {

  def leapprog(eps: R)(z: Z[R, N]): Z[R, N] = {
    val (zp, e) = z.q.branches.map(b => (b, solveForEps(z)(b))).filter(x => x._2 <= eps).toList.sortWith(_._2 < _._2).foldLeft((z, Field[R].zero)) { (ze, bt) =>
      val (z, e) = ze
      val (b, t) = bt
      val zp = leapfrog(t - e)(z)
      (rng.nextInt(3) match {
        case 0 => zp
        case 1 => zp.nni(b, false, U, K)
        case 2 => zp.nni(b, true, U, K)
      }, t)
    }
    leapfrog(eps - e)(zp)
  }

}
