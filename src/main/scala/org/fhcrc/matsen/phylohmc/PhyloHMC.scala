package org.fhcrc.matsen.phylohmc

import spire.algebra._
import spire.random.{Dist, Gaussian, Generator, Uniform}
import spire.syntax.innerProductSpace._
import spire.syntax.order._

abstract class PhyloHMC[R : NRoot : Trig : Uniform : Gaussian, N](val posterior: Tree[R, N] => (R, Tree[R, N]), val eps: R, val alpha: R, val L: Int)(implicit val rng: Generator, implicit val f: Field[R], implicit val s: Signed[R], implicit val o: Order[R]) extends (Z[R, N] => Z[R, N]) {

  val uniform = Dist.uniform(Field[R].zero, Field[R].one)
  val gaussian = Dist.gaussian(Field[R].zero, Field[R].one)
  val sqrtalpha = NRoot[R].sqrt(alpha)
  val sqrt1malpha = NRoot[R].sqrt(1 - alpha)

  def U(q: Tree[R, N]): (R, Tree[R, N]) = {
    val (p, dP) = posterior(q)
    (-p, -dP)
  }

  def K(Minv: Tree[Tree[R, N], N])(p: Tree[R, N]): (R, Tree[R, N]) = {
    val Minvp = p.mapLengths((b, _) => Minv(b) dot p)
    ((p dot Minvp) / 2, Minvp)
  }

  def leapfrog(eps: R)(z: Z[R, N]): Z[R, N] = {
    val halfEps = eps / 2
    val pp = z.p - halfEps *: z.dU
    val qp = z.q + eps *: z.dK
    val Up = U(qp)
    val ppp = qp - halfEps *: Up._2
    Z(qp, ppp, z.Minv, z.L)(Up, K(z.Minv)(ppp))
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
    val pp = z.p.mapLengths((b, l) => l * sqrt1malpha + (z.L(b) dot r) * sqrtalpha)
    z.copy(p = pp)(_K = K(z.Minv)(pp))
  }

  def simulateDynamics(z: Z[R, N]): Z[R, N] = (0 until L).foldLeft(z)((z, _) => leapprog(eps)(z))

  override def apply(z: Z[R, N]): Z[R, N] = {
    val zp = flipMomentum(simulateDynamics(z))
    val a = Trig[R].exp(z.H - zp.H) min 1
    corruptMomentum(flipMomentum(if (rng.next(uniform) < a) zp else z))
  }

}
