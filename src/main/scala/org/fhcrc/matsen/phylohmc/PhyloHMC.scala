package org.fhcrc.matsen.phylohmc

import spire.algebra.{Field, Trig}
import spire.syntax.vectorSpace._

case class Z[N, R : Field](q: Tree[R, N], p: Tree[R, N])(_U: => (R, Tree[R, N]), _K: => (R, Tree[R, N])) {
  lazy val (u, dU) = _U
  lazy val (k, dK) = _K
}

class PhyloHMC[N, R : Field : Trig](val U: Tree[R, N] => (R, Tree[R, N]), val K: Tree[R, N] => (R, Tree[R, N]), val eps: R) extends (Z[N, R] => Z[N, R]) {


  override def apply(z: Z[N, R]): Z[N, R] = {

    def leapfrog(eps: R)(z: Z[N, R]): Z[N, R] = {
      val halfEps = eps / 2
      val pp = z.p - halfEps *: z.dU
      val qp = z.q + eps *: z.dK
      val dUp = U(pp)
      val ppp = qp - halfEps *: dUp._2
      Z(qp, ppp)(dUp, K(ppp))
    }

    // TODO
    null

  }

}
