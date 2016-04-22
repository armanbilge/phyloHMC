package org.fhcrc.matsen.phylohmc

import spire.algebra.Field
import spire.syntax.field._

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
