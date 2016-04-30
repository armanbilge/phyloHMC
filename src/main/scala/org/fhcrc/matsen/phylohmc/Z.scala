package org.fhcrc.matsen.phylohmc

import spire.algebra.Field
import spire.syntax.field._

case class Z[R : Field, N](q: Tree[R, N], p: IndexedSeq[R])(_U: => (R, IndexedSeq[R]), _K: => (R, IndexedSeq[R])) {
  lazy val (u, dU) = _U
  lazy val (k, dK) = _K
  lazy val H = u + k
  def copy(q: Tree[R, N] = this.q, p: IndexedSeq[R] = this.p)(_U: => (R, IndexedSeq[R]) = (u, dU), _K: => (R, IndexedSeq[R]) = (k, dK)): Z[R, N] = Z(q, p)(_U, _K)
}

object Z {

  def apply[R : Field, N](q: Tree[R, N], p: IndexedSeq[R], U: Tree[R, N] => (R, IndexedSeq[R]), K: IndexedSeq[R] => (R, IndexedSeq[R])): Z[R, N] = Z(q, p)(U(q), K(p))

}
