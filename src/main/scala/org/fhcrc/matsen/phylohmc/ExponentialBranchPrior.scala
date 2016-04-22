package org.fhcrc.matsen.phylohmc

import spire.algebra.{Field, Trig}
import spire.syntax.vectorSpace._

class ExponentialBranchPrior[R : Field : Trig, N](val lambda: R) extends (Tree[R, N] => (R, Tree[R, N])) {

  val logLambda = Trig[R].log(lambda)

  override def apply(t: Tree[R, N]): (R, Tree[R, N]) = {
    val tp = -lambda *: t
    (t.branches.size * logLambda + Field[R].sum(tp.lengths.values), tp)
  }

}

object ExponentialBranchPrior {

  def apply[R : Field : Trig, N](lambda: R): ExponentialBranchPrior[R, N] = new ExponentialBranchPrior[R, N](lambda)

}