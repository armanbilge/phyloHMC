package org.fredhutch.matsen.phylohmc

import spire.algebra.{Field, Trig}
import spire.std.seq._
import spire.syntax.vectorSpace._

class ExponentialBranchPrior[R : Field : Trig, N](val lambda: R) extends (Tree[R, N] => (R, IndexedSeq[R])) {

  val logLambda = Trig[R].log(lambda)

  override def apply(t: Tree[R, N]): (R, IndexedSeq[R]) = (t.lengths.size * logLambda + Field[R].sum(-lambda *: t.lengths), IndexedSeq.fill(t.lengths.size)(-lambda))

}

object ExponentialBranchPrior {

  def apply[R : Field : Trig, N](lambda: R): ExponentialBranchPrior[R, N] = new ExponentialBranchPrior[R, N](lambda)

}
