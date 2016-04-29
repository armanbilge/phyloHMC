package org.fhcrc.matsen.phylohmc

import spire.algebra.{Field, Trig}
import spire.syntax.field._

class JC[@specialized(Double) R : Field : Trig] extends SubstitutionModel[R] {

  override val stationaryDistribution: IndexedSeq[R] = IndexedSeq.fill(4)(Field[R].fromDouble(0.25))

  override def apply(t: R): Matrix[R] = {
    val pchange = (1 - Trig[R].exp(- 4.0/3.0 * t)) / 4
    val pnochange = 1 - 3 * pchange
    Matrix(pnochange, pchange, pchange, pchange,
      pchange, pnochange, pchange, pchange,
      pchange, pchange, pnochange, pchange,
      pchange, pchange, pchange, pnochange
    )
  }

}

object JC {

  def apply[R : Field : Trig] = new JC[R]

}
