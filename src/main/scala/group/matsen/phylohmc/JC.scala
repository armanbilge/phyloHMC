package group.matsen.phylohmc

import group.matsen.phylohmc.SubstitutionModel._
import spire.algebra.{Field, Trig}
import spire.syntax.field._

class JC[R : Field : Trig] extends SubstitutionModel[R] {

  override val stationaryDistribution: IndexedSeq[R] = IndexedSeq.fill(4)(Field[R].fromDouble(0.25))

  override val Q: Matrix[_4, R] = {
    val qchange = Field[R].one / 3
    val qnochange = Field[R].one
    Matrix[_4, R](qnochange, qchange, qchange, qchange,
      qchange, qnochange, qchange, qchange,
      qchange, qchange, qnochange, qchange,
      qchange, qchange, qchange, qnochange
    )
  }

  override def apply(t: R): Matrix[_4, R] = {
    val pchange = (1 - Trig[R].exp(- 4.0/3.0 * t)) / 4
    val pnochange = 1 - 3 * pchange
    Matrix[_4, R](pnochange, pchange, pchange, pchange,
      pchange, pnochange, pchange, pchange,
      pchange, pchange, pnochange, pchange,
      pchange, pchange, pchange, pnochange
    )
  }

}

object JC {

  def apply[R : Field : Trig] = new JC[R]

}
