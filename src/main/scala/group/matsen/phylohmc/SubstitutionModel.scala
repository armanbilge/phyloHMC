package group.matsen.phylohmc

import group.matsen.phylohmc.SubstitutionModel._4
import shapeless.Witness

trait SubstitutionModel[R] extends (R => Matrix[_4, R]) {

  val stationaryDistribution: IndexedSeq[R]
  val Q: Matrix[_4, R]

  def apply(t: R): Matrix[_4, R]

}

object SubstitutionModel {
  implicit val w4 = Witness(4)
  type _4 = w4.T
}
