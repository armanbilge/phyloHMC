package group.matsen.phylohmc

import shapeless.Nat

trait SubstitutionModel[R] extends (R => Matrix[Nat._4, R]) {

  val stationaryDistribution: IndexedSeq[R]
  val Q: Matrix[Nat._4, R]

  def apply(t: R): Matrix[Nat._4, R]

}
