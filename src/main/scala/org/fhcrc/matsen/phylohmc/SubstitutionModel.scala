package org.fhcrc.matsen.phylohmc

trait SubstitutionModel[R] extends (R => Matrix[R]) {

  val stationaryDistribution: IndexedSeq[R]

  def apply(t: R): Matrix[R]

}
