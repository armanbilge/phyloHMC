package org.fhcrc.matsen.phylohmc

trait SubstitutionModel[R] extends (R => Matrix[R]) {

  val stationaryDistribution: Vector[R]

  def apply(t: R): Matrix[R]

}
