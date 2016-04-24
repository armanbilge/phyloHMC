package org.fhcrc.matsen.phylohmc

import spire.algebra.Signed

trait ReflectiveLeapFrog[R, N] extends PhyloHMC[R, N] {

  def leapprog(eps: R)(z: Z[R, N]): Z[R, N] = {
    val zp = leapfrog(eps)(z)
    zp.copy(q = zp.q.mapLengths(Signed[R].abs(_)))()
  }

}
