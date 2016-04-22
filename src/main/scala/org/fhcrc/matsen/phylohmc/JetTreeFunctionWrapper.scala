package org.fhcrc.matsen.phylohmc

import spire.algebra._
import spire.math.{Jet, JetDim}

import scala.reflect.ClassTag

class JetTreeFunctionWrapper[R : Eq : Field : NRoot : IsReal : Trig : ClassTag, N](f: Tree[Jet[R], N] => Jet[R]) extends (Tree[R, N] => (R, Tree[R, N])) {

  override def apply(t: Tree[R, N]): (R, Tree[R, N]) = {
    implicit val jd = JetDim(t.branches.size)
    val branch2int = t.branches.zipWithIndex.toMap
    val y = f(t.mapLengths((b, l) => Jet(l, branch2int(b))))
    (y.real, t.mapLengths((b, _) => y.infinitesimal(branch2int(b))))
  }

}
