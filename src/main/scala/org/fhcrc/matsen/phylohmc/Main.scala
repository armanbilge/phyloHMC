package org.fhcrc.matsen.phylohmc

import scala.io.Source
import reflect.runtime.currentMirror
import tools.reflect.ToolBox

object Main extends App {

  val toolBox = currentMirror.mkToolBox()
  toolBox.compile(toolBox.parse("import org.fhcrc.matsen.phylohmc._\n" + Source.fromFile(args(0)).getLines().mkString("\n")))()

}
