package group.matsen.phylohmc

import scala.io.Source
import scala.reflect.runtime.currentMirror
import scala.tools.reflect.ToolBox

object Main extends App {

  val toolBox = currentMirror.mkToolBox()
  toolBox.compile(toolBox.parse("import " + getClass.getPackage.getName + "._\n" + Source.fromFile(args(0)).getLines().mkString("\n")))()

}
