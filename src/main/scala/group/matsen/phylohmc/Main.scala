package group.matsen.phylohmc

import java.io.File

import ammonite.ops.Path

object Main extends App {

  ammonite.Main(predef = "import " + getClass.getPackage.getName + "._", verboseOutput = false).runScript(Path(new File(args(0)).getAbsolutePath), Seq.empty, Seq.empty)

}
