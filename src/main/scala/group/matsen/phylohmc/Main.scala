package group.matsen.phylohmc

import java.io.File

import ammonite.ops.Path
import ammonite.util.Res.{Exception, Failure}

object Main extends App {

  ammonite.Main(predef = "import " + getClass.getPackage.getName + "._", verboseOutput = false).runScript(Path(new File(args(0)).getAbsolutePath), Seq.empty, Seq.empty) match {
    case Failure(ex, msg) =>
      ex match {
        case Some(t) => throw t
        case None => throw new RuntimeException(msg)
      }
    case Exception(t, msg) => throw t
    case _ =>
  }

}
