name := "phyloHMC"
enablePlugins(GitVersioning)
scalaVersion := "2.12.1"
scalacOptions += "-optimize"
libraryDependencies += "com.chuusai" %% "shapeless" % "2.3.2"
libraryDependencies += "com.github.julien-truffaut" %% "monocle-core" % "1.4.0-M2"
libraryDependencies += "com.lihaoyi" % "ammonite" % "0.8.1" cross CrossVersion.full
libraryDependencies += "io.github.nicolasstucki" %% "multisets" % "0.4"
libraryDependencies += "net.java.dev.jna" % "jna" % "4.2.2"
libraryDependencies += "org.apache.commons" % "commons-math3" % "3.6.1"
libraryDependencies += "org.scala-lang.modules" %% "scala-parser-combinators" % "1.0.5"
libraryDependencies += "org.scala-lang" % "scala-compiler" % scalaVersion.value
libraryDependencies += "org.spire-math" %% "spire" % "0.13.0"
