name := "phyloHMC"
enablePlugins(GitVersioning)
scalaVersion := "2.11.8"
scalacOptions += "-optimize"
libraryDependencies += "com.github.julien-truffaut" %% "monocle-core" % "1.2.0"
libraryDependencies += "io.github.nicolasstucki" %% "multisets" % "0.3"
libraryDependencies += "org.apache.commons" % "commons-math3" % "3.6.1"
libraryDependencies += "org.scala-lang" % "scala-compiler" % "2.11.8"
libraryDependencies += "org.spire-math" %% "spire" % "0.11.0"
