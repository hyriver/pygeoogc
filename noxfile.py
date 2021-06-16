import nox


@nox.session(python="3.9")
def tests(session):
    session.install(".[test]")
    session.run("pytest")
    session.run("coverage", "report")
    session.run("coverage", "html")
