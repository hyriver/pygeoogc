import nox


@nox.session(python="3.9")
def tests(session):
    hr_deps = ["async_retriever"]
    for p in hr_deps:
        session.install(f"git+https://github.com/cheginit/{p}.git")
    session.install(".[test]")
    session.run("pytest")
    session.run("coverage", "report")
    session.run("coverage", "html")
