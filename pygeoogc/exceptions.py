"""Customized PyGeoOGC exceptions."""
from typing import Optional

import async_retriever as ar


class ServiceError(ar.ServiceError):
    """Exception raised when the requested data is not available on the server.

    Parameters
    ----------
    err : str
        Service error message.
    """


class InvalidInputValue(ar.InvalidInputValue):
    """Exception raised for invalid input.

    Parameters
    ----------
    inp : str
        Name of the input parameter
    valid_inputs : tuple
        List of valid inputs
    """


class InvalidInputType(ar.InvalidInputType):
    """Exception raised when a function argument type is invalid.

    Parameters
    ----------
    arg : str
        Name of the function argument
    valid_type : str
        The valid type of the argument
    example : str, optional
        An example of a valid form of the argument, defaults to None.
    """


class MissingInputs(ValueError):
    """Exception raised when there are missing function arguments."""


class ServiceUnavailable(Exception):
    """Exception raised when the service is not available.

    Parameters
    ----------
    url : str
        The server url
    """

    def __init__(self, url: str) -> None:
        self.message = f"Service is currently not available, try again later:\n{url}"
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class ZeroMatched(ValueError):
    """Exception raised when a function argument is missing.

    Parameters
    ----------
    msg : str
        The exception error message
    """

    def __init__(self, msg: Optional[str] = None) -> None:
        if msg is None:
            self.message = "Service returned no features."
        else:
            self.message = f"Service returned no features with the following error message:\n{msg}"
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message
